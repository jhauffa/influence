package edu.tum.cs.graph.clique;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Deque;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

import com.carrotsearch.hppc.ShortHashSet;
import com.carrotsearch.hppc.ShortSet;
import com.carrotsearch.hppc.cursors.ShortCursor;

import edu.tum.cs.util.ExperimentConfiguration;

/**
 * A parallelized adaptation of algorithm 1 of Reid et al., 2012. Same limitations as CliquePercolation.
 */
public class CliquePercolationParallel extends CliquePercolationBase<Queue<short[]>> {

	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(CliquePercolationParallel.class);

	public CliquePercolationParallel(int minSize, boolean useLifoOrder) {
		super(minSize, useLifoOrder);
	}

	@Override
	protected Queue<short[]> createCliqueContainer() {
		return new ConcurrentLinkedQueue<short[]>();
	}

	@Override
	protected Queue<short[]> createFrontier() {
		if (useLifoOrder)
			return new LinkedBlockingDeque<short[]>();
		return new PriorityBlockingQueue<short[]>(11, new CliqueSizeComparator());
	}

	private class PercolationWorker extends Thread {
		private final AtomicInteger workCount;
		public final ShortSet curComm = new ShortHashSet();

		public PercolationWorker(AtomicInteger workCount) {
			this.workCount = workCount;
		}

		@Override
		public void run() {
			while (!interrupted()) {
				short[] cliqueJ;
				try {
					cliqueJ = ((BlockingQueue<short[]>) frontier).take();
				} catch (InterruptedException ex) {
					break;
				}

				Iterator<short[]> it = cliques.iterator();
				while (it.hasNext()) {
					short[] cliqueK = it.next();
					if (isOverlapping(cliqueJ, cliqueK, minSize - 1)) {
						it.remove();

						// This code has a race condition by design: cliqueK may be successfully tested for overlap by
						// multiple threads before being removed from the list, and will thus be added to the frontier
						// more than once. This turns out to be rather memory inefficient when using many threads, so
						// we work around this by keeping a hash set of the cliques that were most recently added to the
						// frontier.
						workCount.getAndIncrement();
						if (addCliqueToFrontierCached(cliqueK)) {
							for (short m : cliqueK)
								curComm.add(m);
						} else
							workCount.decrementAndGet();
					}
				}

				if (workCount.decrementAndGet() == 0) {
					synchronized (workCount) {
						workCount.notifyAll();
					}
				}
			}
		}
	}

	private static final int numCacheEntries = 1000;
	protected Set<short[]> frontierCache = new HashSet<short[]>(numCacheEntries);
	protected Deque<short[]> frontierCacheOrder = new LinkedList<short[]>();

	private boolean addCliqueToFrontierCached(short[] clique) {
		synchronized (frontierCache) {
			if (frontierCache.contains(clique))
				return false;
			frontierCache.add(clique);
			frontierCacheOrder.add(clique);
			if (frontierCache.size() >= numCacheEntries)
				frontierCache.remove(frontierCacheOrder.poll());
		}
		addCliqueToFrontier(clique);
		return true;
	}

	@Override
	protected void printProgress() {
		logger.info("remaining cliques = " + cliques.size() + ", frontier size = " + frontier.size());
	}

	@Override
	protected Collection<short[]> findCommunities() {
		Collection<short[]> communities = new ArrayList<short[]>();
		AtomicInteger workCount = new AtomicInteger();
		int numThreads = cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_THREADS,
				Runtime.getRuntime().availableProcessors());
		PercolationWorker[] workers = new PercolationWorker[numThreads];
		for (int j = 0; j < workers.length; j++) {
			workers[j] = new PercolationWorker(workCount);
			workers[j].start();
		}

		while (!cliques.isEmpty()) {
			short[] cliqueI = cliques.remove();
			ShortSet curComm = new ShortHashSet();
			for (short m : cliqueI)
				curComm.add(m);
			workCount.set(2);
			addCliqueToFrontier(cliqueI);

			try {
				synchronized (workCount) {
					if (workCount.decrementAndGet() > 0)
						workCount.wait();
				}
			} catch (InterruptedException ex) {
				throw new RuntimeException("interrupted while waiting for worker threads", ex);
			}

			for (PercolationWorker worker : workers) {
				for (ShortCursor csr : worker.curComm)
					curComm.add(csr.value);
				worker.curComm.clear();
			}
			communities.add(curComm.toArray());
		}

		for (PercolationWorker worker : workers) {
			worker.interrupt();
			try {
				worker.join();
			} catch (InterruptedException ex) {
				throw new RuntimeException("interrupted while waiting for worker threads", ex);
			}
		}
		return communities;
	}

}
