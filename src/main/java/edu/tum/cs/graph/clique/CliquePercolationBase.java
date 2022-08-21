package edu.tum.cs.graph.clique;

import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.Deque;
import java.util.Queue;
import java.util.logging.Logger;

/**
 * Abstract base class to facilitate code sharing between CliquePercolation and CliquePercolationParallel.
 */
public abstract class CliquePercolationBase<T extends Collection<short[]>> {

	protected static final Logger logger = Logger.getLogger(CliquePercolationBase.class.getName());
	protected static final int reportIntervalSeconds = 10;

	protected final int minSize;
	protected final boolean useLifoOrder;
	protected final T cliques;
	protected final Queue<short[]> frontier;

	protected static class CliqueSizeComparator implements Comparator<short[]> {
		/** not consistent with #equals, only for use with PriorityQueue */
		@Override
		public int compare(short[] o1, short[] o2) {
			return o2.length - o1.length;
		}
	}

	private class ProgressMonitor extends Thread {
		@Override
		public void run() {
			while (!interrupted()) {
				printProgress();
				try {
					Thread.sleep(reportIntervalSeconds * 1000);
				} catch (InterruptedException ex) {
					break;
				}
			}
		}
	}

	public CliquePercolationBase(int minSize, boolean useLifoOrder) {
		this.minSize = minSize;
		this.useLifoOrder = useLifoOrder;
		this.cliques = createCliqueContainer();
		this.frontier = createFrontier();
	}

	protected abstract T createCliqueContainer();
	protected abstract Queue<short[]> createFrontier();

	public void addClique(short[] clique) {
		if (clique.length >= minSize) {
			Arrays.sort(clique);	// prepare for binary search
			cliques.add(clique);
		}
	}

	public void addAllCliques(Iterable<short[]> cliques) {
		for (short[] clique : cliques)
			addClique(clique);
	}

	protected static boolean isOverlapping(short[] c1, short[] c2, int minOverlap) {
		if (c1.length > c2.length) {
			short[] tmp = c1;
			c1 = c2;
			c2 = tmp;
		}

		int idx = 0;
		while ((minOverlap > 0) && ((idx + minOverlap - 1) < c1.length)) {
			if (Arrays.binarySearch(c2, c1[idx]) >= 0)
				minOverlap--;
			idx++;
		}
		return (minOverlap == 0);
	}

	protected void addCliqueToFrontier(short[] clique) {
		if (useLifoOrder)
			((Deque<short[]>) frontier).offerFirst(clique);
		else
			frontier.offer(clique);
	}

	protected abstract Collection<short[]> findCommunities();
	protected abstract void printProgress();

	public Collection<short[]> getCommunities() {
		Thread progressMonitor = new ProgressMonitor();
		progressMonitor.start();
		Collection<short[]> communities = findCommunities();
		progressMonitor.interrupt();
		try {
			progressMonitor.join();
		} catch (InterruptedException ex) {
			throw new RuntimeException("interrupted while waiting for progress monitor thread", ex);
		}
		return communities;
	}

}
