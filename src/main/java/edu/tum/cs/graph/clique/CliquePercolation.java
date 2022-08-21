package edu.tum.cs.graph.clique;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;

import com.carrotsearch.hppc.ShortHashSet;
import com.carrotsearch.hppc.ShortSet;

/**
 * An adaptation of algorithm 1 of Reid et al., 2012 with reduced memory requirements. The current implementation
 * cannot handle more than Integer.MAX_VALUE cliques. Cliques have to fit into memory; if this should prove to be
 * overly restrictive, the algorithm "COSpoc" of Gregori et al., 2011, could serve as the basis of an external memory
 * algorithm, at the expense of having to perform more overlap tests.
 */
public class CliquePercolation extends CliquePercolationBase<List<short[]>> {

	private int nextFinishedCliqueIdx;

	public CliquePercolation(int minSize, boolean useLifoOrder) {
		super(minSize, useLifoOrder);
	}

	@Override
	protected List<short[]> createCliqueContainer() {
		return new ArrayList<short[]>();	// need random access
	}

	@Override
	protected Queue<short[]> createFrontier() {
		if (useLifoOrder)
			return new LinkedList<short[]>();
		return new PriorityQueue<short[]>(11, new CliqueSizeComparator());
	}

	@Override
	protected void printProgress() {
		logger.info("remaining cliques = " + (nextFinishedCliqueIdx + 1) + ", frontier size = " + frontier.size());
	}

	// TODO: clique list should be shrunk in regular intervals
	private int moveClique(int idx, short[] clique, int nextFinishedCliqueIdx) {
		cliques.set(idx, cliques.get(nextFinishedCliqueIdx));
		cliques.set(nextFinishedCliqueIdx, clique);
		return --nextFinishedCliqueIdx;
	}

	@Override
	protected Collection<short[]> findCommunities() {
		Collection<short[]> communities = new ArrayList<short[]>();

		nextFinishedCliqueIdx = cliques.size() - 1;
		while (nextFinishedCliqueIdx >= 0) {
			short[] cliqueI = cliques.get(0);
			nextFinishedCliqueIdx = moveClique(0, cliqueI, nextFinishedCliqueIdx);

			ShortSet curComm = new ShortHashSet();
			for (short m : cliqueI)
				curComm.add(m);
			addCliqueToFrontier(cliqueI);

			while (!frontier.isEmpty()) {
				short[] cliqueJ = frontier.poll();

				int k = 0;
				while (k <= nextFinishedCliqueIdx) {
					short[] cliqueK = cliques.get(k);
					if (isOverlapping(cliqueJ, cliqueK, minSize - 1)) {
						nextFinishedCliqueIdx = moveClique(k, cliqueK, nextFinishedCliqueIdx);
						for (short m : cliqueK)
							curComm.add(m);
						addCliqueToFrontier(cliqueK);
					} else
						k++;
				}
			}

			communities.add(curComm.toArray());
		}

		return communities;
	}

}
