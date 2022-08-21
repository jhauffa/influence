package edu.tum.cs.graph.clique;

import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Phaser;
import java.util.concurrent.RecursiveAction;

import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedGraph;

public class MaximalPlex<V, E> extends MaximalClique<V, E> {

	private final int missingEdges;

	private class FindPlexesTask extends RecursiveAction {
		private static final long serialVersionUID = 1L;

		private final Phaser sync;
		private final SubgraphOutput<V> output;
		private final LinkedList<Integer> potentialPlex;
		private final int[] subset;
		private int notEndIdx;
		private final int candEndIdx;

		public FindPlexesTask(Phaser sync, SubgraphOutput<V> output) {
			this(sync, output, new LinkedList<Integer>(), new int[adjacency.length], 0, adjacency.length);
			for (int i = 0; i < adjacency.length; i++)
				subset[i] = i;
		}

		private FindPlexesTask(Phaser sync, SubgraphOutput<V> output, LinkedList<Integer> potentialPlex, int[] subset,
				int notEndIdx, int candEndIdx) {
			this.sync = sync;
			this.output = output;
			this.potentialPlex = potentialPlex;
			this.subset = subset;
			this.notEndIdx = notEndIdx;
			this.candEndIdx = candEndIdx;

			sync.register();
		}

		@Override
		protected void compute() {
			int candIdx = notEndIdx;
			while ((candIdx < candEndIdx) && (notEndIdx < candEndIdx)) {
				// move candidate to the beginning of the candidate list
				int candVertex = subset[candIdx];
				subset[candIdx] = subset[notEndIdx];
				subset[notEndIdx] = candVertex;

				potentialPlex.addLast(candVertex);
				int[] newSubset = new int[candEndIdx];

				// build new set of excluded vertices
				int newNotEndIdx = 0;
				for (int i = 0; i < notEndIdx; i++)
					if (isPlex(potentialPlex, subset[i]))
						newSubset[newNotEndIdx++] = subset[i];

				// build new set of candidates
				int newCandEndIdx = newNotEndIdx;
				for (int i = notEndIdx + 1; i < candEndIdx; i++)	// exclude current candidate
					if (isPlex(potentialPlex, subset[i]))
						newSubset[newCandEndIdx++] = subset[i];

				if (newCandEndIdx == 0) {	// implies newNotEndIdx == 0 -> found maximal plex
					if (potentialPlex.size() >= minSize)
						output.writeSubgraph(potentialPlex);
				} else if (newNotEndIdx < newCandEndIdx) {
					// recursive call
					Phaser nextLevelSync = new Phaser(sync);
					FindPlexesTask recursiveTask = new FindPlexesTask(nextLevelSync, output,
							new LinkedList<Integer>(potentialPlex), newSubset, newNotEndIdx, newCandEndIdx);
					recursiveTask.fork();
				}
				potentialPlex.removeLast();
				notEndIdx++;	// add candidate to set of excluded vertices

				// choose next candidate
				candIdx = notEndIdx;
			}

			sync.arrive();
		}

		private boolean isPlex(List<Integer> potentialPlex, int candidate) {
			int edgesRequired = potentialPlex.size() - missingEdges;
			if (edgesRequired <= 0)
				return true;

			int candidateEdgesRequired = edgesRequired;
			for (int v1 : potentialPlex) {
				if (adjacency[v1][candidate] == 0) {
					// no edge to candidate -> verify that v1 has sufficient edges to other vertices in potentialPlex
					int innerEdgesRequired = edgesRequired;
					for (int v2 : potentialPlex) {
						if ((v1 != v2) && (adjacency[v1][v2] != 0))
							if (--innerEdgesRequired == 0)
								break;
					}
					if (innerEdgesRequired > 0)
						return false;
				} else
					candidateEdgesRequired--;
			}
			return (candidateEdgesRequired <= 0);
		}
	}

	public MaximalPlex(UndirectedGraph<V, E> graph, VertexIndex<V> vertexIndex, int missingEdges) {
		super((Graph<V, E>) graph, vertexIndex, false, 3);
		this.missingEdges = missingEdges;
	}

	public MaximalPlex(DirectedGraph<V, E> dgraph, VertexIndex<V> vertexIndex, int missingEdges,
			boolean requireBidirectionalEdge) {
		super((Graph<V, E>) dgraph, vertexIndex, requireBidirectionalEdge, 3);
		this.missingEdges = missingEdges;
	}

	public MaximalPlex(DirectedGraph<V, E> dgraph, VertexIndex<V> vertexIndex, int missingEdges,
			boolean requireBidirectionalEdge, int minSize) {
		super((Graph<V, E>) dgraph, vertexIndex, requireBidirectionalEdge, minSize);
		this.missingEdges = missingEdges;
	}

	@Override
	protected RecursiveAction createInitialTask(Phaser sync, SubgraphOutput<V> output) {
		return new FindPlexesTask(sync, output);
	}

}
