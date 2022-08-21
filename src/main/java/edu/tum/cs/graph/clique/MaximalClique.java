package edu.tum.cs.graph.clique;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Phaser;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

import edu.tum.cs.util.ExperimentConfiguration;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedGraph;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * Implementation of version 2 of the Bron-Kerbosch algorithm (with pivoting) as described in Bron and Kerbosch, 1971.
 * @param <V> vertex class of graph
 * @param <E> edge class of graph
 */
public class MaximalClique<V, E> {

	private static final Logger logger = Logger.getLogger(MaximalClique.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(MaximalClique.class);
	private static final AtomicLong seenMaximalCliques = new AtomicLong(0);

	private final VertexIndex<V> vertexIndex;
	protected final int minSize;
	protected final int[][] adjacency;
	private SubgraphOutput<V> lastOutput;

	private static class WriterThread<V> extends Thread {
		private final BlockingQueue<Collection<V>> workQueue = new ArrayBlockingQueue<Collection<V>>(20000);
		private final File file;

		public WriterThread(File file) {
			this.file = file;
		}

		@Override
		public void run() {
			try {
				Writer writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(
						new FileOutputStream(file))), 32 * 1024 * 1024);
				try {
					while (!interrupted()) {
						Collection<V> subgraph;
						try {
							subgraph = workQueue.take();
						} catch (InterruptedException ex) {
							break;
						}

						for (V v : subgraph)
							writer.write(v.toString() + ";");
						writer.write("\n");
					}
				} finally {
					writer.close();
				}
			} catch (IOException ex) {
				logger.log(Level.SEVERE, "error closing output file", ex);
			}
		}

		public void enqueueSubgraph(Collection<V> subgraph) {
			try {
				workQueue.put(subgraph);
			} catch (InterruptedException ex) {
				logger.log(Level.SEVERE, "interrupted while waiting for writer thread", ex);
			}
		}
	}

	protected static class SubgraphOutput<V> {
		private static final int reportInterval = 10000000;

		private final Collection<short[]> subgraphs;
		private final VertexIndex<V> vertexIndex;
		private final Map<V, VertexSubgraphStats<V>> largestSubgraphForNodes = new HashMap<V, VertexSubgraphStats<V>>();
		private final WriterThread<V> writerThread;
		private long numSubgraphsFound = 0;

		private SubgraphOutput(VertexIndex<V> vertexIndex, Collection<short[]> subgraphs,
				WriterThread<V> writerThread) {
			this.vertexIndex = vertexIndex;
			this.subgraphs = subgraphs;
			this.writerThread = writerThread;
		}

		public SubgraphOutput(VertexIndex<V> vertexIndex, Collection<short[]> subgraphs) {
			this(vertexIndex, subgraphs, null);
		}

		public SubgraphOutput(VertexIndex<V> vertexIndex, File file) {
			this(vertexIndex, null, new WriterThread<V>(file));
			writerThread.start();
		}

		public synchronized void writeSubgraph(Collection<Integer> subgraph) {
			Collection<V> unpackedSubgraph = vertexIndex.unpackVertices(subgraph);
			if (writerThread != null)
				writerThread.enqueueSubgraph(unpackedSubgraph);
			else {
				short[] packedSubgraph = new short[subgraph.size()];
				int idx = 0;
				for (Integer v : subgraph)
					packedSubgraph[idx++] = v.shortValue();
				subgraphs.add(packedSubgraph);
			}

			for (V v : unpackedSubgraph) {
				VertexSubgraphStats<V> stat = largestSubgraphForNodes.get(v);
				if (stat == null) {
					stat = new VertexSubgraphStats<V>(1, unpackedSubgraph);
					largestSubgraphForNodes.put(v, stat);
				} else {
					stat.numSubgraphs++;
					if (subgraph.size() > stat.maxSizeSubgraph.size())
						stat.maxSizeSubgraph = unpackedSubgraph;
				}
			}

			// report progress
			if ((++numSubgraphsFound % reportInterval) == 0)
				logger.info(numSubgraphsFound + " subgraphs found");
		}

		public void flush() {
			if (writerThread != null)
				writerThread.interrupt();
		}
	}

	private class FindCliquesTask extends RecursiveAction {
		private static final long serialVersionUID = 1L;

		private final Phaser sync;
		private final SubgraphOutput<V> output;
		private final LinkedList<Integer> potentialClique;
		private final int[] subset;
		private int notEndIdx;
		private final int candEndIdx;

		public FindCliquesTask(Phaser sync, SubgraphOutput<V> output) {
			this(sync, output, new LinkedList<Integer>(), new int[adjacency.length], 0, adjacency.length);
			for (int i = 0; i < adjacency.length; i++)
				subset[i] = i;
		}

		private FindCliquesTask(Phaser sync, SubgraphOutput<V> output, LinkedList<Integer> potentialClique,
				int[] subset, int notEndIdx, int candEndIdx) {
			this.sync = sync;
			this.output = output;
			this.potentialClique = potentialClique;
			this.subset = subset;
			this.notEndIdx = notEndIdx;
			this.candEndIdx = candEndIdx;

			sync.register();
		}

		@Override
		protected void compute() {
			// Find vertex in union of candidates and excluded vertices (== array "subset") with maximum number of
			// connections to other candidates (== minimum number of missing edges), store index in refIdx.
			int pos = -1;
			int refVertex = -1;
			int candIdx = -1;
			int minNumMissing = candEndIdx;
			int missingCounter = 0;
			for (int i = 0; i < candEndIdx; i++) {
				if (minNumMissing == 0)
					break;
				int n = 0;
				for (int j = notEndIdx; j < candEndIdx; j++) {
					if (n >= minNumMissing)
						break;
					if (adjacency[subset[i]][subset[j]] == 0) {
						n++;
						pos = j;
					}
				}
				if (n < minNumMissing) {	// found new minimum
					refVertex = subset[i];
					minNumMissing = n;
					// select first candidate
					if (i < notEndIdx) {
						candIdx = pos;
					} else {
						candIdx = i;
						missingCounter = 1;	// reference vertex in candidates -> pre-increment counter by 1
					}
				}
			}

			for (missingCounter += minNumMissing; missingCounter >= 1; missingCounter--) {
				// move candidate to the beginning of the candidate list
				int candVertex = subset[candIdx];
				subset[candIdx] = subset[notEndIdx];
				subset[notEndIdx] = candVertex;

				int[] newSubset = new int[candEndIdx];

				// build new set of excluded vertices
				int newNotEndIdx = 0;
				for (int i = 0; i < notEndIdx; i++)
					if (adjacency[candVertex][subset[i]] != 0)
						newSubset[newNotEndIdx++] = subset[i];

				// build new set of candidates
				int newCandEndIdx = newNotEndIdx;
				for (int i = notEndIdx + 1; i < candEndIdx; i++)	// exclude current candidate
					if (adjacency[candVertex][subset[i]] != 0)
						newSubset[newCandEndIdx++] = subset[i];

				potentialClique.addLast(candVertex);
				if (newCandEndIdx == 0) {	// implies newNotEndIdx == 0 -> found maximal clique
					seenMaximalCliques.incrementAndGet();
					if (potentialClique.size() >= minSize)
						output.writeSubgraph(potentialClique);
				} else if (newNotEndIdx < newCandEndIdx) {
					// recursive call
					Phaser nextLevelSync = new Phaser(sync);
					FindCliquesTask recursiveTask = new FindCliquesTask(nextLevelSync, output,
							new LinkedList<Integer>(potentialClique), newSubset, newNotEndIdx, newCandEndIdx);
					recursiveTask.fork();
				}
				potentialClique.removeLast();
				notEndIdx++;	// add candidate to set of excluded vertices

				if (missingCounter > 1) {
					// choose next candidate
					candIdx = notEndIdx;
					while (adjacency[refVertex][subset[candIdx]] != 0)
						candIdx++;
				}
			}

			sync.arrive();
		}
	}

	protected MaximalClique(Graph<V, E> graph, VertexIndex<V> vertexIndex, boolean requireBidirectionalEdge,
			int minSize) {
		this.vertexIndex = vertexIndex;
		this.minSize = minSize;

		int n = graph.getVertexCount();
		adjacency = new int[n][n];
		for (E edge : graph.getEdges()) {
			Pair<V> endpoints = graph.getEndpoints(edge);
			if (!requireBidirectionalEdge || (graph.findEdge(endpoints.getSecond(), endpoints.getFirst()) != null)) {
				short idx1 = vertexIndex.get(endpoints.getFirst());
				short idx2 = vertexIndex.get(endpoints.getSecond());
				adjacency[idx1][idx2] = 1;
				adjacency[idx2][idx1] = 1;
			}
		}
		for (int i = 0; i < n; i++)
			adjacency[i][i] = 1;	// required by this implementation of BK
	}

	public MaximalClique(UndirectedGraph<V, E> graph, VertexIndex<V> vertexIndex) {
		this((Graph<V, E>) graph, vertexIndex, false, 2);
	}

	public MaximalClique(DirectedGraph<V, E> dgraph, VertexIndex<V> vertexIndex, boolean requireBidirectionalEdge) {
		this((Graph<V, E>) dgraph, vertexIndex, requireBidirectionalEdge, 2);
	}

	public MaximalClique(DirectedGraph<V, E> dgraph, VertexIndex<V> vertexIndex, boolean requireBidirectionalEdge,
			int minSize) {
		this((Graph<V, E>) dgraph, vertexIndex, requireBidirectionalEdge, minSize);
	}

	protected RecursiveAction createInitialTask(Phaser sync, SubgraphOutput<V> output) {
		return new FindCliquesTask(sync, output);
	}

	private void computeMaximalCliques(SubgraphOutput<V> output) {
		ForkJoinPool pool = new ForkJoinPool(cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_THREADS,
				Runtime.getRuntime().availableProcessors()));
		Phaser sync = new Phaser();
		sync.register();
		pool.invoke(createInitialTask(sync, output));
		sync.arriveAndAwaitAdvance();
		pool.shutdown();
		try {
			if (!pool.awaitTermination(1, TimeUnit.SECONDS))
				logger.warning("thread pool did not shut down orderly");
		} catch (InterruptedException ex) {
			throw new RuntimeException("interrupted while waiting for thread pool to shut down", ex);
		}
		output.flush();
		lastOutput = output;
		logger.info("seen " + seenMaximalCliques.get() + " maximal cliques");
	}

	public ArrayList<short[]> getAllMaximalSubgraphs() {
		ArrayList<short[]> cliques = new ArrayList<short[]>();
		computeMaximalCliques(new SubgraphOutput<V>(vertexIndex, cliques));
		return cliques;
	}

	public void saveAllMaximalSubgraphs(File file) {
		computeMaximalCliques(new SubgraphOutput<V>(vertexIndex, file));
	}

	public Map<V, VertexSubgraphStats<V>> getLargestMaximalSubgraphForNodes() {
		return lastOutput.largestSubgraphForNodes;
	}

}
