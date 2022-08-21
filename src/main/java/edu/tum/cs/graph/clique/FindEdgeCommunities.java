package edu.tum.cs.graph.clique;

import java.io.File;
import java.util.Collection;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Logger;

import edu.tum.cs.graph.GraphCreator;
import edu.tum.cs.util.ExperimentConfiguration;
import edu.tum.cs.util.io.Serializer;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.UndirectedGraph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;

public class FindEdgeCommunities {

	private static final Logger logger = Logger.getLogger(FindEdgeCommunities.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(FindEdgeCommunities.class);

	private static class PartitionDensityTask implements Runnable {
		private static final int reportInterval = 100;

		private final EdgeClustering<Long, Integer> clust;
		private final AtomicInteger workCount;
		private final double[] score;
		private final int start, end;

		public PartitionDensityTask(EdgeClustering<Long, Integer> clust, AtomicInteger workCount, double[] score,
				int start, int end) {
			this.clust = clust;
			this.workCount = workCount;
			this.score = score;
			this.start = start;
			this.end = end;
		}

		@Override
		public void run() {
			for (int i = start; i < end; i++) {
				score[i] = clust.getPartitionDensity(i);

				int c = workCount.incrementAndGet();
				if ((c % reportInterval) == 0)
					logger.info("processed " + c + " levels");
			}
		}
	}

	public static void main(String[] args) throws Exception {
		boolean requireComm = false;
		boolean requireBidir = false;
		String pathPrefix = ".";
		for (String arg : args) {
			if (arg.equals("comm"))
				requireComm = true;
			else if (arg.equals("bidir"))
				requireBidir = true;
			else
				pathPrefix = arg;
		}

		logger.info("generating social network graph");
		DirectedGraph<Long, Integer> graph = GraphCreator.createGraph(requireComm);
		UndirectedGraph<Long, Integer> undirectedGraph = new UndirectedSparseGraph<Long, Integer>();
		for (Long vertex : graph.getVertices())
			undirectedGraph.addVertex(vertex);
		int idx = 0;
		for (Integer edge : graph.getEdges()) {
			Pair<Long> endpoints = graph.getEndpoints(edge);
			if (endpoints.getFirst().equals(endpoints.getSecond())) {
				logger.severe("vertex " + endpoints.getFirst() + " has self edge, ignoring");
				continue;
			}
			if (!requireBidir || (graph.findEdge(endpoints.getSecond(), endpoints.getFirst()) != null))
				undirectedGraph.addEdge(idx++, endpoints);
		}
		graph = null;

		logger.info("extracting communities (" + undirectedGraph.getEdgeCount() + " undirected edges)");
		long startTime = System.currentTimeMillis();

		EdgeClustering<Long, Integer> clust = new EdgeClustering<Long, Integer>(undirectedGraph, true);
		clust.process();

		int maxLevel = clust.getMaxLevel();
		logger.info("maximum level = " + maxLevel);
		AtomicInteger workCount = new AtomicInteger(0);
		double[] score = new double[maxLevel + 1];
		int numThreads = cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_THREADS,
				Runtime.getRuntime().availableProcessors());
		ExecutorService pool = Executors.newFixedThreadPool(numThreads);
		Future<?>[] tasks = new Future[numThreads];
		int start = 0;
		int sliceSize = (maxLevel + 1) / numThreads;
		for (int i = 0; i < (numThreads - 1); i++) {
			int end = start + sliceSize;
			tasks[i] = pool.submit(new PartitionDensityTask(clust, workCount, score, start, end));
			start = end;
		}
		tasks[numThreads - 1] = pool.submit(new PartitionDensityTask(clust, workCount, score, start, maxLevel + 1));
		for (Future<?> task : tasks)
			task.get();
		pool.shutdown();

		int optLevel = 0;
		double optScore = 0.0;
		for (int i = 0; i < score.length; i++) {
			System.out.println(i + ";" + score[i]);
			if (score[i] > optScore) {
				optScore = score[i];
				optLevel = i;
			}
		}

		long endTime = System.currentTimeMillis();
		Collection<Collection<Long>> communities = clust.getVertexClusters(optLevel);
		logger.info("found " + communities.size() + " communities in " + (endTime - startTime) + " milliseconds");

		File f = new File(pathPrefix, "communities.ser");
		logger.info("writing communities to '" + f.getPath() + "'...");
		Serializer.saveObjectToFile(communities, f);
		FindPercolatedCommunities.saveCommunitiesAsText(communities, new File(pathPrefix, "communities-onmi.txt"));
		FindPercolatedCommunities.printCoverage(communities, undirectedGraph.getVertexCount());
		logger.info("done");
	}

}
