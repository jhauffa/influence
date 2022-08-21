package edu.tum.cs.graph.clique;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collection;

import org.junit.Test;

import edu.tum.cs.math.cluster.SingleLinkageClusteringTest;
import edu.uci.ics.jung.graph.UndirectedGraph;

public class EdgeClusteringTest {

	// Zachary's Karate Club network from http://konect.cc/networks/ucidata-zachary/
	private static final Long[][] karateClubEdgeList = {
		{ 1L, 2L }, { 1L, 3L }, { 2L, 3L }, { 1L, 4L }, { 2L, 4L }, { 3L, 4L }, { 1L, 5L }, { 1L, 6L }, { 1L, 7L },
		{ 5L, 7L }, { 6L, 7L }, { 1L, 8L }, { 2L, 8L }, { 3L, 8L }, { 4L, 8L }, { 1L, 9L }, { 3L, 9L }, { 3L, 10L },
		{ 1L, 11L }, { 5L, 11L }, { 6L, 11L }, { 1L, 12L }, { 1L, 13L }, { 4L, 13L }, { 1L, 14L }, { 2L, 14L },
		{ 3L, 14L }, { 4L, 14L }, { 6L, 17L }, { 7L, 17L }, { 1L, 18L }, { 2L, 18L }, { 1L, 20L }, { 2L, 20L },
		{ 1L, 22L }, { 2L, 22L }, { 24L, 26L }, { 25L, 26L }, { 3L, 28L }, { 24L, 28L }, { 25L, 28L }, { 3L, 29L },
		{ 24L, 30L }, { 27L, 30L }, { 2L, 31L }, { 9L, 31L }, { 1L, 32L }, { 25L, 32L }, { 26L, 32L }, { 29L, 32L },
		{ 3L, 33L }, { 9L, 33L }, { 15L, 33L }, { 16L, 33L }, { 19L, 33L }, { 21L, 33L }, { 23L, 33L }, { 24L, 33L },
		{ 30L, 33L }, { 31L, 33L }, { 32L, 33L }, { 9L, 34L }, { 10L, 34L }, { 14L, 34L }, { 15L, 34L }, { 16L, 34L },
		{ 19L, 34L }, { 20L, 34L }, { 21L, 34L }, { 23L, 34L }, { 24L, 34L }, { 27L, 34L }, { 28L, 34L }, { 29L, 34L },
		{ 30L, 34L }, { 31L, 34L }, { 32L, 34L }, { 33L, 34L }
	};
	private static final Long[][] refClusters = {
		{ 1L, 2L, 18L, 3L, 4L, 20L, 22L, 8L, 13L, 14L }, { 34L, 16L, 19L, 33L, 21L, 23L, 15L }, { 34L, 20L, 14L },
		{ 24L, 26L }, { 3L, 28L }, { 24L, 28L }, { 34L, 33L, 24L, 27L, 28L, 30L }, { 25L, 28L }, { 34L, 33L, 9L, 31L },
		{ 3L, 10L, 29L }, { 24L, 30L }, { 27L, 30L }, { 2L, 31L }, { 1L, 32L }, { 34L, 32L, 33L, 10L, 29L }, { 1L, 9L },
		{ 32L, 25L, 26L }, { 32L, 29L }, { 1L, 5L, 6L, 7L, 11L }, { 33L, 3L, 9L }, { 1L, 12L }, { 17L, 5L, 6L, 7L, 11L }
	};

	@Test
	public void testEdgeClustering() {
		UndirectedGraph<Long, Integer> graph = SubgraphTest.buildGraph(karateClubEdgeList);
		EdgeClustering<Long, Integer> clust = new EdgeClustering<Long, Integer>(graph, false);
		clust.process();
		SingleLinkageClusteringTest.verifyClustering(clust, graph.getEdges());

		int optLevel = -1;
		double maxScore = 0.0;
		for (int i = 0; i <= clust.getMaxLevel(); i++) {
			double score = clust.getPartitionDensity(i);
			if (score > maxScore) {
				maxScore = score;
				optLevel = i;
			}
		}
		assertTrue(optLevel >= 0);

		Collection<Collection<Long>> optClusters = clust.getVertexClusters(optLevel);
		assertEquals(refClusters.length, optClusters.size());
		SubgraphTest.verifySubgraphs(optClusters, refClusters);
	}

}
