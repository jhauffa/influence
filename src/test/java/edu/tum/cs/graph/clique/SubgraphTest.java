package edu.tum.cs.graph.clique;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

import edu.uci.ics.jung.graph.UndirectedGraph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;

public class SubgraphTest {

	public static <V> UndirectedGraph<V, Integer> buildGraph(V[][] edgeList) {
		UndirectedGraph<V, Integer> graph = new UndirectedSparseGraph<>();
		int edgeIdx = 0;
		for (V[] edge : edgeList) {
			graph.addVertex(edge[0]);
			graph.addVertex(edge[1]);
			graph.addEdge(edgeIdx++, edge[0], edge[1]);
		}
		return graph;
	}

	public static void verifySubgraphs(Collection<Collection<Long>> subgraphs, Long[][] refSubgraphs) {
		Set<Set<Long>> subgraphsSet = new HashSet<Set<Long>>();
		for (Collection<Long> subgraph : subgraphs)
			subgraphsSet.add(new HashSet<Long>(subgraph));

		assertEquals(subgraphsSet.size(), subgraphs.size());
		for (Long[] refSubgraph : refSubgraphs)
			assertTrue(subgraphsSet.contains(new HashSet<Long>(Arrays.asList(refSubgraph))));
	}

	@Test
	public void testMaximalClique() {
		// 2 \     / 6
		// |  3 - 4  |
		// 1 /     \ 5
		Long[][] edgeList = new Long[][] {
			{ 1L, 2L }, { 1L, 3L }, { 2L, 3L }, { 3L, 4L }, { 4L, 5L }, { 4L, 6L }, { 5L, 6L }
		};
		Long[][] refCliques = new Long[][] {
			{ 1L, 2L, 3L }, { 3L, 4L }, { 4L, 5L, 6L }
		};

		UndirectedGraph<Long, Integer> graph = buildGraph(edgeList);
		VertexIndex<Long> index = new VertexIndex<Long>(graph);
		MaximalClique<Long, Integer> mc = new MaximalClique<Long, Integer>(graph, index);
		verifySubgraphs(index.unpackAllVertices(mc.getAllMaximalSubgraphs()), refCliques);

		// 1   2
		// | X |
		// 3 - 4
		edgeList = new Long[][] {
			{ 1L, 3L }, { 3L, 4L }, { 3L, 2L }, { 4L, 1L }, { 4L, 2L }
		};
		refCliques = new Long[][] {
			{ 1L, 3L, 4L }, { 2L, 3L, 4L }
		};

		graph = buildGraph(edgeList);
		index = new VertexIndex<Long>(graph);
		mc = new MaximalClique<Long, Integer>(graph, index);
		verifySubgraphs(index.unpackAllVertices(mc.getAllMaximalSubgraphs()), refCliques);

		graph.addEdge(5, 1L, 2L);
		mc = new MaximalClique<Long, Integer>(graph, index);
		verifySubgraphs(index.unpackAllVertices(mc.getAllMaximalSubgraphs()), new Long[][] { { 1L, 2L, 3L, 4L } });
	}

	@Test
	public void testMaximalPlex() {
		// 2 \     / 6
		//    3 - 4
		// 1 /     \ 5
		Long[][] edgeList = new Long[][] {
			{ 1L, 3L }, { 2L, 3L }, { 3L, 4L }, { 4L, 5L }, { 4L, 6L }
		};
		Long[][] refPlexes = new Long[][] {
			{ 1L, 2L, 3L }, { 1L, 3L, 4L }, { 2L, 3L, 4L }, { 3L, 4L, 5L }, { 3L, 4L, 6L }, { 4L, 5L, 6L }
		};

		UndirectedGraph<Long, Integer> graph = buildGraph(edgeList);
		VertexIndex<Long> index = new VertexIndex<Long>(graph);
		MaximalPlex<Long, Integer> mp = new MaximalPlex<Long, Integer>(graph, index, 1);
		verifySubgraphs(index.unpackAllVertices(mp.getAllMaximalSubgraphs()), refPlexes);

		edgeList = new Long[][] {
			{ 1L, 2L }, { 1L, 3L }, { 1L, 4L }, { 1L, 6L }, { 2L, 4L }, { 2L, 5L }, { 2L, 6L }, { 2L, 7L }, { 2L, 8L },
			{ 3L, 4L }, { 3L, 5L }, { 3L, 6L }, { 4L, 5L }, { 5L, 6L }, { 5L, 7L }, { 5L, 8L }, { 7L, 8L }, { 7L, 9L },
			{ 8L, 9L }
		};
		refPlexes = new Long[][] {
			{ 2L, 4L, 5L, 8L }, { 2L, 4L, 5L, 7L }, { 1L, 2L, 3L, 4L, 5L, 6L }, { 2L, 5L, 6L, 8L }, { 2L, 5L, 6L, 7L },
			{ 2L, 5L, 7L, 8L }, { 1L, 2L, 7L }, { 1L, 2L, 8L }, { 2L, 7L, 8L, 9L }, { 5L, 7L, 8L, 9L }, { 3L, 5L, 7L },
			{ 3L, 5L, 8L }
		};

		graph = buildGraph(edgeList);
		index = new VertexIndex<Long>(graph);
		mp = new MaximalPlex<Long, Integer>(graph, index, 1);
		verifySubgraphs(index.unpackAllVertices(mp.getAllMaximalSubgraphs()), refPlexes);
	}

	@Test
	public void testMultipleMissingEdges() {
		// 1   2
		// | X
		// 3 - 4
		Long[][] edgeList = new Long[][] {
			{ 1L, 3L }, { 3L, 4L }, { 3L, 2L }, { 4L, 1L }
		};
		Long[][] refPlexes1 = new Long[][] {
			{ 1L, 2L, 3L }, { 1L, 3L, 4L }, { 2L, 3L, 4L }
		};
		Long[][] refPlexes2 = new Long[][] {
			{ 1L, 2L, 3L, 4L }
		};

		UndirectedGraph<Long, Integer> graph = buildGraph(edgeList);
		VertexIndex<Long> index = new VertexIndex<Long>(graph);
		MaximalPlex<Long, Integer> mp = new MaximalPlex<Long, Integer>(graph, index, 1);
		verifySubgraphs(index.unpackAllVertices(mp.getAllMaximalSubgraphs()), refPlexes1);
		mp = new MaximalPlex<Long, Integer>(graph, index, 2);
		verifySubgraphs(index.unpackAllVertices(mp.getAllMaximalSubgraphs()), refPlexes2);
	}

}
