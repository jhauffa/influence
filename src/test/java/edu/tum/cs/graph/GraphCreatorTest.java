package edu.tum.cs.graph;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.Random;

import org.junit.Test;

import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;

public class GraphCreatorTest {

	private static final int numVertices = 10000;
	private static final double p = 0.01;
	private static long seed = 34877L;

	@Test
	public void testComplementGraph() {
		// construct random graph
		Random rng = new Random(seed);
		DirectedSparseGraph<Long, Integer> graph = new DirectedSparseGraph<Long, Integer>();
		for (Long i = 0L; i < numVertices; i++)
			graph.addVertex(i);
		Integer e = 0;
		for (Long i = 0L; i < numVertices; i++) {
			for (Long j = 0L; j < numVertices; j++) {
				if (rng.nextDouble() < p)
					graph.addEdge(e++, i, j);
			}
		}

		DirectedSparseGraph<Long, Integer> compGraph = GraphCreator.sampleComplementGraph(graph, graph.getEdgeCount());
		assertEquals(graph.getVertexCount(), compGraph.getVertexCount());
		assertEquals(graph.getEdgeCount(), compGraph.getEdgeCount());
		for (Integer ce : compGraph.getEdges()) {
			Pair<Long> v = compGraph.getEndpoints(ce);
			assertNull(graph.findEdge(v.getFirst(), v.getSecond()));
		}
	}

}
