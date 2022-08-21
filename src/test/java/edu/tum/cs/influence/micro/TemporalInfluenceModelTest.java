package edu.tum.cs.influence.micro;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

public class TemporalInfluenceModelTest {

	private static class TestCase {
		public final double[] past, present, direction;
		public final double magnitude;

		public TestCase(double[] past, double[] present, double[] direction, double magnitude) {
			this.past = past;
			this.present = present;
			this.direction = direction;
			this.magnitude = magnitude;
		}
	}

	private static final TestCase[] testCases = new TestCase[] {
		new TestCase(new double[] { 1.0, 0.0, 0.0, 0.0, 0.0 }, new double[] { 0.0, 0.0, 0.0, 0.0, 1.0 },
				new double[] { 0.0, 0.0, 0.0, 0.0, 1.0 }, 1.0),
		new TestCase(new double[] { 1.0, 0.0, 0.0, 0.0, 0.0 }, new double[] { 1.0, 0.0, 0.0, 0.0, 0.0 },
				new double[] { 1.0 / 5, 1.0 / 5, 1.0 / 5, 1.0 / 5, 1.0 / 5 }, 0.0),
		new TestCase(new double[] { 1.0, 0.0, 0.0, 0.0, 0.0 }, new double[] { 0.5, 0.0, 0.0, 0.0, 0.5 },
				new double[] { 0.0, 0.0, 0.0, 0.0, 1.0 }, 0.5),
		new TestCase(new double[] { 1.0 / 5, 1.0 / 5, 1.0 / 5, 1.0 / 5, 1.0 / 5 },
				new double[] { 1.0, 0.0, 0.0, 0.0, 0.0 }, new double[] { 1.0, 0.0, 0.0, 0.0, 0.0 }, 1.0)
	};
	private static final double eps = 1E-14;

	@Test
	public void testChangeMeasurement() {
		for (TestCase c : testCases) {
			ChangeMeasurement m = ChangeMeasurement.computeChange(c.past, c.present);
			assertArrayEquals(c.direction, m.direction, eps);
			assertEquals(c.magnitude, m.magnitude, eps);
		}
	}

	private static final int numDim = 10;
	private static final int numTimeSteps = 20;

	@Test
	public void testModels() {
		// Create synthetic data: four nodes (A, B, C, D) connected by two edges (A -> B, C -> D), with A having
		// complete influence in the Granger sense on B (B can be fully predicted from A), and C having no influence on
		// D.
		TimeSeriesData[] data = new TimeSeriesData[4];
		for (int i = 0; i < data.length; i++)
			data[i] = new TimeSeriesData(numTimeSteps);
		double[] v2 = new double[numDim];
		double[] v3 = new double[numDim];
		int numActive = (numDim / 2);
		double vActive = 1.0 / numActive;
		for (int i = 0; i < numDim; i++) {
			if (i < numActive)
				v2[i] = vActive;
			else
				v3[i] = vActive;
		}
		for (int i = 0; i < numTimeSteps; i++) {
			double[] v = new double[numDim];
			v[i % numDim] = 1.0;
			data[0].addValue(v);	// A: one-hot, different active component at each time step
			if (i == 0) {
				v = new double[numDim];
				Arrays.fill(v, 1.0 / numDim);
				data[1].addValue(v);
			} else
				data[1].addValue(data[0].getValue(i - 1).clone());	// B: time-lagged copy of A
			data[2].addValue(v2);	// C: uniform distribution over first half of components
			data[3].addValue(v3);	// D: uniform distribution over second half of components
		}
		Map<Long, TimeSeriesData> actions = new HashMap<Long, TimeSeriesData>();
		for (long i = 0L; i < data.length; i++)
			actions.put(i, new TimeSeriesData(data[(int) i]));

		Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> networks =
				new HashMap<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>>();
		DirectedSparseGraph<Long, Integer> network = new DirectedSparseGraph<Long, Integer>();
		network.addEdge(0, 0L, 1L, EdgeType.DIRECTED);
		network.addEdge(1, 2L, 3L, EdgeType.DIRECTED);
		networks.put(ModelVariant.NetworkType.FOLLOWING, network);
		// fake network is only used for evaluation, can be empty
		Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> fakeNetworks =
				new HashMap<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>>();
		fakeNetworks.put(ModelVariant.NetworkType.FOLLOWING, new DirectedSparseGraph<Long, Integer>());

		double[] topicPriorParam = new double[numDim];
		Arrays.fill(topicPriorParam, 1.0 / 3);

		for (ModelVariant.ModelType type : ModelVariant.ModelType.values()) {
			TemporalInfluenceModel<?> model = TemporalInfluenceModel.createModel(type, actions, topicPriorParam);
			Map<ModelVariant, WeightedDirectedGraph<Long, Double>> influenceGraphs =
					model.buildInfluenceGraphs(networks, fakeNetworks, null);

			for (ModelVariant.InfluenceType influence : ModelVariant.InfluenceType.values()) {
				ModelVariant variant = new ModelVariant(type, influence, ModelVariant.NetworkType.FOLLOWING);
				WeightedDirectedGraph<Long, Double> graph = influenceGraphs.get(variant);
				if ((graph == null) && (influence == ModelVariant.InfluenceType.EDGE))
					continue;
				assertNotNull(variant.toString(), graph);
				assertEquals(variant.toString(), 1, graph.getEdgeCount());
				assertEquals(variant.toString(), 2, graph.getVertexCount());
				WeightedEdge<Double> e = graph.getEdges().iterator().next();
				assertTrue(variant.toString() + ", weight = " + e.weight, e.weight > 0.5);
				Pair<Long> v = graph.getEndpoints(e);
				assertEquals(variant.toString(), 0L, (long) v.getFirst());
				assertEquals(variant.toString(), 1L, (long) v.getSecond());
			}
		}
	}

}
