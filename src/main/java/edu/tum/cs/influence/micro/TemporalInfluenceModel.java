package edu.tum.cs.influence.micro;

import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;

public abstract class TemporalInfluenceModel<T extends TimeSeriesData> {

	public abstract Map<ModelVariant, WeightedDirectedGraph<Long, Double>> buildInfluenceGraphs(
			Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> networks,
			Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> fakeNetworks,
			Map<ModelVariant, WeightedDirectedGraph<Long, Explanation>> explanations);

	public abstract void writeStatistics(PrintWriter w);

	public abstract Map<String, long[]> reportProgress();

	// TODO: This hack should eventually be replaced with a unified mechanism to retrieve per-variant statistics.
	public final Map<ModelVariant, Integer> numEarlyRejections = new HashMap<ModelVariant, Integer>();

	protected static final double baseSignificanceLevel = 0.05;

	protected static DirectedSparseGraph<Long, Integer> buildUnionNetwork(
			Iterable<DirectedSparseGraph<Long, Integer>> networks) {
		DirectedSparseGraph<Long, Integer> unionNetwork = new DirectedSparseGraph<Long, Integer>();
		int idx = 0;
		for (DirectedSparseGraph<Long, Integer> network : networks) {
			for (Integer e : network.getEdges())
				unionNetwork.addEdge(idx++, network.getEndpoints(e), EdgeType.DIRECTED);
		}
		return unionNetwork;
	}

	protected final Map<Long, T> actionsCache =
			Collections.synchronizedMap(new LinkedHashMap<Long, T>() {
				private static final long serialVersionUID = 1L;
				private static final int cacheSize = 20000;

				@Override
				protected boolean removeEldestEntry(Map.Entry<Long, T> e) {
					return size() >= cacheSize;
				}
			});

	public static TemporalInfluenceModel<?> createModel(ModelVariant.ModelType type, Map<Long, TimeSeriesData> actions,
			double[] topicPriorParam) {
		if (type == ModelVariant.ModelType.SIMILARITY)
			return new SimilarityInfluenceModel(actions, topicPriorParam);
		if (type == ModelVariant.ModelType.CALIBRATED)
			return new CalibratedInfluenceModel(actions);
		return new GrangerInfluenceModel(type, actions, topicPriorParam.length);
	}

}
