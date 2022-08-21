package edu.tum.cs.influence.micro;

import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import edu.tum.cs.graph.GraphCreator;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

public class CalibratedInfluenceModel extends TemporalInfluenceModel<ChangeTimeSeriesData> {

	private static class EdgeData {
		public final double similarity;
		public final boolean insufficientData;

		public EdgeData(double similarity) {
			if (Double.isNaN(similarity)) {
				this.similarity = 0.0;
				this.insufficientData = true;
			} else {
				this.similarity = similarity;
				this.insufficientData = false;
			}
		}
	}

	private static final int numSd = 3;

	private final Map<Long, TimeSeriesData> actions;
	private final Map<ModelVariant, SummaryStatistics> variantStats = new HashMap<ModelVariant, SummaryStatistics>();

	public CalibratedInfluenceModel(Map<Long, TimeSeriesData> actions) {
		this.actions = actions;
	}

	private ChangeTimeSeriesData getTransformedData(Long id, TimeSeriesData data) {
		ChangeTimeSeriesData tfData = actionsCache.get(id);
		if (tfData == null) {
			tfData = ChangeTimeSeriesData.transformData(data, null);
			actionsCache.put(id, tfData);
		}
		return tfData;
	}

	private int buildSimilarityGraph(DirectedSparseGraph<Long, Integer> network,
			WeightedDirectedGraph<Long, EdgeData> simGraph, int baseIdx, boolean useSimInfluencer) {
		for (Integer idx : network.getEdges()) {
			Pair<Long> v = network.getEndpoints(idx);
			if (simGraph.findEdge(v.getFirst(), v.getSecond()) != null)
				continue;

			TimeSeriesData srcData = actions.get(v.getFirst());
			if (useSimInfluencer) {
				ChangeTimeSeriesData tfSrcData = getTransformedData(v.getFirst(), srcData);
				srcData = tfSrcData;
			}
			ChangeTimeSeriesData dstData = getTransformedData(v.getSecond(), actions.get(v.getSecond()));
			double similarity = dstData.computeSimilarity(srcData, 1, null);
			simGraph.addEdge(new WeightedEdge<EdgeData>(baseIdx++, new EdgeData(similarity)), v, EdgeType.DIRECTED);
		}
		return baseIdx;
	}

	private double computeThreshold(WeightedDirectedGraph<Long, EdgeData> simGraph,
			DirectedSparseGraph<Long, Integer> network, ModelVariant variant) {
		SummaryStatistics stats = new SummaryStatistics();
		variantStats.put(variant, stats);
		for (Integer idx : network.getEdges()) {
			Pair<Long> v = network.getEndpoints(idx);
			WeightedEdge<EdgeData> e = simGraph.findEdge(v.getFirst(), v.getSecond());
			stats.addValue(e.weight.similarity);
		}
		return stats.getMean() + (numSd * stats.getStandardDeviation());
	}

	private static WeightedDirectedGraph<Long, Double> filterSimilarityGraph(
			WeightedDirectedGraph<Long, EdgeData> simGraph, DirectedSparseGraph<Long, Integer> network,
			double threshold) {
		WeightedDirectedGraph<Long, Double> influenceGraph = new WeightedDirectedGraph<Long, Double>();
		for (Integer idx : network.getEdges()) {
			Pair<Long> v = network.getEndpoints(idx);
			WeightedEdge<EdgeData> e = simGraph.findEdge(v.getFirst(), v.getSecond());
			if (e.weight.similarity >= threshold)
				influenceGraph.addEdge(new WeightedEdge<Double>(e.id, e.weight.similarity), v, EdgeType.DIRECTED);
		}
		return influenceGraph;
	}

	private static int countEarlyRejections(DirectedSparseGraph<Long, Integer> network,
			WeightedDirectedGraph<Long, EdgeData> resultGraph) {
		int n = 0;
		for (Integer e : network.getEdges()) {
			Pair<Long> v = network.getEndpoints(e);
			WeightedEdge<EdgeData> we = resultGraph.findEdge(v.getFirst(), v.getSecond());
			if ((we != null) && we.weight.insufficientData)
				n++;
		}
		return n;
	}

	@Override
	public Map<ModelVariant, WeightedDirectedGraph<Long, Double>> buildInfluenceGraphs(
			Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> networks,
			Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> fakeNetworks,
			Map<ModelVariant, WeightedDirectedGraph<Long, Explanation>> explanations) {
		Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> calNetworks =
				new HashMap<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>>();
		for (Map.Entry<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> e : networks.entrySet()) {
			int numCalEdges = Math.max(e.getValue().getEdgeCount(), 5);
			calNetworks.put(e.getKey(), GraphCreator.sampleComplementGraph(e.getValue(), numCalEdges));
		}

		Map<ModelVariant, WeightedDirectedGraph<Long, Double>> influenceGraphs =
				new HashMap<ModelVariant, WeightedDirectedGraph<Long, Double>>();
		for (ModelVariant.InfluenceType influence : ModelVariant.InfluenceType.values()) {
			boolean useSimInfluencer = (influence == ModelVariant.InfluenceType.EDGE);

			// generate union graph of all real and fake networks weighted with similarity
			WeightedDirectedGraph<Long, EdgeData> similarityGraph = new WeightedDirectedGraph<Long, EdgeData>();
			int edgeIdx = 0;
			for (DirectedSparseGraph<Long, Integer> g : networks.values())
				edgeIdx = buildSimilarityGraph(g, similarityGraph, edgeIdx, useSimInfluencer);
			for (DirectedSparseGraph<Long, Integer> g : calNetworks.values())
				edgeIdx = buildSimilarityGraph(g, similarityGraph, edgeIdx, useSimInfluencer);
			for (DirectedSparseGraph<Long, Integer> g : fakeNetworks.values())
				edgeIdx = buildSimilarityGraph(g, similarityGraph, edgeIdx, useSimInfluencer);

			// compute threshold and apply it to the similarity-weighted networks
			for (ModelVariant.NetworkType network : networks.keySet()) {
				ModelVariant variant = new ModelVariant(ModelVariant.ModelType.CALIBRATED, influence, network);
				double threshold = computeThreshold(similarityGraph, calNetworks.get(network), variant);
				WeightedDirectedGraph<Long, Double> influenceGraph = filterSimilarityGraph(similarityGraph,
						networks.get(network), threshold);
				influenceGraphs.put(variant, influenceGraph);
				numEarlyRejections.put(variant, countEarlyRejections(networks.get(network), similarityGraph));
				influenceGraph = filterSimilarityGraph(similarityGraph, fakeNetworks.get(network), threshold);
				influenceGraphs.put(variant.fake(), influenceGraph);
			}
		}
		return influenceGraphs;
	}

	@Override
	public void writeStatistics(PrintWriter w) {
		w.println(ModelVariant.getCsvHeader() + ";mean;sd");
		for (Map.Entry<ModelVariant, SummaryStatistics> e : variantStats.entrySet())
			w.println(e.getKey().toString() + ";" + e.getValue().getMean() + ";" + e.getValue().getStandardDeviation());
	}

	@Override
	public Map<String, long[]> reportProgress() {
		return Collections.emptyMap();
	}

}
