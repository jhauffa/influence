package edu.tum.cs.influence.micro;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.atomic.AtomicLong;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;

import com.carrotsearch.hppc.DoubleArrayList;

import edu.tum.cs.math.MultipleHypothesisTest;
import edu.tum.cs.math.dist.BetaDistribution;
import edu.tum.cs.math.dist.ConstrainedDirichletDistribution;
import edu.tum.cs.math.dist.DirichletDistribution;
import edu.tum.cs.math.dist.DiscreteDistributionSampler;
import edu.tum.cs.math.dist.HistogramBin;
import edu.tum.cs.math.dist.HistogramImpl.QuantizedProbabilityHistogram;
import edu.tum.cs.util.ExperimentConfiguration;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

public class SimilarityInfluenceModel extends TemporalInfluenceModel<ChangeTimeSeriesData> {

	private static class ParameterHistogram extends QuantizedProbabilityHistogram {
		private final Map<Integer, Mean> binAvgStdDev = new HashMap<Integer, Mean>();

		public ParameterHistogram(double quantizationFactor) {
			super(quantizationFactor);
		}

		public int addValue(double value, double stdDev) {
			int binIdx = super.addValue(value);
			Mean avgStdDev = binAvgStdDev.get(binIdx);
			if (avgStdDev == null) {
				avgStdDev = new Mean();
				binAvgStdDev.put(binIdx, avgStdDev);
			}
			avgStdDev.increment(stdDev);
			return binIdx;
		}

		public String toString(String agg) {
			StringBuilder sb = new StringBuilder();
			sb.append("aggregation;L < x;x <= U;count;prob;stddev\n");
			for (HistogramBin<Double> bin : this) {
				sb.append(agg).append(';').append(bin.lowerLimit).append(';').append(bin.upperLimit).append(';');
				sb.append(bin.count).append(';').append(bin.p).append(';').append(binAvgStdDev.get(bin.index));
				sb.append('\n');
			}
			return sb.toString();
		}
	}

	private static class ParameterStatistics {
		private static final int numBins = 25;
		private static final int numFitTestSamples = 250;

		private final ParameterHistogram histogram = new ParameterHistogram(numBins);
		private final DoubleArrayList values = new DoubleArrayList();	// for Beta fitting & goodness-of-fit test
		private BetaDistribution dist;
		private double fitP;

		public ParameterStatistics() {
		}

		public ParameterStatistics(Map<Long, SummaryStatistics> agg) {
			for (SummaryStatistics userStats : agg.values())
				addValue(userStats.getMean(), userStats.getStandardDeviation());
		}

		public void addValue(double value, double aggStdDev) {
			histogram.addValue(value, aggStdDev);
			values.add(value);
		}

		public void addMagnitude(ChangeTimeSeriesData data) {
			for (int i = 0; i < data.getLength(); i++)
				addValue(data.getMagnitude(i), 0.0);
		}

		public BetaDistribution fitBetaDistribution() {
			if (dist == null) {
				double[] data = values.toArray();
				values.clear();
				dist = BetaDistribution.fit(data, true);
				fitP = dist.testGoodnessOfFit(data, numFitTestSamples);
			}
			return dist;
		}

		public void printHisto(PrintWriter w, String agg) {
			w.println(histogram.toString(agg));
		}

		public void printStats(PrintWriter w, String agg) {
			BetaDistribution dist = fitBetaDistribution();
			w.println("aggregation;N;alpha;beta;mean;variance;p");
			w.println(agg + ";" + histogram.getSampleSize() + ";" + dist.getAlpha() + ";" + dist.getBeta() + ";" +
					dist.getNumericalMean() + ";" + dist.getNumericalVariance() + ";" + fitP);
		}
	}

	private static class ModelStatistics {
		private final ParameterStatistics statsSimilarity = new ParameterStatistics();
		private final Map<Long, SummaryStatistics> aggSimilarityInfluencer = new HashMap<Long, SummaryStatistics>();
		private final Map<Long, SummaryStatistics> aggSimilarityInfluencee = new HashMap<Long, SummaryStatistics>();

		private static void addToAggregation(Map<Long, SummaryStatistics> agg, long id, double value) {
			SummaryStatistics stats = agg.get(id);
			if (stats == null) {
				stats = new SummaryStatistics();
				agg.put(id, stats);
			}
			stats.addValue(value);
		}

		public void fitBeta() {
			// fit Beta distribution to reclaim space
			statsSimilarity.fitBetaDistribution();
		}

		public void addSimilarity(double similarity, long fromId, long toId) {
			statsSimilarity.addValue(similarity, 0.0);
			addToAggregation(aggSimilarityInfluencer, fromId, similarity);
			addToAggregation(aggSimilarityInfluencee, toId, similarity);
		}

		public void writeStatistics(PrintWriter w, String id) {
			// aggregate similarity statistics by influencer/influencee
			ParameterStatistics statsInfluencer = new ParameterStatistics(aggSimilarityInfluencer);
			ParameterStatistics statsInfluencee = new ParameterStatistics(aggSimilarityInfluencee);

			// write histograms
			statsSimilarity.printHisto(w, id + "-all");
			statsInfluencer.printHisto(w, id + "-influencers");
			statsInfluencee.printHisto(w, id + "-influencees");

			// write beta distribution parameters
			statsSimilarity.printStats(w, id + "-all");
			statsInfluencer.printStats(w, id + "-influencers");
			statsInfluencee.printStats(w, id + "-influencees");
		}
	}

	private static enum TestType {
		PERMUTATION (0), DRIFT_DIRECTION (1), DRIFT_CHANGE (2), DRIFT_UNCONSTRAINED (3);

		public final int index;

		private TestType(int index) {
			this.index = index;
		}

		@Override
		public String toString() {
			switch (this) {
			case PERMUTATION:
				return "permutation";
			case DRIFT_DIRECTION:
				return "drift-direction";
			case DRIFT_CHANGE:
				return "drift-change";
			case DRIFT_UNCONSTRAINED:
				return "drift-unconstrained";
			default:
				return "unknown";
			}
		}
	}
	private static final TestType defaultDriftTestType = TestType.DRIFT_UNCONSTRAINED;

	private static class EdgeData {
		private final double similarity;
		private final Explanation explanation;
		private double[] p;
		public boolean insufficientData;

		public EdgeData(double similarity, Explanation explanation) {
			this.similarity = similarity;
			this.explanation = explanation;
			p = new double[TestType.values().length];
		}

		public double getPValue(TestType type) {
			return p[type.index];
		}

		public void setPValue(TestType type, double p) {
			this.p[type.index] = p;
		}
	}

	private static final Logger logger = Logger.getLogger(SimilarityInfluenceModel.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(SimilarityInfluenceModel.class);

	private final Map<Long, TimeSeriesData> actions;
	private final double[] topicPriorParam;
	private final ParameterStatistics statsMagnitudeLevel = new ParameterStatistics();
	private final ParameterStatistics statsMagnitudeEdge = new ParameterStatistics();
	private final Map<ModelVariant, ModelStatistics> variantStats = new HashMap<ModelVariant, ModelStatistics>();
	private final AtomicLong numPermutationTests = new AtomicLong(0);
	private final AtomicLong numDriftTests = new AtomicLong(0);
	private long numEdges, totalDriftTests;

	public SimilarityInfluenceModel(Map<Long, TimeSeriesData> actions, double[] topicPriorParam) {
		this.actions = actions;
		this.topicPriorParam = topicPriorParam;
	}

	private ChangeTimeSeriesData getTransformedData(Long id, TimeSeriesData data) {
		ChangeTimeSeriesData tfData = actionsCache.get(id);
		if (tfData == null) {
			tfData = ChangeTimeSeriesData.transformData(data, null);
			// We can get away without synchronization, because we visit all edges of the union real and fake graphs in
			// buildSimilarityGraphs, before any subtasks are spawned. After that, the cache is only read from.
			actionsCache.put(id, tfData);
		}
		return tfData;
	}

	private class PermutationTest
			extends MultipleHypothesisTest<WeightedEdge<EdgeData>, Queue<Integer>> {
		private static final int numMcSamples = 2000;
		private static final int minPeriodsReplSampling = 9;
		private static final int minActivePeriods = 6;	// minimum number for permutation test at alpha = 0.05
		private final UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S);
		private final DiscreteDistributionSampler samp = new DiscreteDistributionSampler();
		private final WeightedDirectedGraph<Long, EdgeData> graph;
		private final boolean useSimInfluencer;

		public PermutationTest(double alpha, WeightedDirectedGraph<Long, EdgeData> graph, boolean useSimInfluencer,
				PermutationTest prevInst) {
			super(alpha, numMcSamples, true);
			this.graph = graph;
			this.useSimInfluencer = useSimInfluencer;
		}

		@Override
		protected Map<WeightedEdge<EdgeData>, Queue<Integer>> preparePerHypothesisData(
				List<WeightedEdge<EdgeData>> edges, int maxIter) {
			Map<WeightedEdge<EdgeData>, Queue<Integer>> data = new HashMap<WeightedEdge<EdgeData>, Queue<Integer>>();
			for (WeightedEdge<EdgeData> edge : edges) {
				Pair<Long> v = graph.getEndpoints(edge);
				TimeSeriesData srcData = actions.get(v.getFirst());
				if (useSimInfluencer)
					srcData = getTransformedData(v.getFirst(), srcData);
				int srcLen = srcData.getLength() - srcData.getNumMissing();
				if (srcLen < minPeriodsReplSampling) {
					Queue<Integer> permIdx = new LinkedList<Integer>();
					// assume that minPresentPeriods is chosen so that srcLen! <= Integer.MAX_VALUE
					int numPerm = (int) CombinatoricsUtils.factorial(srcLen);
					int maxSamples = Math.min(maxIter, numPerm);
					for (int idx : samp.sample1DUniformWithoutReplacement(numPerm, maxSamples))
						permIdx.add(idx);
					data.put(edge, permIdx);
				}
			}
			return data;
		}

		@Override
		protected int generateSamples(int n, WeightedEdge<EdgeData> edge, Queue<Integer>[] unused,
				Queue<Integer> permIdx) {
			Pair<Long> v = graph.getEndpoints(edge);
			TimeSeriesData srcData = actions.get(v.getFirst());
			if (useSimInfluencer)
				srcData = getTransformedData(v.getFirst(), srcData);
			ChangeTimeSeriesData dstData = getTransformedData(v.getSecond(), actions.get(v.getSecond()));

			int numExceedance = 0;
			int[] srcIndices = srcData.getIndices();
			TimeSeriesData permSrcData = new TimeSeriesData(srcData);
			for (int i = 0; i < n; i++) {
				if (permIdx != null) {
					// cannot permute in place, as the permutations reference the original order of elements
					Integer idx = permIdx.poll();
					if (idx == null)	// exhausted all permutations
						return -1;
					srcData.permute(srcIndices, idx, permSrcData);
				} else
					permSrcData.permute(srcIndices, rand);

				double similarity = dstData.computeSimilarity(permSrcData, minActivePeriods, null);
				if (Double.isNaN(similarity))
					return -1;
				if (similarity > edge.weight.similarity)
					numExceedance++;
			}
			return numExceedance;
		}
	}

	private class TopicDriftTest extends MultipleHypothesisTest<WeightedEdge<EdgeData>, ChangeTimeSeriesData> {
		private static final int numMcSamples = 2000;
		private final TestType samplingType;
		private final BetaDistribution distMagnitude;
		private final WeightedDirectedGraph<Long, EdgeData> graph;
		private final boolean useSimInfluencer;
		private final List<ChangeTimeSeriesData> cachedSyntheticActions;

		protected TopicDriftTest(TestType samplingType, double alpha, BetaDistribution distMagnitude,
				WeightedDirectedGraph<Long, EdgeData> graph, boolean useSimInfluencer, TopicDriftTest prevInst) {
			super(alpha, numMcSamples, true);
			this.samplingType = samplingType;
			this.distMagnitude = distMagnitude;
			this.graph = graph;
			this.useSimInfluencer = useSimInfluencer;

			if (prevInst != null)
				cachedSyntheticActions = prevInst.cachedSyntheticActions;
			else
				cachedSyntheticActions = new ArrayList<ChangeTimeSeriesData>();
		}

		@Override
		protected ChangeTimeSeriesData[] preparePerIterationData(int n, int baseIter,
				List<WeightedEdge<EdgeData>> edges) {
			ChangeTimeSeriesData[] syntheticActions = new ChangeTimeSeriesData[n];
			Long v2 = graph.getEndpoints(edges.get(0)).getSecond();
			TimeSeriesData dstData = actions.get(v2);
			ChangeTimeSeriesData dstDataTf = null;
			if (samplingType == TestType.DRIFT_DIRECTION)
				dstDataTf = getTransformedData(v2, dstData);

			for (int i = 0; i < n; i++) {
				if (cachedSyntheticActions.size() > (baseIter + i)) {
					syntheticActions[i] = cachedSyntheticActions.get(baseIter + i);
				} else {
					// generate synthetic "present" topic distributions for influencee (same for all edges!)
					TimeSeriesData syntheticDstData = new TimeSeriesData(dstData.getLength());
					syntheticDstData.addValue(null);
					for (int j = 0; j < dstData.getLength() - 1; j++) {
						double[] syntheticDstT2 = null;
						if (!dstData.isMissing(j + 1)) {	// do not generate synthetic data if real data is missing
							if ((samplingType == TestType.DRIFT_UNCONSTRAINED) || dstData.isMissing(j)) {
								DirichletDistribution dist = new DirichletDistribution(topicPriorParam);
								syntheticDstT2 = dist.sample();
							} else {
								double[] dstT1 = dstData.getValue(j);
								double magnitude;
								if (samplingType == TestType.DRIFT_CHANGE)
									magnitude = distMagnitude.sample();
								else	// TestType.DRIFT_DIRECTION
									magnitude = dstDataTf.getMagnitude(j + 1);

								double minComp = 1.0;
								for (double v : dstT1) {
									if (v < minComp)
										minComp = v;
								}
								magnitude *= 1.0 - minComp;

								ConstrainedDirichletDistribution dist =
										new ConstrainedDirichletDistribution(topicPriorParam, dstT1, 2.0 * magnitude);
								syntheticDstT2 = dist.sample();
							}
						}
						syntheticDstData.addValue(syntheticDstT2);	// setting the j+1'th value
					}
					syntheticActions[i] = ChangeTimeSeriesData.transformData(dstData, syntheticDstData);
					cachedSyntheticActions.add(syntheticActions[i]);
				}
			}
			return syntheticActions;
		}

		@Override
		protected int generateSamples(int n, WeightedEdge<EdgeData> edge, ChangeTimeSeriesData[] syntheticActions,
				ChangeTimeSeriesData unused) {
			Pair<Long> v = graph.getEndpoints(edge);
			TimeSeriesData srcData = actions.get(v.getFirst());
			if (useSimInfluencer)
				srcData = getTransformedData(v.getFirst(), srcData);

			int numExceedance = 0;
			for (int i = 0; i < n; i++) {
				double similarity = syntheticActions[i].computeSimilarity(srcData, 1, null);
				if (Double.isNaN(similarity))
					return -1;
				if (similarity > edge.weight.similarity)
					numExceedance++;
			}
			return numExceedance;
		}
	}

	private int buildSimilarityGraph(DirectedSparseGraph<Long, Integer> network,
			WeightedDirectedGraph<Long, EdgeData> simGraph, int baseIdx, boolean useSimInfluencer,
			ParameterStatistics statsMagnitude, boolean generateExplanations) {
		for (Integer idx : network.getEdges()) {
			Pair<Long> v = network.getEndpoints(idx);
			if (simGraph.findEdge(v.getFirst(), v.getSecond()) != null)
				continue;

			TimeSeriesData srcData = actions.get(v.getFirst());
			if (useSimInfluencer) {
				ChangeTimeSeriesData tfSrcData = getTransformedData(v.getFirst(), srcData);
				if (statsMagnitude != null)
					statsMagnitude.addMagnitude(tfSrcData);
				srcData = tfSrcData;
			}
			ChangeTimeSeriesData dstData = getTransformedData(v.getSecond(), actions.get(v.getSecond()));
			if (statsMagnitude != null)
				statsMagnitude.addMagnitude(dstData);

			Explanation explanation = null;
			if (generateExplanations)
				explanation = new Explanation();
			double similarity = dstData.computeSimilarity(srcData, 0, explanation);
			simGraph.addEdge(new WeightedEdge<EdgeData>(baseIdx++, new EdgeData(similarity, explanation)), v,
					EdgeType.DIRECTED);
		}
		return baseIdx;
	}

	private void runTests(final TestType type,
			final Map<ModelVariant.InfluenceType, WeightedDirectedGraph<Long, EdgeData>> unfilteredGraphs,
			final BetaDistribution distMagnitudeLevel, final BetaDistribution distMagnitudeEdge, ForkJoinPool pool,
			List<ForkJoinTask<?>> tasks, final AtomicLong numTestsCompleted) {
		// executing hypothesis tests for all influence types in a single task to be able to share data (in this case:
		// synthetic observations for topic drift test)
		Set<Long> vertices = new HashSet<Long>();
		for (WeightedDirectedGraph<Long, EdgeData> graph : unfilteredGraphs.values())
			vertices.addAll(graph.getVertices());
		for (final Long v2 : vertices) {
			tasks.add(pool.submit(new Runnable() {
				@Override
				public void run() {
					try {
						MultipleHypothesisTest<WeightedEdge<EdgeData>, ?> test = null;
						for (Map.Entry<ModelVariant.InfluenceType, WeightedDirectedGraph<Long, EdgeData>> e :
							unfilteredGraphs.entrySet()) {
							List<WeightedEdge<EdgeData>> curEdges = new ArrayList<WeightedEdge<EdgeData>>(
									e.getValue().getInEdges(v2));
							if (curEdges.isEmpty())
								continue;

							boolean useSimInfluencer = false;
							BetaDistribution distMagnitude = distMagnitudeLevel;
							if (e.getKey() == ModelVariant.InfluenceType.EDGE) {
								useSimInfluencer = true;
								distMagnitude = distMagnitudeEdge;
							}
							if (type == TestType.PERMUTATION) {
								test = new PermutationTest(baseSignificanceLevel, e.getValue(), useSimInfluencer,
										(PermutationTest) test);
							} else {
								test = new TopicDriftTest(type, baseSignificanceLevel, distMagnitude, e.getValue(),
										useSimInfluencer, (TopicDriftTest) test);
							}

							boolean[] insufficientData = new boolean[curEdges.size()];
							double[] pValues = test.performTests(curEdges, numTestsCompleted, insufficientData);
							for (int i = 0; i < pValues.length; i++) {
								EdgeData weight = curEdges.get(i).weight;
								weight.setPValue(type, pValues[i]);
								weight.insufficientData |= insufficientData[i];
							}
						}
					} catch (Throwable t) {
						logger.log(Level.SEVERE, "error in test task", t);
						throw t;
					}
				}
			}));
		}
	}

	private static WeightedDirectedGraph<Long, Double> filterInfluenceGraph(DirectedSparseGraph<Long, Integer> network,
			WeightedDirectedGraph<Long, EdgeData> resultGraph, WeightedDirectedGraph<Long, Explanation> expGraph,
			TestType type, boolean correctFdr) {
		WeightedDirectedGraph<Long, Double> graph = new WeightedDirectedGraph<Long, Double>();
		for (Long v2 : network.getVertices()) {
			List<WeightedEdge<EdgeData>> edges = null;
			if (correctFdr)
				edges = new ArrayList<WeightedEdge<EdgeData>>(network.inDegree(v2));

			for (Integer e : network.getInEdges(v2)) {
				Long v1 = network.getOpposite(v2, e);
				WeightedEdge<EdgeData> we = resultGraph.findEdge(v1, v2);
				if (we != null) {
					if (correctFdr) {
						edges.add(we);
					} else {
						if (we.weight.getPValue(type) < baseSignificanceLevel) {
							Pair<Long> v = new Pair<Long>(v1, v2);
							graph.addEdge(new WeightedEdge<Double>(we.id, we.weight.similarity), v, EdgeType.DIRECTED);
							if (expGraph != null) {
								expGraph.addEdge(new WeightedEdge<Explanation>(we.id,  we.weight.explanation), v,
										EdgeType.DIRECTED);
							}
						}
					}
				}
			}

			if (correctFdr) {
				double[] p = new double[edges.size()];
				for (int i = 0; i < p.length; i++)
					p[i] = edges.get(i).weight.getPValue(type);
				for (Integer idx : MultipleHypothesisTest.testSignificanceBH(p, baseSignificanceLevel)) {
					WeightedEdge<EdgeData> we = edges.get(idx);
					Long v1 = resultGraph.getOpposite(v2, we);
					Pair<Long> v = new Pair<Long>(v1, v2);
					graph.addEdge(new WeightedEdge<Double>(we.id, we.weight.similarity), v, EdgeType.DIRECTED);
					if (expGraph != null) {
						expGraph.addEdge(new WeightedEdge<Explanation>(we.id, we.weight.explanation), v,
								EdgeType.DIRECTED);
					}
				}
			}
		}
		return graph;
	}

	private static ModelStatistics filterInfluenceStats(DirectedSparseGraph<Long, Integer> network,
			WeightedDirectedGraph<Long, EdgeData> resultGraph) {
		ModelStatistics stats = new ModelStatistics();
		for (Integer e : network.getEdges()) {
			Pair<Long> v = network.getEndpoints(e);
			WeightedEdge<EdgeData> we = resultGraph.findEdge(v.getFirst(), v.getSecond());
			if (we != null)
				stats.addSimilarity(we.weight.similarity, v.getFirst(), v.getSecond());
		}
		stats.fitBeta();
		return stats;
	}

	private static WeightedDirectedGraph<Long, Double> weightedGraphIntersection(
			WeightedDirectedGraph<Long, Double> graph1, WeightedDirectedGraph<Long, Double> graph2) {
		WeightedDirectedGraph<Long, Double> intersection = new WeightedDirectedGraph<Long, Double>();
		for (WeightedEdge<Double> edge : graph1.getEdges()) {
			if (graph2.containsEdge(edge))	// assumes consistent edge indices
				intersection.addEdge(edge, graph1.getEndpoints(edge), EdgeType.DIRECTED);
		}
		return intersection;
	}

	private static void buildVariantInfluenceGraphs(Map<ModelVariant, WeightedDirectedGraph<Long, Double>> graphs,
			Map<ModelVariant, WeightedDirectedGraph<Long, Explanation>> explanations, ModelVariant variant,
			DirectedSparseGraph<Long, Integer> network, EnumSet<TestType> driftTestTypes,
			WeightedDirectedGraph<Long, EdgeData> resultGraph) {
		WeightedDirectedGraph<Long, Explanation> expGraph = null;
		if (explanations != null) {
			expGraph = new WeightedDirectedGraph<Long, Explanation>();
			explanations.put(variant, expGraph);
		}

		WeightedDirectedGraph<Long, Double> permSubgraph = filterInfluenceGraph(network, resultGraph, null,
				TestType.PERMUTATION, false);
		graphs.put(variant.subModel(TestType.PERMUTATION.toString()), permSubgraph);
		for (TestType type : driftTestTypes) {
			WeightedDirectedGraph<Long, Double> driftSubgraph = filterInfluenceGraph(network, resultGraph, null,
							type, false);
			graphs.put(variant.subModel(type.toString()), driftSubgraph);
			graphs.put(variant.subModel(type + "-uncorrected"), weightedGraphIntersection(permSubgraph, driftSubgraph));
		}

		permSubgraph = filterInfluenceGraph(network, resultGraph, expGraph, TestType.PERMUTATION, true);
		for (TestType type : driftTestTypes) {
			WeightedDirectedGraph<Long, Double> driftSubgraph = filterInfluenceGraph(network, resultGraph,
						expGraph, type, true);
			ModelVariant correctedVariant = variant;
			if (type != defaultDriftTestType)
				correctedVariant = variant.subModel(type + "-corrected");
			graphs.put(correctedVariant, weightedGraphIntersection(permSubgraph, driftSubgraph));
		}
	}

	private void buildSimilarityGraphs(
			Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> networks,
			Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> fakeNetworks,
			Map<ModelVariant.InfluenceType, WeightedDirectedGraph<Long, EdgeData>> unfilteredGraphs,
			boolean generateExplanations) {
		DirectedSparseGraph<Long, Integer> unionNetwork = buildUnionNetwork(networks.values());
		DirectedSparseGraph<Long, Integer> unionFakeNetwork = buildUnionNetwork(fakeNetworks.values());
		numEdges = 0;
		for (ModelVariant.InfluenceType influence : ModelVariant.InfluenceType.values()) {
			boolean useSimInfluencer = false;
			ParameterStatistics statsMagnitude = statsMagnitudeLevel;
			if (influence == ModelVariant.InfluenceType.EDGE) {
				useSimInfluencer = true;
				statsMagnitude = statsMagnitudeEdge;
			}

			WeightedDirectedGraph<Long, EdgeData> simGraph = new WeightedDirectedGraph<Long, EdgeData>();
			int edgeIdx = 0;
			edgeIdx = buildSimilarityGraph(unionNetwork, simGraph, edgeIdx, useSimInfluencer, statsMagnitude,
					generateExplanations);
			edgeIdx = buildSimilarityGraph(unionFakeNetwork, simGraph, edgeIdx, useSimInfluencer, null, false);
			numEdges += edgeIdx;
			unfilteredGraphs.put(influence, simGraph);

			for (Map.Entry<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> e : networks.entrySet()) {
				ModelVariant variant = new ModelVariant(ModelVariant.ModelType.SIMILARITY, influence, e.getKey());
				variantStats.put(variant, filterInfluenceStats(e.getValue(), simGraph));
			}
		}
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
		// compute similarity for edges of union graph
		Map<ModelVariant.InfluenceType, WeightedDirectedGraph<Long, EdgeData>> unfilteredGraphs =
				new HashMap<ModelVariant.InfluenceType, WeightedDirectedGraph<Long, EdgeData>>();
		buildSimilarityGraphs(networks, fakeNetworks, unfilteredGraphs, (explanations != null));

		// fit Beta distribution to measured influence magnitudes for generation of synthetic data in topic drift test
		BetaDistribution distMagnitudeLevel = statsMagnitudeLevel.fitBetaDistribution();
		BetaDistribution distMagnitudeEdge = statsMagnitudeEdge.fitBetaDistribution();

		// run MC hypothesis tests
		List<ForkJoinTask<?>> tasks = new ArrayList<ForkJoinTask<?>>();
		ForkJoinPool pool = new ForkJoinPool(cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_THREADS,
				Runtime.getRuntime().availableProcessors()));
		runTests(TestType.PERMUTATION, unfilteredGraphs, null, null, pool, tasks, numPermutationTests);
		EnumSet<TestType> driftTestTypes;
		if (cfg.getLocalBooleanProperty("performAllDriftTests", false)) {
			driftTestTypes = EnumSet.allOf(TestType.class);
			driftTestTypes.remove(TestType.PERMUTATION);
		} else {
			driftTestTypes = EnumSet.noneOf(TestType.class);
			driftTestTypes.add(defaultDriftTestType);
		}
		totalDriftTests = numEdges * driftTestTypes.size();
		for (TestType type : driftTestTypes)
			runTests(type, unfilteredGraphs, distMagnitudeLevel, distMagnitudeEdge, pool, tasks, numDriftTests);

		// Wait for tests to complete. Cannot use a CompletionService here, because we need all subtasks for a
		// particular graph to be complete before further processing is possible.
		for (ForkJoinTask<?> task : tasks)
			task.join();
		pool.shutdown();

		// generate filtered graphs according to corrected/uncorrected p-values
		Map<ModelVariant, WeightedDirectedGraph<Long, Double>> graphs =
				new HashMap<ModelVariant, WeightedDirectedGraph<Long, Double>>();
		for (Map.Entry<ModelVariant.InfluenceType, WeightedDirectedGraph<Long, EdgeData>> e :
				unfilteredGraphs.entrySet()) {
			for (Map.Entry<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> ee : networks.entrySet()) {
				ModelVariant variant = new ModelVariant(ModelVariant.ModelType.SIMILARITY, e.getKey(), ee.getKey());
				buildVariantInfluenceGraphs(graphs, explanations, variant, ee.getValue(), driftTestTypes, e.getValue());
				numEarlyRejections.put(variant, countEarlyRejections(ee.getValue(), e.getValue()));
			}
			for (Map.Entry<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> ee : fakeNetworks.entrySet()) {
				ModelVariant variant = new ModelVariant(ModelVariant.ModelType.SIMILARITY, e.getKey(), ee.getKey());
				buildVariantInfluenceGraphs(graphs, explanations, variant.fake(), ee.getValue(), driftTestTypes,
						e.getValue());
			}
		}
		return graphs;
	}

	@Override
	public void writeStatistics(PrintWriter w) {
		// write histogram and beta distribution parameters of similarity
		for (Map.Entry<ModelVariant, ModelStatistics> e : variantStats.entrySet())
			e.getValue().writeStatistics(w, e.getKey().getFilePrefix());

		// write histogram and beta distribution parameters of magnitude
		statsMagnitudeLevel.printHisto(w, "magnitude-level");
		statsMagnitudeEdge.printHisto(w, "magnitude-edge");
		statsMagnitudeLevel.printStats(w, "magnitude-level");
		statsMagnitudeEdge.printStats(w, "magnitude-edge");
	}

	@Override
	public Map<String, long[]> reportProgress() {
		Map<String, long[]> progress = new HashMap<String, long[]>();
		progress.put("permutation tests done", new long[] { numPermutationTests.get(), numEdges });
		progress.put("drift tests done", new long[] { numDriftTests.get(), totalDriftTests });
		return progress;
	}

}
