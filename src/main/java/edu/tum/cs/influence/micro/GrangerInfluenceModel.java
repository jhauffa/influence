package edu.tum.cs.influence.micro;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.atomic.AtomicLong;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.stat.descriptive.SynchronizedSummaryStatistics;

import com.google.common.collect.Iterables;

import edu.tum.cs.math.dist.DiscreteDistribution;
import edu.tum.cs.util.ExperimentConfiguration;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

public class GrangerInfluenceModel extends TemporalInfluenceModel<TimeSeriesData> {

	private static final Logger logger = Logger.getLogger(GrangerInfluenceModel.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(GrangerInfluenceModel.class);

	private final ModelVariant.ModelType type;
	private final Map<Long, TimeSeriesData> actions;
	private final boolean transform;
	private final boolean useSimInfluencee;
	private final int numDim;
	private final int absMinObservations;
	private final double relMinObservations;

	private final AtomicLong numTests = new AtomicLong(0);
	private final AtomicLong numSkippedTests = new AtomicLong(0);
	private final AtomicLong numTestsRowsEliminated = new AtomicLong(0);
	private final AtomicLong numTestErrors = new AtomicLong(0);
	private long totalTests;
	private final AtomicLong numZeroRows = new AtomicLong(0);
	private final AtomicLong totalRows = new AtomicLong(0);
	private final AtomicLong numZeroesReplaced = new AtomicLong(0);
	private final AtomicLong totalValues = new AtomicLong(0);
	private final SynchronizedSummaryStatistics survivingTimeSteps = new SynchronizedSummaryStatistics();

	public GrangerInfluenceModel(ModelVariant.ModelType type, Map<Long, TimeSeriesData> actions, int numDim) {
		this.type = type;
		this.actions = actions;
		this.numDim = numDim;
		absMinObservations = cfg.getLocalIntProperty("absMinObservations");
		relMinObservations = cfg.getLocalDoubleProperty("relMinObservations");
		transform = (type != ModelVariant.ModelType.GRANGER);
		useSimInfluencee = (type == ModelVariant.ModelType.GRANGER_SIM_TF);
	}

	private TimeSeriesData getTransformedTimeSeries(Long id, boolean similarity) {
		if (!similarity)
			return actions.get(id);

		TimeSeriesData data = actionsCache.get(id);
		if (data == null) {
			// Transform time series by computing the additive change between two successive time periods.
			TimeSeriesData origData = actions.get(id);
			data = new TimeSeriesData(origData.getLength());
			data.addValue(null);
			for (int i = 0; i < origData.getLength() - 1; i++) {
				if (!origData.isMissing(i) && !origData.isMissing(i + 1)) {
					ChangeMeasurement change = ChangeMeasurement.computeChange(origData.getValue(i),
							origData.getValue(i + 1));
					data.addValue(change.direction);
				} else
					data.addValue(null);
			}
			actionsCache.put(id, data);
		}

		return data;
	}

	private static boolean[] findLongestSubseq(TimeSeriesData src, TimeSeriesData dst) {
		boolean[] mask = new boolean[dst.getLength()];
		int bestLen = 0, bestStart = 0, bestEnd = 0;
		int curLen = 0, curStart = 0;
		for (int i = 0; i < src.getLength(); i++) {
			if (!src.isMissing(i) && !dst.isMissing(i)) {
				curLen++;
			} else {
				int curEnd = i;
				if ((curLen > 0) && !dst.isMissing(i)) {
					curLen++;
					curEnd++;
				}
				if (curLen > bestLen) {
					bestLen = curLen;
					bestStart = curStart;
					bestEnd = curEnd;
				}
				curStart = i + 1;
				curLen = 0;
			}
		}
		if (curLen > bestLen) {
			bestStart = curStart;
			bestEnd = src.getLength();
		}

		for (int i = 0; i < bestStart; i++)
			mask[i] = true;
		for (int i = bestEnd; i < src.getLength(); i++)
			mask[i] = true;
		return mask;
	}

	private class ProcessVertex extends RecursiveAction {
		private static final long serialVersionUID = 1L;

		private final DirectedSparseGraph<Long, Integer> unionNetwork;
		private final Long v2;
		private final ModelVariant.InfluenceType influence;
		private final Map<ModelVariant, WeightedDirectedGraph<Long, Double>> graphs;
		private final Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> networks;
		private final Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> fakeNetworks;
		private final Map<ModelVariant, WeightedDirectedGraph<Long, Explanation>> explanations;

		public ProcessVertex(DirectedSparseGraph<Long, Integer> unionNetwork, Long v2,
				ModelVariant.InfluenceType influence, Map<ModelVariant, WeightedDirectedGraph<Long, Double>> graphs,
				Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> networks,
				Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> fakeNetworks,
				Map<ModelVariant, WeightedDirectedGraph<Long, Explanation>> explanations) {
			this.unionNetwork = unionNetwork;
			this.v2 = v2;
			this.influence = influence;
			this.graphs = graphs;
			this.networks = networks;
			this.fakeNetworks = fakeNetworks;
			this.explanations = explanations;
		}

		@Override
		protected void compute() {
			boolean useSimInfluencer = (influence == ModelVariant.InfluenceType.EDGE);
			Map<Integer, GrangerCausality> results = new HashMap<Integer, GrangerCausality>();
			for (Integer e : unionNetwork.getInEdges(v2)) {
				Long v1 = unionNetwork.getOpposite(v2, e);
				TimeSeriesData srcData = getTransformedTimeSeries(v1, useSimInfluencer);
				TimeSeriesData dstData = getTransformedTimeSeries(v2, useSimInfluencee);

				// eliminate missing values if possible
				boolean[] mask = findLongestSubseq(srcData, dstData);
				int numRemainingSteps = 0;
				for (boolean isEliminated : mask) {
					if (!isEliminated)
						numRemainingSteps++;
				}
				TimeSeriesData tfSrcData = new TimeSeriesData(numRemainingSteps);
				TimeSeriesData tfDstData = new TimeSeriesData(numRemainingSteps);
				for (int i = 0; i < mask.length; i++) {
					if (!mask[i]) {
						tfSrcData.addValue(srcData.getValue(i));
						tfDstData.addValue(dstData.getValue(i));
					}
				}
				// final value of influencer is irrelevant and may be missing -> replace
				if ((numRemainingSteps > 1) && tfSrcData.isMissing(numRemainingSteps - 1))
					tfSrcData.setValue(numRemainingSteps - 1, tfSrcData.getValue(numRemainingSteps - 2));

				int numSteps = srcData.getLength();
				if (useSimInfluencer || useSimInfluencee)
					numSteps--;
				double relRemainingSteps = (double) numRemainingSteps / numSteps;
				survivingTimeSteps.addValue(relRemainingSteps);

				GrangerCausality result = new GrangerCausality();
				if ((numRemainingSteps > absMinObservations) && (relRemainingSteps > relMinObservations)) {
					try {
						result = GrangerCausality.compute(tfSrcData.toRowMajorArray(), tfDstData.toRowMajorArray(),
								baseSignificanceLevel, numDim, transform, (explanations != null));
					} catch (RuntimeException ex) {
						numTestErrors.incrementAndGet();
						logger.log(Level.WARNING, "GC hypothesis test (" + type + ") failed for edge " +
								v1 + " -> " + v2 + ": " + ex.getMessage());
					}
				}
				results.put(e, result);

				if (result.isEmpty())
					numSkippedTests.incrementAndGet();
				numTests.incrementAndGet();

				GrangerCausality.PreprocessingStatistics stats = result.getPreprocessingStatistics();
				if (stats.allRowsEliminated)
					numTestsRowsEliminated.incrementAndGet();
				numZeroRows.addAndGet(stats.numZeroRows);
				totalRows.addAndGet(stats.totalRows);
				numZeroesReplaced.addAndGet(stats.numZeroesReplaced);
				totalValues.addAndGet(stats.totalValues);
			}

			synchronized (graphs) {
				for (Map.Entry<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> ee : networks.entrySet())
					filterInfluenceGraph(results, ee.getValue(), new ModelVariant(type, influence, ee.getKey()));
				for (Map.Entry<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> ee :
						fakeNetworks.entrySet())
					filterInfluenceGraph(results, ee.getValue(), new ModelVariant(type, influence, ee.getKey()).fake());
			}
		}

		private WeightedDirectedGraph<Long, Double> getGraph(ModelVariant variant) {
			WeightedDirectedGraph<Long, Double> graph = graphs.get(variant);
			if (graph == null) {
				graph = new WeightedDirectedGraph<Long, Double>();
				graphs.put(variant, graph);
			}
			return graph;
		}

		private WeightedDirectedGraph<Long, Explanation> getExplanationGraph(ModelVariant variant) {
			if (explanations == null)
				return null;
			WeightedDirectedGraph<Long, Explanation> graph = explanations.get(variant);
			if (graph == null) {
				graph = new WeightedDirectedGraph<Long, Explanation>();
				explanations.put(variant, graph);
			}
			return graph;
		}

		private Explanation buildExplanation(GrangerCausality gc) {
			Explanation exp = new Explanation();
			exp.timePeriodIdx = 0;
			exp.topicsInfluencer = gc.getSrcAbsCoeff().clone();
			DiscreteDistribution.normalize(exp.topicsInfluencer);
			exp.topicsInfluencee = gc.getDstAbsCoeff().clone();
			DiscreteDistribution.normalize(exp.topicsInfluencee);
			exp.localMagnitude = 0.0;
			exp.globalMagnitude = gc.getMagnitude();
			return exp;
		}

		private void filterInfluenceGraph(Map<Integer, GrangerCausality> results,
				DirectedSparseGraph<Long, Integer> network, ModelVariant variant) {
			WeightedDirectedGraph<Long, Double> graph = getGraph(variant);
			WeightedDirectedGraph<Long, Explanation> expGraph = getExplanationGraph(variant);
			WeightedDirectedGraph<Long, Double> uncorrGraph = getGraph(variant.subModel("uncorrected"));
			Map<Long, GrangerCausality> inEdges = new HashMap<Long, GrangerCausality>();
			for (Map.Entry<Integer, GrangerCausality> e : results.entrySet()) {
				Pair<Long> v = unionNetwork.getEndpoints(e.getKey());
				if (network.findEdge(v.getFirst(), v.getSecond()) != null) {
					GrangerCausality result = e.getValue();
					if (result.isEmpty()) {
						Integer n = numEarlyRejections.get(variant);
						if (n == null)
							n = 0;
						numEarlyRejections.put(variant, n + 1);
					}

					inEdges.put(v.getFirst(), result);
					if (result.isSignificantJJ()) {
						uncorrGraph.addEdge(new WeightedEdge<Double>(e.getKey(), result.getMagnitude()), v,
								EdgeType.DIRECTED);
						if (result.isSignificantJJ(baseSignificanceLevel / network.inDegree(v.getSecond()))) {
							graph.addEdge(new WeightedEdge<Double>(e.getKey(), result.getMagnitude()), v,
									EdgeType.DIRECTED);
							if (expGraph != null) {
								expGraph.addEdge(new WeightedEdge<Explanation>(e.getKey(), buildExplanation(result)), v,
										EdgeType.DIRECTED);
							}
						}
					}
				}
			}

			WeightedDirectedGraph<Long, Double> nodewiseGraph = getGraph(variant.subModel("nodewise"));
			expGraph = getExplanationGraph(variant.subModel("nodewise"));
			GrangerCausality.testGroupSignificanceJJ(inEdges, baseSignificanceLevel);
			for (Map.Entry<Long, GrangerCausality> e : inEdges.entrySet()) {
				Pair<Long> v = new Pair<Long>(e.getKey(), v2);
				Integer edge = network.findEdge(v.getFirst(), v.getSecond());
				nodewiseGraph.addEdge(new WeightedEdge<Double>(edge, e.getValue().getMagnitude()), v,
						EdgeType.DIRECTED);
				if (expGraph != null) {
					expGraph.addEdge(new WeightedEdge<Explanation>(edge, buildExplanation(e.getValue())), v,
							EdgeType.DIRECTED);
				}
			}
		}
	}

	@Override
	public Map<ModelVariant, WeightedDirectedGraph<Long, Double>> buildInfluenceGraphs(
			Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> networks,
			Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> fakeNetworks,
			Map<ModelVariant, WeightedDirectedGraph<Long, Explanation>> explanations) {
		Map<ModelVariant, WeightedDirectedGraph<Long, Double>> graphs =
				new HashMap<ModelVariant, WeightedDirectedGraph<Long, Double>>();
		DirectedSparseGraph<Long, Integer> unionNetwork = buildUnionNetwork(Iterables.concat(networks.values(),
				fakeNetworks.values()));
		totalTests += unionNetwork.getEdgeCount() * (useSimInfluencee ? 2 : 1);

		ForkJoinPool pool = new ForkJoinPool(cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_THREADS,
				Runtime.getRuntime().availableProcessors()));
		List<ForkJoinTask<?>> tasks = new ArrayList<ForkJoinTask<?>>();
		for (ModelVariant.InfluenceType influence : ModelVariant.InfluenceType.values()) {
			if ((influence == ModelVariant.InfluenceType.EDGE) && !useSimInfluencee)
				continue;
			for (final Long v2 : unionNetwork.getVertices()) {
				tasks.add(pool.submit(new ProcessVertex(unionNetwork, v2, influence, graphs, networks, fakeNetworks,
						explanations)));
			}
		}
		for (ForkJoinTask<?> t : tasks)
			t.join();
		pool.shutdown();
		return graphs;
	}

	@Override
	public void writeStatistics(PrintWriter w) {
		w.println(";count;total;percent;");
		w.println("rows discarded;" + numZeroRows.get() + ";" + totalRows.get() + ";" +
				(((double) numZeroRows.get() / totalRows.get()) * 100.0) + ";");
		w.println("zeroes replaced;" + numZeroesReplaced.get() + ";" + totalValues.get() + ";" +
				(((double) numZeroesReplaced.get() / totalValues.get()) * 100.0) + ";");
		w.println("edges skipped;" + numSkippedTests.get() + ";" + totalTests + ";" +
				(((double) numSkippedTests.get() / totalTests) * 100.0) + ";");
		w.println("edges eliminated;" + numTestsRowsEliminated.get() + ";" + totalTests + ";" +
				(((double) numTestsRowsEliminated.get() / totalTests) * 100.0) + ";");
		w.println("errors;" + numTestErrors.get() + ";" + totalTests + ";" +
				(((double) numTestErrors.get() / totalTests) * 100.0) + ";");
		w.println(";;;;");
		w.println(";min;max;mean;sd");
		w.println("surviving time steps;" + survivingTimeSteps.getMin() + ";" + survivingTimeSteps.getMax() + ";" +
				survivingTimeSteps.getMean() + ";" + survivingTimeSteps.getStandardDeviation());
	}

	@Override
	public Map<String, long[]> reportProgress() {
		Map<String, long[]> progress = new HashMap<String, long[]>();
		progress.put(type.toString() + " tests done", new long[] { numTests.get(), totalTests });
		return progress;
	}

}
