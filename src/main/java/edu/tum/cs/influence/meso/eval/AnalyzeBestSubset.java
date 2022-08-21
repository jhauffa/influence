package edu.tum.cs.influence.meso.eval;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.StorelessUnivariateStatistic;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import edu.tum.cs.influence.meso.EntityType;
import edu.tum.cs.influence.meso.ExperimentSeries;
import edu.tum.cs.influence.meso.IndicatorFunction;
import edu.tum.cs.influence.meso.InfluenceExperimentsAnova;
import edu.tum.cs.influence.meso.ModelVariant;
import edu.tum.cs.influence.meso.ParameterFittingResult;
import edu.tum.cs.influence.meso.PredictionResult;
import edu.tum.cs.influence.meso.ResultFilter;
import edu.tum.cs.influence.meso.WeightFunction;
import edu.tum.cs.util.ExperimentConfiguration;

public class AnalyzeBestSubset {

	private static class SankeyAggregator implements ResultFilter<ParameterFittingResult> {
		private static final String[] coeffGroupsScope = { "personal", "relationship", "neighborhood", "medium" };
		private static final String[] coeffGroupsRole = { "inertia", "direct exposure", "indirect exposure" };

		private static final String[] coeffsNode = {
			"sent by ego (non-addr.)", "sent by ego (addr.)", "received by ego (non-addr.)", "received by ego (addr.)",
			"from neighbors to ego (addr.)", "from ego to neighbors (addr.)", "sent by neighbors (non-addr.)",
			"sent by neighbors (addr.)", "all (non-addr.)", "all (addr.)", "data-driven"
		};
		private static final int[] coeffScopeMapNode = { 0, 0, 0, 0, 2, 2, 2, 2, 3, 3, 3 };
		private static final int[] coeffRoleMapNode  = { 0, 0, 1, 1, 1, 0, 2, 2, 2, 2, 2 };
		private static final String[] coeffsEdge = {
			"from ego to alter (addr.)", "from alter to ego (addr.)", "sent by ego (non-addr.)", "sent by ego (addr.)",
			"received by ego (non-addr.)", "received by ego (addr.)", "sent by alter (non-addr.)",
			"sent by alter (addr.)", "from neighbors to ego (addr.)", "from ego to neighbors (addr.)",
			"sent by neighbors (non-addr.)", "sent by neighbors (addr.)", "all (non-addr.)", "all (addr.)",
			"data-driven"
		};
		private static final int[] coeffScopeMapEdge = { 1, 1, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3 };
		private static final int[] coeffRoleMapEdge  = { 0, 1, 0, 0, 1, 1, 2, 2, 1, 0, 2, 2, 2, 2, 2 };

		private final boolean nodes;
		private final StorelessUnivariateStatistic errorMean;
		private final StorelessUnivariateStatistic errorStdDev;
		private final int numCoeff;
		private final StorelessUnivariateStatistic[] coeffMean;
		private final StorelessUnivariateStatistic[] coeffVariance;
		private final StorelessUnivariateStatistic[] coeffStdDev;

		public SankeyAggregator(boolean nodes) {
			this.nodes = nodes;
			this.errorMean = new Mean();
			this.errorStdDev = new StandardDeviation();
			this.numCoeff = nodes ? coeffsNode.length : coeffsEdge.length;
			this.coeffMean = new StorelessUnivariateStatistic[numCoeff];
			this.coeffVariance = new StorelessUnivariateStatistic[numCoeff];
			this.coeffStdDev = new StorelessUnivariateStatistic[numCoeff];
			for (int i = 0; i < numCoeff; i++) {
				coeffMean[i] = new Mean();
				coeffVariance[i] = new Variance();
				coeffStdDev[i] = new StandardDeviation();
			}
		}

		@Override
		public boolean evaluate(Date timePeriodBase, ModelVariant variant, ParameterFittingResult result) {
			if (variant.getType() != ModelVariant.ModelType.SCIM)
				return false;
			if ((nodes == (variant.getEntity() == EntityType.NODE)) && IndicatorWeightPair.inBestSubset(variant)) {
				errorMean.increment(result.getTrainErrorMean());
				errorStdDev.increment(result.getTrainErrorMean());

				int i;
				double coeffSum = 0.0;
				for (i = 0; i < numCoeff - 1; i++) {
					double c = result.getModel().getCoefficients()[i];
					coeffMean[i].increment(c);
					coeffVariance[i].increment(c);
					coeffStdDev[i].increment(c);
					coeffSum += c;
				}
				double cd = 1.0 - coeffSum;
				coeffMean[i].increment(cd);
				coeffVariance[i].increment(cd);
				coeffStdDev[i].increment(cd);
			}
			return false;
		}

		private void writeCsvForGroup(File baseName, String[] coeffNames, String[] groupNames, int[] coeffGroupMap)
				throws IOException {
			double[] groupMeanSums = new double[groupNames.length];
			// compute sum of variances as per
			// http://stats.stackexchange.com/questions/25848/how-to-sum-a-standard-deviation
			double[] groupVarSums = new double[groupNames.length];
			for (int i = 0; i < numCoeff; i++) {
				int groupIdx = coeffGroupMap[i];
				groupMeanSums[groupIdx] += coeffMean[i].getResult();
				groupVarSums[groupIdx] += coeffVariance[i].getResult();
			}
			double overallMeanSum = 0.0;
			double overallVarSum = 0.0;
			for (int i = 0; i < groupNames.length; i++) {
				overallMeanSum += groupMeanSums[i];
				overallVarSum += groupVarSums[i];
			}

			PrintWriter w = new PrintWriter(new FileWriter(new File(baseName.getAbsolutePath() + ".csv")));
			try {
				w.println("variable;group;mean;std.dev");
				w.println("MSE;;" + errorMean.getResult() + ";" + errorStdDev.getResult());
				for (int i = 0; i < numCoeff; i++) {
					int groupIdx = coeffGroupMap[i];
					w.println(coeffNames[i] + ";" + groupNames[groupIdx] + ";" + coeffMean[i].getResult() + ";" +
							coeffStdDev[i].getResult());
				}
				for (int i = 0; i < groupNames.length; i++) {
					if (groupMeanSums[i] != 0.0)
						w.println(groupNames[i] + ";all;" + groupMeanSums[i] + ";" + Math.sqrt(groupVarSums[i]));
				}
				w.println("overall;;" + overallMeanSum + ";" + Math.sqrt(overallVarSum));
			} finally {
				w.close();
			}

			w = new PrintWriter(new FileWriter(new File(baseName.getAbsolutePath() + ".json")));
			try {
				w.println("var treeData =");
				w.printf(Locale.US, "  { \"name\": \"\", \"mean\": %g, \"stddev\": %g, \"children\": [\n",
						overallMeanSum, Math.sqrt(overallVarSum));
				for (int i = 0; i < groupNames.length; i++) {
					if (groupMeanSums[i] == 0.0)
						continue;
					w.printf(Locale.US, "    { \"name\": \"%s\", \"mean\": %g, \"stddev\": %g, \"children\": [\n",
							groupNames[i], groupMeanSums[i], Math.sqrt(groupVarSums[i]));
					for (int j = 0; j < numCoeff; j++) {
						if (coeffGroupMap[j] != i)
							continue;
						w.printf(Locale.US, "      { \"name\": \"%s\", \"mean\": %g, \"stddev\": %g },\n",
								coeffNames[j], coeffMean[j].getResult(), coeffStdDev[j].getResult());
					}
					w.println("    ] },");
				}
				w.println("  ] };");
			} finally {
				w.close();
			}
		}

		@Override
		public void writeCsv(String path) throws IOException {
			if (nodes) {
				writeCsvForGroup(new File(path, "coeff-nodes-scope"), coeffsNode, coeffGroupsScope, coeffScopeMapNode);
				writeCsvForGroup(new File(path, "coeff-nodes-role"), coeffsNode, coeffGroupsRole, coeffRoleMapNode);
			} else {
				writeCsvForGroup(new File(path, "coeff-edges-scope"), coeffsEdge, coeffGroupsScope, coeffScopeMapEdge);
				writeCsvForGroup(new File(path, "coeff-edges-role"), coeffsEdge, coeffGroupsRole, coeffRoleMapEdge);
			}
		}
	}

	private static class CompletenessRecord {
		public final double fromCompleteness;
		public final double toCompleteness;
		public final double jsd;
		public final boolean inBestSubset;

		public CompletenessRecord(double fromCompleteness, double toCompleteness, double jsd, boolean inBestSubset) {
			this.fromCompleteness = fromCompleteness;
			this.toCompleteness = toCompleteness;
			this.jsd = jsd;
			this.inBestSubset = inBestSubset;
		}
	}

	private static class CompletenessAggregator implements ResultFilter<PredictionResult> {
		private final boolean nodes;
		private final Map<Long, Double> completeness;
		private final List<CompletenessRecord> data = new ArrayList<CompletenessRecord>();
		private int numInBestSubset = 0;

		public CompletenessAggregator(boolean nodes, Map<Long, Double> completeness) {
			this.nodes = nodes;
			this.completeness = completeness;
		}

		@Override
		public boolean evaluate(Date timePeriodBase, ModelVariant variant, PredictionResult result) {
			boolean matches = false;
			if (nodes == (variant.getEntity() == EntityType.NODE))
				matches = IndicatorWeightPair.inBestSubset(variant);

			Double completenessFrom = completeness.get(result.getFromUserId());
			Double completenessTo = completeness.get(result.getToUserId());
			if ((completenessFrom != null) && (completenessTo != null)) {
				data.add(new CompletenessRecord(completenessFrom, completenessTo, result.getDistJensenShannon(),
						matches));
				if (matches)
					numInBestSubset++;
			}
			return false;
		}

		@Override
		public void writeCsv(String path) throws IOException {
			int numRecords = data.size();
			double[][] dataAll = new double[numRecords][3];
			double[][] dataBest = new double[numInBestSubset][3];
			int idx = 0;
			int idxBest = 0;
			for (CompletenessRecord r : data) {
				dataAll[idx][0] = r.fromCompleteness;
				dataAll[idx][1] = r.toCompleteness;
				dataAll[idx][2] = r.jsd;
				if (r.inBestSubset) {
					dataBest[idxBest][0] = r.fromCompleteness;
					dataBest[idxBest][1] = r.toCompleteness;
					dataBest[idxBest][2] = r.jsd;
					idxBest++;
				}
				idx++;
			}

			PearsonsCorrelation pc = new PearsonsCorrelation(dataAll);
			RealMatrix corAll = pc.getCorrelationMatrix();
			RealMatrix sigAll = pc.getCorrelationPValues();
			RealMatrix corBest = null, sigBest = null;
			if (numInBestSubset > 0) {
				pc = new PearsonsCorrelation(dataBest);
				corBest = pc.getCorrelationMatrix();
				sigBest = pc.getCorrelationPValues();
			}

			PrintWriter w = new PrintWriter(new FileWriter(new File(path,
					"completeness-" + (nodes ? "nodes" : "edges") + ".txt")));
			try {
				if (numInBestSubset > 0) {
					w.println("best subset: from = " + corBest.getEntry(0, 2) + " (p = " + sigBest.getEntry(0, 2) +
							"), to = " + corBest.getEntry(1, 2) + " (p = " + sigBest.getEntry(1, 2) + ")");
				}
				w.println("all: from = " + corAll.getEntry(0, 2) + " (p = " + sigAll.getEntry(0, 2) +
						"), to = " + corAll.getEntry(1, 2) + " (p = " + sigAll.getEntry(1, 2) + ")");
			} finally {
				w.close();
			}
		}
	}

	private static class PerformanceRecord {
		public final StorelessUnivariateStatistic statJsd = new Mean();
		public final StorelessUnivariateStatistic statZeroOne = new Mean();
	}

	private static class BaselineAggregator implements ResultFilter<PredictionResult> {
		private final Map<ModelVariant.ModelType, PerformanceRecord> statsNodes =
				new HashMap<ModelVariant.ModelType, PerformanceRecord>();
		private final Map<ModelVariant.ModelType, PerformanceRecord> statsEdges =
				new HashMap<ModelVariant.ModelType, PerformanceRecord>();

		@Override
		public boolean evaluate(Date timePeriodBase, ModelVariant variant, PredictionResult result) {
			if (IndicatorWeightPair.inBestSubset(variant) ||
				(variant.getType() == ModelVariant.ModelType.BASELINE_PREVIOUS) ||
				(variant.getType() == ModelVariant.ModelType.BASELINE_UNIFORM) ||
				(variant.getType() == ModelVariant.ModelType.BASELINE_RANDOM)) {
				Map<ModelVariant.ModelType, PerformanceRecord> stats = (variant.getEntity() == EntityType.NODE) ?
						statsNodes : statsEdges;
				PerformanceRecord r = stats.get(variant.getType());
				if (r == null) {
					r = new PerformanceRecord();
					stats.put(variant.getType(), r);
				}
				r.statJsd.increment(result.getDistJensenShannon());
				r.statZeroOne.increment(result.getDistZeroOne());
			}
			return false;
		}

		@Override
		public void writeCsv(String path) throws IOException {
			PrintWriter w = new PrintWriter(new FileWriter(new File(path, "baseline.csv")));
			try {
				w.println("entity;type;mean jsd;mean zero-one");
				for (Map.Entry<ModelVariant.ModelType, PerformanceRecord> e : statsNodes.entrySet()) {
					w.println("NODE;" + e.getKey() + ";" + e.getValue().statJsd.getResult() + ";" +
							e.getValue().statZeroOne.getResult());
				}
				for (Map.Entry<ModelVariant.ModelType, PerformanceRecord> e : statsEdges.entrySet()) {
					w.println("COMMUNICATION;" + e.getKey() + ";" + e.getValue().statJsd.getResult() + ";" +
							e.getValue().statZeroOne.getResult());
				}
			} finally {
				w.close();
			}
		}
	}

	private static class TrainingErrorStatistics {
		private double sumPerUserErrorMean = 0.0;
		private double sumPerUserErrorVar = 0.0;
		private double sumErrorMean = 0.0;
		private double sumErrorVar = 0.0;
		private int numVariants = 0;

		public void add(ParameterFittingResult result) {
			sumPerUserErrorMean += result.getTrainErrorUserMean();
			sumPerUserErrorVar += Math.pow(result.getTrainErrorUserStdDev(), 2.0);
			sumErrorMean += result.getTrainErrorMean();
			sumErrorVar += Math.pow(result.getTrainErrorStdDev(), 2.0);
			numVariants++;
		}

		@Override
		public String toString() {
			return (sumPerUserErrorMean / numVariants) + ";" + (Math.sqrt(sumPerUserErrorVar) / numVariants) + ";" +
					(sumErrorMean / numVariants) + ";" + (Math.sqrt(sumErrorVar) / numVariants);
		}
	}

	private static class TrainingErrorAggregator implements ResultFilter<ParameterFittingResult> {
		private final TrainingErrorStatistics nodeStats = new TrainingErrorStatistics();
		private final TrainingErrorStatistics edgeStats = new TrainingErrorStatistics();

		@Override
		public boolean evaluate(Date timePeriodBase, ModelVariant variant, ParameterFittingResult result) {
			if (variant.getType() == ModelVariant.ModelType.SCIM) {
				if (variant.getEntity() == EntityType.NODE)
					nodeStats.add(result);
				else
					edgeStats.add(result);
			}
			return false;
		}

		@Override
		public void writeCsv(String path) throws IOException {
			PrintWriter w = new PrintWriter(new FileWriter(new File(path, "training-error.csv")));
			try {
				w.println("entity;per-user error mean;per-user error std.dev.;error mean;error std.dev.");
				w.println("NODE;" + nodeStats);
				w.println("COMMUNICATION;" + edgeStats);
			} finally {
				w.close();
			}
		}
	}

	private static class IndicatorWeightPair {
		private static IndicatorWeightPair[] bestPairsNode;
		private static IndicatorWeightPair[] bestPairsEdge;

		public final IndicatorFunction indicator;
		public final WeightFunction weight;

		public IndicatorWeightPair(IndicatorFunction indicator, WeightFunction weight) {
			this.indicator = indicator;
			this.weight = weight;
		}

		public boolean matches(ModelVariant variant) {
			return ((variant.getIndicator() == indicator) && (variant.getWeight() == weight));
		}

		private static IndicatorWeightPair[] readBestSubset(String fileName) throws IOException {
			List<IndicatorWeightPair> bestPairs = new ArrayList<IndicatorWeightPair>();
			BufferedReader r = new BufferedReader(new FileReader(fileName));
			try {
				String line = r.readLine();	// skip header
				line = r.readLine();	// first row contains the members of the best subset
				String[] parts = line.split(";");
				for (int i = 4; i < parts.length; i += 2) {
					String part = parts[i];
					if (part.isEmpty() || part.startsWith("BASELINE_"))
						continue;
					// assume that no indicator function name is a prefix of another
					boolean found = false;
					for (IndicatorFunction indicator : IndicatorFunction.values()) {
						if (part.startsWith(indicator.name())) {
							part = part.substring(indicator.name().length() + 1);
							WeightFunction weight = WeightFunction.valueOf(part);
							bestPairs.add(new IndicatorWeightPair(indicator, weight));
							found = true;
							break;
						}
					}
					if (!found)
						System.err.println("could not parse '" + parts[i] + "'");
				}
			} finally {
				r.close();
			}
			return bestPairs.toArray(new IndicatorWeightPair[bestPairs.size()]);
		}

		public static void loadBestSubsets(String fileNameNode, String fileNameEdge) throws IOException {
			bestPairsNode = readBestSubset(fileNameNode);
			bestPairsEdge = readBestSubset(fileNameEdge);
		}

		public static boolean inBestSubset(ModelVariant variant) {
			IndicatorWeightPair[] bestPairs = (variant.getEntity() == EntityType.NODE) ? bestPairsNode : bestPairsEdge;
			for (IndicatorWeightPair p : bestPairs) {
				if (p.matches(variant))
					return true;
			}
			return false;
		}
	}

	private static Map<Long, Double> loadCompleteness(File f) throws IOException {
		Map<Long, Double> completeness = new HashMap<Long, Double>();
		BufferedReader r = new BufferedReader(new FileReader(f));
		try {
			String line;
			while ((line = r.readLine()) != null) {
				String[] parts = line.split(";");
				completeness.put(Long.parseLong(parts[0]), Double.parseDouble(parts[1]));
			}
		} finally {
			r.close();
		}
		return completeness;
	}

	// Operates on output of ComputeCompleteness, HomogeneousSubsets, PrepareAnova and InfluenceExperimentsAnova. The
	// output of HomogeneousSubsets requires manual post-processing, so the file names are specified on the command
	// line. Take the first part of the output of HomogeneousSubsets ("homogeneous groups"), sort by asc. mean, save as
	// CSV.
	public static void main(String[] args) throws Exception {
		ExperimentConfiguration cfg = new ExperimentConfiguration(InfluenceExperimentsAnova.class);
		String path = cfg.getLocalProperty(InfluenceExperimentsAnova.PROP_OUT_PATH, ".");

		if (args.length < 2) {
			System.err.println("usage: " + AnalyzeBestSubset.class.getName() + " <subsets nodes> <subset edges>");
			return;
		}
		IndicatorWeightPair.loadBestSubsets(args[0], args[1]);

		Map<Long, Double> completeness = loadCompleteness(new File(path, ComputeCompleteness.FILE_COMPLETENESS));

		Set<ResultFilter<?>> aggregators = new HashSet<ResultFilter<?>>();
		BaselineAggregator baselineAggregator = new BaselineAggregator();

		List<ResultFilter<PredictionResult>> aggregatorsNode = new ArrayList<ResultFilter<PredictionResult>>();
		aggregatorsNode.add(new CompletenessAggregator(true, completeness));
		aggregatorsNode.add(baselineAggregator);
		aggregators.addAll(aggregatorsNode);
		ExperimentSeries.readPredictionResultsCsv(new File(path, PrepareAnova.FILE_EVAL_ANOVA_NODES), aggregatorsNode);

		List<ResultFilter<PredictionResult>> aggregatorsEdge = new ArrayList<ResultFilter<PredictionResult>>();
		aggregatorsEdge.add(new CompletenessAggregator(false, completeness));
		aggregatorsEdge.add(baselineAggregator);
		aggregators.addAll(aggregatorsEdge);
		ExperimentSeries.readPredictionResultsCsv(new File(path, PrepareAnova.FILE_EVAL_ANOVA_EDGES), aggregatorsEdge);

		List<ResultFilter<ParameterFittingResult>> aggregatorsParam =
				new ArrayList<ResultFilter<ParameterFittingResult>>();
		aggregatorsParam.add(new SankeyAggregator(true));
		aggregatorsParam.add(new SankeyAggregator(false));
		aggregatorsParam.add(new TrainingErrorAggregator());
		aggregators.addAll(aggregatorsParam);
		ExperimentSeries.readParameterFittingResultsCsv(new File(path, InfluenceExperimentsAnova.FILE_EVAL_PARAMETERS),
				aggregatorsParam);

		for (ResultFilter<?> aggregator : aggregators)
			aggregator.writeCsv(path);
	}

}
