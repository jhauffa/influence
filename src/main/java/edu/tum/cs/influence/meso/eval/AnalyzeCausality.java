package edu.tum.cs.influence.meso.eval;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import edu.tum.cs.influence.meso.EntityType;
import edu.tum.cs.influence.meso.ExperimentSeries;
import edu.tum.cs.influence.meso.InfluenceExperimentsAnova;
import edu.tum.cs.influence.meso.ModelVariant;
import edu.tum.cs.influence.meso.ParameterFittingResult;
import edu.tum.cs.influence.meso.ResultFilter;
import edu.tum.cs.influence.meso.ModelVariant.ModelType;
import edu.tum.cs.math.dist.DirichletDistribution;
import edu.tum.cs.util.ExperimentConfiguration;

public class AnalyzeCausality {

	private static class CausalityAggregator implements ResultFilter<ParameterFittingResult> {
		private final boolean nodes;
		private final Map<ModelType, List<double[]>> coeff = new HashMap<ModelType, List<double[]>>();

		public CausalityAggregator(boolean nodes) {
			this.nodes = nodes;
			coeff.put(ModelType.SCIM, new ArrayList<double[]>());
			coeff.put(ModelType.SCIM_REVERSE, new ArrayList<double[]>());
			coeff.put(ModelType.SCIM_SHUFFLE, new ArrayList<double[]>());
		}

		@Override
		public boolean evaluate(Date timePeriodBase, ModelVariant variant, ParameterFittingResult result) {
			if (nodes == (variant.getEntity() == EntityType.NODE)) {
				ModelType type = variant.getType();
				if ((type == ModelType.SCIM) || (type == ModelType.SCIM_REVERSE) || (type == ModelType.SCIM_SHUFFLE)) {
					double[] modelParams = result.getModel().getCoefficients();
					int n = (nodes ? 10 : 14) + 1;
					double[] actualCoefficients = new double[n];
					double sum = 0.0;
					for (int i = 0; i < (n - 1); i++) {
						actualCoefficients[i] = modelParams[i];
						sum += modelParams[i];
					}
					actualCoefficients[n - 1] = 1.0 - sum;
					coeff.get(type).add(actualCoefficients);
				}
			}
			return false;
		}

		@Override
		public void writeCsv(String path) throws IOException {
			// perform likelihood-ratio test for "regular vs. reversed" and "regular vs. shuffled"
			List<double[]> coeffList = coeff.get(ModelType.SCIM);
			double[][] dataRef = coeffList.toArray(new double[coeffList.size()][]);

			PrintWriter w = new PrintWriter(new FileWriter(
					new File(path, "causality-" + (nodes ? "nodes" : "edges") + ".csv")));
			try {
				w.print("type;test statistic;p");
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < dataRef[0].length; j++)
						w.print(";a" + i + "_" + j);
				w.println();

				for (Map.Entry<ModelType, List<double[]>> e : coeff.entrySet()) {
					if (e.getKey() == ModelType.SCIM)
						continue;
					coeffList = e.getValue();
					double[][] data = coeffList.toArray(new double[coeffList.size()][]);

					w.print(e.getKey());
					for (double v : dirichletLikelihoodRatioTest(dataRef, data))
						w.print(";" + v);
					w.println();
				}
			} finally {
				w.close();
			}
		}

		/**
		 * The likelihood-ratio test compares how two models fit a particular sample. To use this as a two-sample test
		 * of homogeneity, we compare a joint model consisting of a separate model for each sample to a single model fit
		 * to the union of both samples. In this case our model is the Dirichlet distribution.
		 * High test statistic, low p: can choose to reject the null hypothesis of homogeneity.
		 * (see https://en.wikipedia.org/wiki/Likelihood-ratio_test and https://en.wikipedia.org/wiki/Wilks%27_theorem)
		 */
		private static double[] dirichletLikelihoodRatioTest(double[][] s1, double[][] s2) {
			double[][] sBoth = new double[s1.length + s2.length][];
			System.arraycopy(s1, 0, sBoth, 0, s1.length);
			System.arraycopy(s2, 0, sBoth, s1.length, s2.length);

			int n = sBoth[0].length;
			double[] stats = new double[2 + (3 * n)];
			DirichletDistribution dist1 = DirichletDistribution.fit(s1);
			System.arraycopy(dist1.getAlpha(), 0, stats, 2, n);
			DirichletDistribution dist2 = DirichletDistribution.fit(s2);
			System.arraycopy(dist2.getAlpha(), 0, stats, 2 + n, n);
			DirichletDistribution distBoth = DirichletDistribution.fit(sBoth);
			System.arraycopy(distBoth.getAlpha(), 0, stats, 2 + (2 * n), n);

			stats[0] = 2.0 * (dist1.logLikelihood(s1) + dist2.logLikelihood(s2) - distBoth.logLikelihood(sBoth));
			ChiSquaredDistribution distTest = new ChiSquaredDistribution(n);
			stats[1] = 1.0 - distTest.cumulativeProbability(stats[0]);
			return stats;
		}
	}

	public static void main(String[] args) throws Exception {
		ExperimentConfiguration cfg = new ExperimentConfiguration(InfluenceExperimentsAnova.class);
		String path = cfg.getLocalProperty(InfluenceExperimentsAnova.PROP_OUT_PATH, ".");

		List<ResultFilter<ParameterFittingResult>> aggregatorsParam =
				new ArrayList<ResultFilter<ParameterFittingResult>>();
		aggregatorsParam.add(new CausalityAggregator(true));
		aggregatorsParam.add(new CausalityAggregator(false));
		ExperimentSeries.readParameterFittingResultsCsv(new File(path, InfluenceExperimentsAnova.FILE_EVAL_PARAMETERS),
				aggregatorsParam);

		for (ResultFilter<?> aggregator : aggregatorsParam)
			aggregator.writeCsv(path);
	}

}
