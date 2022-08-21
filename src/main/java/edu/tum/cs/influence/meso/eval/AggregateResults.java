package edu.tum.cs.influence.meso.eval;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.StorelessUnivariateStatistic;
import org.apache.commons.math3.stat.descriptive.moment.Mean;

import edu.tum.cs.db.SocialGraphDao;
import edu.tum.cs.influence.meso.EntityType;
import edu.tum.cs.influence.meso.ExperimentResult;
import edu.tum.cs.influence.meso.ExperimentSeries;
import edu.tum.cs.influence.meso.IndicatorFunction;
import edu.tum.cs.influence.meso.InfluenceExperimentsAnova;
import edu.tum.cs.influence.meso.ModelVariant;
import edu.tum.cs.influence.meso.ModelVariantExperiment;
import edu.tum.cs.influence.meso.PredictionResult;
import edu.tum.cs.influence.meso.ResultFilter;
import edu.tum.cs.influence.meso.WeightFunction;
import edu.tum.cs.influence.meso.ModelVariant.ModelType;
import edu.tum.cs.util.ExperimentConfiguration;

public class AggregateResults {

	private static class ExperimentAverageAggregator implements ResultFilter<PredictionResult> {
		private final Map<Date, Map<ModelVariant, List<StorelessUnivariateStatistic>>> stat =
				new HashMap<Date, Map<ModelVariant, List<StorelessUnivariateStatistic>>>();
		private final Map<Date, Map<ModelVariant, Double>> bestBaselines =
				new HashMap<Date, Map<ModelVariant, Double>>();

		private final boolean isExplicit;
		private final boolean isSymmetric;

		public ExperimentAverageAggregator(boolean isExplicit, boolean isSymmetric) {
			this.isExplicit = isExplicit;
			this.isSymmetric = isSymmetric;
		}

		@Override
		public boolean evaluate(Date timePeriodBase, ModelVariant variant, PredictionResult result) {
			// Aggregate all experiment results (prediction accuracy of individual users) for each combination of
			// experiment parameters: timePeriodBase + ModelVariant (timePeriodLength, entity, type, indicator, weight)
			Map<ModelVariant, List<StorelessUnivariateStatistic>> variants = stat.get(timePeriodBase);
			if (variants == null) {
				variants = new HashMap<ModelVariant, List<StorelessUnivariateStatistic>>();
				stat.put(timePeriodBase, variants);
			}

			List<StorelessUnivariateStatistic> variantStat = variants.get(variant);
			if (variantStat == null) {
				variantStat = new ArrayList<StorelessUnivariateStatistic>(3);
				for (int i = 0; i < 3; i++)
					variantStat.add(new Mean());
				variants.put(variant, variantStat);
			}

			variantStat.get(0).increment(result.getDistJensenShannon());
			variantStat.get(1).increment(result.getDistCosine());
			variantStat.get(2).increment(result.getDistZeroOne());
			return false;
		}

		private void registerBaseline(Date timePeriodBase, int timePeriodLength, EntityType entity,
				IndicatorFunction indicator, WeightFunction weight, double jsd) {
			Map<ModelVariant, Double> baselineMap = bestBaselines.get(timePeriodBase);
			if (baselineMap == null) {
				baselineMap = new HashMap<ModelVariant, Double>();
				bestBaselines.put(timePeriodBase, baselineMap);
			}
			ModelVariant variant = new ModelVariant(timePeriodLength, entity, ModelVariant.ModelType.SCIM,
					indicator, weight);
			Double curJsd = baselineMap.get(variant);
			if ((curJsd == null) || (jsd < curJsd))
				baselineMap.put(variant, jsd);
		}

		private double getBestBaseline(Date timePeriodBase, ModelVariant variant) {
			Double bestBaseline = null;
			Map<ModelVariant, Double> baselineMap = bestBaselines.get(timePeriodBase);
			if (baselineMap != null)
				bestBaseline = baselineMap.get(variant);
			return (bestBaseline != null) ? bestBaseline : Double.POSITIVE_INFINITY;
		}

		private void writeVariantsBelowBaseline(File fileVariants, File fileNeighborhoodsDiscarded,
				File fileNeighborhoodStats, ExperimentSeries avgSeries) throws IOException {
			int numRead = 0;
			int numDiscarded = 0;
			Set<ModelVariant> neighborhoods = new HashSet<ModelVariant>();
			PrintWriter w = new PrintWriter(new FileWriter(fileVariants));
			try {
				for (ModelVariantExperiment exp : avgSeries) {
					for (Map.Entry<ModelVariant, ExperimentResult> e : exp.getResults().entrySet()) {
						ModelVariant variant = e.getKey();
						PredictionResult result = e.getValue().iterator().next();
						if (variant.getType() == ModelVariant.ModelType.SCIM) {
							numRead++;
							if (getBestBaseline(exp.getTimePeriodBase(), variant) <= result.getDistJensenShannon()) {
								w.println(variant);
								neighborhoods.add(new ModelVariant(0, (variant.getEntity() == EntityType.NODE ?
										EntityType.NODE : EntityType.EDGE_COMMUNICATION), ModelVariant.ModelType.SCIM,
										variant.getIndicator(), variant.getWeight()));
								numDiscarded++;
							}
						}
					}
				}
			} finally {
				w.close();
			}
			System.out.println("discarded " + numDiscarded + " of " + numRead + " variants");

			Set<IndicatorFunction> indicators = IndicatorFunction.getApplicable(isExplicit, isSymmetric);
			Set<WeightFunction> weights = WeightFunction.getApplicable(isExplicit, isSymmetric);

			PrintWriter wDiscarded = new PrintWriter(new FileWriter(fileNeighborhoodsDiscarded));
			PrintWriter wStats = new PrintWriter(new FileWriter(fileNeighborhoodStats));
			try {
				Map<WeightFunction, Set<IndicatorFunction>> keptIndicatorsNodes =
						new HashMap<WeightFunction, Set<IndicatorFunction>>();
				Map<IndicatorFunction, Set<WeightFunction>> keptWeightsNodes =
						new HashMap<IndicatorFunction, Set<WeightFunction>>();
				Map<WeightFunction, Set<IndicatorFunction>> keptIndicatorsEdges =
						new HashMap<WeightFunction, Set<IndicatorFunction>>();
				Map<IndicatorFunction, Set<WeightFunction>> keptWeightsEdges =
						new HashMap<IndicatorFunction, Set<WeightFunction>>();
				for (IndicatorFunction indicator : indicators) {
					keptWeightsNodes.put(indicator, new HashSet<WeightFunction>());
					keptWeightsEdges.put(indicator, new HashSet<WeightFunction>());
				}
				for (WeightFunction weight : weights) {
					keptIndicatorsNodes.put(weight, new HashSet<IndicatorFunction>());
					keptIndicatorsEdges.put(weight, new HashSet<IndicatorFunction>());
				}

				wStats.println("kept neighborhoods:");
				for (EntityType entity : new EntityType[] { EntityType.NODE, EntityType.EDGE_COMMUNICATION }) {
					for (IndicatorFunction indicator : indicators) {
						for (WeightFunction weight : weights) {
							ModelVariant variant = new ModelVariant(0, entity, ModelVariant.ModelType.SCIM,
									indicator, weight);
							if (!neighborhoods.contains(variant)) {
								wStats.println(variant);
								if (entity == EntityType.NODE) {
									keptIndicatorsNodes.get(weight).add(indicator);
									keptWeightsNodes.get(indicator).add(weight);
								} else {
									keptIndicatorsEdges.get(weight).add(indicator);
									keptWeightsEdges.get(indicator).add(weight);
								}
							} else
								wDiscarded.println(variant);
						}
					}
				}

				wStats.println("survival rates:");
				for (Map.Entry<IndicatorFunction, Set<WeightFunction>> e : keptWeightsNodes.entrySet()) {
					wStats.println("NODE;" + e.getKey().name() + ";;" +
							((double) e.getValue().size() / weights.size()));
				}
				for (Map.Entry<WeightFunction, Set<IndicatorFunction>> e : keptIndicatorsNodes.entrySet()) {
					wStats.println("NODE;;" + e.getKey().name() + ";" +
							((double) e.getValue().size() / indicators.size()));
				}
				for (Map.Entry<IndicatorFunction, Set<WeightFunction>> e : keptWeightsEdges.entrySet()) {
					wStats.println("COMMUNICATION;" + e.getKey().name() + ";;" +
							((double) e.getValue().size() / weights.size()));
				}
				for (Map.Entry<WeightFunction, Set<IndicatorFunction>> e : keptIndicatorsEdges.entrySet()) {
					wStats.println("COMMUNICATION;;" + e.getKey().name() + ";" +
							((double) e.getValue().size() / indicators.size()));
				}
			} finally {
				wStats.close();
				wDiscarded.close();
			}
			System.out.println("discarded " + neighborhoods.size() + " of " + (2 * indicators.size() * weights.size()) +
					" neighborhoods");
		}

		@Override
		public void writeCsv(String path) throws IOException {
			ExperimentSeries series = new ExperimentSeries();
			for (Map.Entry<Date, Map<ModelVariant, List<StorelessUnivariateStatistic>>> e1 : stat.entrySet()) {
				Date timePeriodBase = e1.getKey();
				ModelVariantExperiment exp = new ModelVariantExperiment(timePeriodBase);
				for (Map.Entry<ModelVariant, List<StorelessUnivariateStatistic>> e2 : e1.getValue().entrySet()) {
					List<StorelessUnivariateStatistic> variantStat = e2.getValue();
					PredictionResult result = new PredictionResult(-1, -1, variantStat.get(0).getResult(),
							variantStat.get(1).getResult(), variantStat.get(2).getResult());

					ModelVariant variant = e2.getKey();
					ExperimentResult resultContainer = new ExperimentResult(variant);
					resultContainer.addPredictionResult(result);
					exp.addResult(resultContainer);

					// update best baseline map
					if (variant.getType() == ModelType.CONST_COEFF) {
						registerBaseline(timePeriodBase, variant.getTimePeriodLength(), variant.getEntity(),
								variant.getIndicator(), variant.getWeight(), result.getDistJensenShannon());
					} else if ((variant.getType() == ModelType.BASELINE_PREVIOUS) ||
							   (variant.getType() == ModelType.BASELINE_RANDOM) ||
							   (variant.getType() == ModelType.BASELINE_UNIFORM) ||
							   (variant.getType() == ModelType.EMPTY_NEIGHBORHOOD)) {
						for (IndicatorFunction indicator : IndicatorFunction.getApplicable(isExplicit, isSymmetric)) {
							for (WeightFunction weight : WeightFunction.getApplicable(isExplicit, isSymmetric)) {
								registerBaseline(timePeriodBase, variant.getTimePeriodLength(), variant.getEntity(),
										indicator, weight, result.getDistJensenShannon());
							}
						}
					}
				}
				series.addExperiment(exp);
			}
			series.writePredictionResultsCsv(new File(path, "eval-prediction-aggregate.csv"));
			writeVariantsBelowBaseline(new File(path, "discarded-variants.csv"),
					new File(path, "discarded-neighborhoods.csv"), new File(path, "neighborhood-stats.csv"), series);
		}
	}

	public static void main(String[] args) throws Exception {
		ExperimentConfiguration cfg = new ExperimentConfiguration(InfluenceExperimentsAnova.class);
		String path = cfg.getLocalProperty(InfluenceExperimentsAnova.PROP_OUT_PATH, ".");

		boolean isExplicit, isSymmetric;
		if (args.length > 1) {
			isExplicit = Boolean.parseBoolean(args[0]);
			isSymmetric = Boolean.parseBoolean(args[1]);
		} else {
			SocialGraphDao dao = new SocialGraphDao();
			isExplicit = dao.hasExplicitEdges();
			isSymmetric = dao.isExplicitSymmetric();
		}

		List<ResultFilter<PredictionResult>> aggregators = new ArrayList<ResultFilter<PredictionResult>>();
		aggregators.add(new ExperimentAverageAggregator(isExplicit, isSymmetric));
		ExperimentSeries.readPredictionResultsCsv(new File(path, InfluenceExperimentsAnova.FILE_EVAL_PREDICTION),
				aggregators);
		aggregators.get(0).writeCsv(path);
	}

}
