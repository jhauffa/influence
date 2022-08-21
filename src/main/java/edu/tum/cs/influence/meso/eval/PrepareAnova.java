package edu.tum.cs.influence.meso.eval;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

import edu.tum.cs.influence.meso.EntityType;
import edu.tum.cs.influence.meso.ExperimentResult;
import edu.tum.cs.influence.meso.ExperimentSeries;
import edu.tum.cs.influence.meso.InfluenceExperimentsAnova;
import edu.tum.cs.influence.meso.ModelVariant;
import edu.tum.cs.influence.meso.ModelVariantExperiment;
import edu.tum.cs.influence.meso.PredictionResult;
import edu.tum.cs.influence.meso.ResultFilter;
import edu.tum.cs.influence.meso.ModelVariant.ModelType;
import edu.tum.cs.util.ExperimentConfiguration;

public class PrepareAnova {

	public static final String FILE_EVAL_ANOVA_NODES = "eval-anova-nodes.csv";
	public static final String FILE_EVAL_ANOVA_EDGES = "eval-anova-edges.csv";

	private static class AnovaAggregator implements ResultFilter<PredictionResult> {
		private final Set<ModelVariant> filteredVariants;
		private final ModelType targetModelType;
		private final boolean nodes;
		private final Map<Date, ModelVariantExperiment> experimentContainer =
				new HashMap<Date, ModelVariantExperiment>();

		public AnovaAggregator(Set<ModelVariant> filteredVariants, ModelType targetModelType, boolean nodes) {
			this.filteredVariants = filteredVariants;
			this.targetModelType = targetModelType;
			this.nodes = nodes;
		}

		private boolean isFiltered(ModelVariant variant) {
			if (filteredVariants == null)
				return false;
			return (filteredVariants.contains(variant) ||
					filteredVariants.contains(new ModelVariant(0, (variant.getEntity() == EntityType.NODE ?
							EntityType.NODE : EntityType.EDGE_COMMUNICATION), ModelVariant.ModelType.SCIM,
							variant.getIndicator(), variant.getWeight())));
		}

		@Override
		public boolean evaluate(Date timePeriodBase, ModelVariant variant, PredictionResult result) {
			if (((variant.getType() == targetModelType) ||
				 (variant.getType() == ModelType.BASELINE_PREVIOUS) ||
				 (variant.getType() == ModelType.BASELINE_RANDOM) ||
				 (variant.getType() == ModelType.BASELINE_UNIFORM)) &&
				(nodes == (variant.getEntity() == EntityType.NODE))) {
				// map to SCIM if necessary
				if (variant.getType() == targetModelType)
					variant.setType(ModelType.SCIM);

				// filter out variants that perform worse than the best baseline method
				if (isFiltered(variant))
					return false;

				// save remaining results
				ModelVariantExperiment exp = experimentContainer.get(timePeriodBase);
				if (exp == null) {
					exp = new ModelVariantExperiment(timePeriodBase);
					experimentContainer.put(timePeriodBase, exp);
				}
				ExperimentResult resultContainer = exp.getResult(variant);
				if (resultContainer == null) {
					resultContainer = new ExperimentResult(variant);
					exp.addResult(resultContainer);
				}
				resultContainer.addPredictionResult(result);
			}
			return false;
		}

		@Override
		public void writeCsv(String path) throws IOException {
			ExperimentSeries filteredSeries = new ExperimentSeries();
			for (ModelVariantExperiment exp : experimentContainer.values())
				filteredSeries.addExperiment(exp);
			filteredSeries.writePredictionResultsCsv(
					new File(path, nodes ? FILE_EVAL_ANOVA_NODES : FILE_EVAL_ANOVA_EDGES));
		}
	}

	private static Set<ModelVariant> readModelVariants(File f) throws IOException {
		Set<ModelVariant> variants = new HashSet<ModelVariant>();
		BufferedReader r = new BufferedReader(new FileReader(f));
		try {
			String line;
			while ((line = r.readLine()) != null) {
				String[] parts = line.split(";");
				Queue<String> remainingParts = new LinkedList<String>(Arrays.asList(parts));
				variants.add(ModelVariant.readCsv(remainingParts));
			}
		} finally {
			r.close();
		}
		return variants;
	}

	public static void main(String[] args) throws Exception {
		ExperimentConfiguration cfg = new ExperimentConfiguration(InfluenceExperimentsAnova.class);
		String path = cfg.getLocalProperty(InfluenceExperimentsAnova.PROP_OUT_PATH, ".");

		ModelType targetModelType = ModelType.SCIM;
		Set<ModelVariant> filteredVariants = null;
		if (args.length > 0) {
			File filterList = new File(args[0]);
			filteredVariants = readModelVariants(filterList);

			if (args.length > 1)
				targetModelType = ModelType.valueOf(args[1]);
		}

		List<ResultFilter<PredictionResult>> aggregators = new ArrayList<ResultFilter<PredictionResult>>();
		aggregators.add(new AnovaAggregator(filteredVariants, targetModelType, true));
		aggregators.add(new AnovaAggregator(filteredVariants, targetModelType, false));
		ExperimentSeries.readPredictionResultsCsv(new File(path, InfluenceExperimentsAnova.FILE_EVAL_PREDICTION),
				aggregators);
		for (ResultFilter<PredictionResult> aggregator : aggregators)
			aggregator.writeCsv(path);
	}

}
