package edu.tum.cs.influence.meso;

import java.io.PrintWriter;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Queue;

public class ModelVariantExperiment {

	private static final SimpleDateFormat df = new SimpleDateFormat("yyyy.MM.dd");

	private final Date timePeriodBase;
	private final Map<ModelVariant, ExperimentResult> variantResult;

	public ModelVariantExperiment(Date timePeriodBase) {
		this.timePeriodBase = timePeriodBase;
		this.variantResult = new HashMap<ModelVariant, ExperimentResult>();
	}

	public Date getTimePeriodBase() {
		return timePeriodBase;
	}

	public void addResult(ExperimentResult result) {
		variantResult.put(result.getModelVariant(), result);
	}

	public ExperimentResult getResult(ModelVariant variant) {
		return variantResult.get(variant);
	}

	public Map<ModelVariant, ExperimentResult> getResults() {
		return variantResult;
	}

	public void merge(ModelVariantExperiment other) {
		for (Map.Entry<ModelVariant, ExperimentResult> e : other.variantResult.entrySet()) {
			ExperimentResult result = variantResult.get(e.getKey());
			if (result == null)
				variantResult.put(e.getKey(), e.getValue());
			else
				result.merge(e.getValue());
		}
	}

	public static void writePredictionCsvHeader(PrintWriter w) {
		w.println("timePeriodBase;" + ModelVariant.getCsvHeader() + ";" + PredictionResult.getCsvHeader());
	}

	public void writePredictionCsv(PrintWriter w) {
		for (Map.Entry<ModelVariant, ExperimentResult> e : variantResult.entrySet()) {
			String variantColumns = e.getKey().toString();
			for (PredictionResult result : e.getValue())
				w.println(df.format(timePeriodBase) + ";" + variantColumns + ";" + result);
		}
	}

	public static ModelVariantExperiment readPredictionCsv(Queue<String> parts,
			List<? extends ResultFilter<PredictionResult>> filters) throws ParseException {
		Date timePeriodBase = df.parse(parts.poll());
		ModelVariant variant = ModelVariant.readCsv(parts);
		PredictionResult predictionResult = PredictionResult.readCsv(parts);

		ModelVariantExperiment exp = new ModelVariantExperiment(timePeriodBase);
		boolean keepResult = false;
		if (filters != null) {
			for (ResultFilter<PredictionResult> filter : filters)
				keepResult |= filter.evaluate(timePeriodBase, variant, predictionResult);
		}
		if (keepResult) {
			ExperimentResult result = new ExperimentResult(variant);
			result.addPredictionResult(predictionResult);
			exp.addResult(result);
		}
		return exp;
	}

	public static void writeParameterFittingCsvHeader(PrintWriter w) {
		w.println("timePeriodBase;" + ModelVariant.getCsvHeader() + ";" + ParameterFittingResult.getCsvHeader());
	}

	public void writeParameterFittingCsv(PrintWriter w) {
		for (Map.Entry<ModelVariant, ExperimentResult> e : variantResult.entrySet()) {
			ParameterFittingResult result = e.getValue().getParameterFittingResult();
			if (result != null)
				w.println(df.format(timePeriodBase) + ";" + e.getKey() + ";" + result);
		}
	}

	public static ModelVariantExperiment readParameterFittingCsv(Queue<String> parts,
			List<? extends ResultFilter<ParameterFittingResult>> filters) throws ParseException {
		Date timePeriodBase = df.parse(parts.poll());
		ModelVariant variant = ModelVariant.readCsv(parts);
		ParameterFittingResult parameterFittingResult = ParameterFittingResult.readCsv(parts);

		ModelVariantExperiment exp = new ModelVariantExperiment(timePeriodBase);
		boolean keepResult = false;
		if (filters != null) {
			for (ResultFilter<ParameterFittingResult> filter : filters)
				keepResult |= filter.evaluate(timePeriodBase, variant, parameterFittingResult);
		}
		if (keepResult) {
			ExperimentResult result = new ExperimentResult(variant);
			result.setParameterFittingResult(parameterFittingResult);
			exp.addResult(result);
		}
		return exp;
	}

}
