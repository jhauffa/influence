package edu.tum.cs.influence.meso;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class ExperimentResult implements Iterable<PredictionResult> {

	private final ModelVariant modelVariant;
	private final List<PredictionResult> predictionResults = new ArrayList<PredictionResult>();
	private ParameterFittingResult parameterFittingResult;

	public ExperimentResult(ModelVariant modelVariant) {
		this.modelVariant = modelVariant;
	}

	public ModelVariant getModelVariant() {
		return modelVariant;
	}

	public void addPredictionResult(PredictionResult result) {
		predictionResults.add(result);
	}

	public Iterator<PredictionResult> iterator() {
		return predictionResults.iterator();
	}

	public void setParameterFittingResult(ParameterFittingResult parameterFittingResult) {
		this.parameterFittingResult = parameterFittingResult;
	}

	public ParameterFittingResult getParameterFittingResult() {
		return parameterFittingResult;
	}

	public void merge(ExperimentResult other) {
		predictionResults.addAll(other.predictionResults);
	}

}
