package edu.tum.cs.influence.meso;

import java.util.Queue;

public class ParameterFittingResult {

	private final PredictiveModel model;
	private final double trainErrorMean;
	private final double trainErrorStdDev;
	private final double trainErrorUserMean;
	private final double trainErrorUserStdDev;

	public ParameterFittingResult(PredictiveModel model, double trainErrorMean, double trainErrorStdDev,
			double trainErrorUserMean, double trainErrorUserStdDev) {
		this.model = model;
		this.trainErrorMean = trainErrorMean;
		this.trainErrorStdDev = trainErrorStdDev;
		this.trainErrorUserMean = trainErrorUserMean;
		this.trainErrorUserStdDev = trainErrorUserStdDev;
	}

	public PredictiveModel getModel() {
		return model;
	}

	public double getTrainErrorMean() {
		return trainErrorMean;
	}

	public double getTrainErrorStdDev() {
		return trainErrorStdDev;
	}

	public double getTrainErrorUserMean() {
		return trainErrorUserMean;
	}

	public double getTrainErrorUserStdDev() {
		return trainErrorUserStdDev;
	}

	public static String getCsvHeader() {
		StringBuilder sb = new StringBuilder(PredictiveModel.getCsvHeader());
		sb.append(";trainErrorMean;trainErrorStdDev;trainErrorUserMean;trainErrorUserStdDev");
		return sb.toString();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder(model.toString());
		sb.append(';').append(trainErrorMean).append(';').append(trainErrorStdDev);
		sb.append(';').append(trainErrorUserMean).append(';').append(trainErrorUserStdDev);
		return sb.toString();
	}

	public static ParameterFittingResult readCsv(Queue<String> parts) {
		PredictiveModel model = PredictiveModel.readCsv(parts);
		double trainErrorMean = Double.parseDouble(parts.poll());
		double trainErrorStdDev = Double.parseDouble(parts.poll());
		double trainErrorUserMean = Double.parseDouble(parts.poll());
		double trainErrorUserStdDev = Double.parseDouble(parts.poll());
		return new ParameterFittingResult(model, trainErrorMean, trainErrorStdDev,
				trainErrorUserMean, trainErrorUserStdDev);
	}

}
