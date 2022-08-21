package edu.tum.cs.influence.meso;

import java.util.Queue;

public class PredictionResult {

	private final long fromUserId;
	private final long toUserId;	// == fromUserId for nodes
	private final double distJensenShannon;
	private final double distCosine;
	private final double distZeroOne;

	public PredictionResult(long fromUserId, long toUserId, double distJensenShannon, double distCosine,
			double distZeroOne) {
		this.fromUserId = fromUserId;
		this.toUserId = toUserId;
		this.distJensenShannon = distJensenShannon;
		this.distCosine = distCosine;
		this.distZeroOne = distZeroOne;
	}

	public long getFromUserId() {
		return fromUserId;
	}

	public long getToUserId() {
		return toUserId;
	}

	public double getDistJensenShannon() {
		return distJensenShannon;
	}

	public double getDistCosine() {
		return distCosine;
	}

	public double getDistZeroOne() {
		return distZeroOne;
	}

	@Override
	public String toString() {
		return fromUserId + ";" + toUserId + ";" + distJensenShannon + ";" + distCosine + ";" + distZeroOne;
	}

	public static String getCsvHeader() {
		return "fromUserId;toUserId;distJensenShannon;distCosine;distZeroOne";
	}

	public static PredictionResult readCsv(Queue<String> parts) {
		long fromUserId = Long.parseLong(parts.poll());
		long toUserId = Long.parseLong(parts.poll());
		double distJensenShannon = Double.parseDouble(parts.poll());
		double distCosine = Double.parseDouble(parts.poll());
		double distZeroOne = Double.parseDouble(parts.poll());
		return new PredictionResult(fromUserId, toUserId, distJensenShannon, distCosine, distZeroOne);
	}

}
