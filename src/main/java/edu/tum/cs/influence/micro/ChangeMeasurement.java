package edu.tum.cs.influence.micro;

import java.util.Arrays;

public class ChangeMeasurement {

	public double magnitude;
	public final double[] direction;

	public ChangeMeasurement(double[] direction, double magnitude) {
		this.direction = direction;
		this.magnitude = magnitude;
	}

	public static ChangeMeasurement computeChange(double[] tPast, double[] tPresent) {
		double[] direction = new double[tPast.length];
		double magnitude = 0.0;
		double minComp = 1.0;
		for (int i = 0; i < tPast.length; i++) {
			direction[i] = Math.max(0.0, tPresent[i] - tPast[i]);
			magnitude += direction[i];
			if (tPast[i] < minComp)
				minComp = tPast[i];
		}
		if (magnitude > 0.0) {
			for (int i = 0; i < direction.length; i++)
				direction[i] /= magnitude;
		} else
			Arrays.fill(direction, 1.0 / direction.length);
		magnitude /= 1.0 - minComp;	// normalize to [0,1]
		return new ChangeMeasurement(direction, magnitude);
	}

}
