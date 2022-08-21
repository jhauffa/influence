package edu.tum.cs.influence.micro;

import edu.tum.cs.math.dist.DiscreteDistribution;

public class ChangeTimeSeriesData extends TimeSeriesData {

	private static final long serialVersionUID = 1879791318601849794L;

	private final double[] magnitudes;

	public ChangeTimeSeriesData(int maxValues) {
		super(maxValues);
		magnitudes = new double[maxValues];
	}

	public int addValues(double[] change, double magnitude) {
		int idx = super.addValue(change);
		magnitudes[idx - 1] = magnitude;
		return idx;
	}

	public double getMagnitude(int t) {
		return magnitudes[t];
	}

	public double computeSimilarity(TimeSeriesData srcData, int minActivePeriods, Explanation explanation) {
		double similarity = 0.0;
		int numActivePeriods = 0;

		double maxSimilarity = 0.0;
		int maxSimilarityIdx = 0;
		for (int i = 1; i < srcData.getLength(); i++) {
			if (!srcData.isMissing(i - 1) && !isMissing(i)) {
				// compute similarity for current time interval
				double curSimilarity = getMagnitude(i) * Math.max(0.0, (DiscreteDistribution.LOG2 -
						DiscreteDistribution.distJSe(srcData.getValue(i - 1), getValue(i))));
				if (curSimilarity >= maxSimilarity) {
					maxSimilarity = curSimilarity;
					maxSimilarityIdx = i;
				}

				// update compound similarity (no temporal decay: uniform weight for each time step)
				similarity += curSimilarity;
				numActivePeriods++;
			}
		}

		if (numActivePeriods < minActivePeriods)
			similarity = Double.NaN;
		else if (numActivePeriods == 0)
			similarity = 0.0;
		else
			similarity = (similarity / DiscreteDistribution.LOG2) / (srcData.getLength() - 1);

		if (explanation != null) {
			if (numActivePeriods > 0) {
				explanation.timePeriodIdx = maxSimilarityIdx;
				explanation.topicsInfluencer = srcData.getValue(maxSimilarityIdx - 1);
				explanation.topicsInfluencee = getValue(maxSimilarityIdx);
				explanation.localMagnitude = maxSimilarity;
				explanation.globalMagnitude = similarity;
			} else
				explanation.haveData = false;
		}
		return similarity;
	}

	public static ChangeTimeSeriesData transformData(TimeSeriesData data, TimeSeriesData syntheticData) {
		int n = data.getLength();
		ChangeTimeSeriesData tfData = new ChangeTimeSeriesData(n);
		if (n > 0) {
			tfData.addValues(null, 1.0);

			if (n > 1) {
				for (int i = 1; i < n; i++) {
					double[] t1 = data.getValue(i - 1);
					double[] t2;
					if (syntheticData != null)
						t2 = syntheticData.getValue(i);
					else
						t2 = data.getValue(i);

					if ((t1 != null) && (t2 != null)) {
						ChangeMeasurement change = ChangeMeasurement.computeChange(t1, t2);
						tfData.addValues(change.direction, change.magnitude);
					} else
						tfData.addValues(null, 1.0);
				}
			}
		}
		return tfData;
	}

}
