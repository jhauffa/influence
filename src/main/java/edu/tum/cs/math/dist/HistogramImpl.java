package edu.tum.cs.math.dist;

public class HistogramImpl {

	public static class IdentityHistogram extends Histogram<Integer> {
		private final int maxValue;

		public IdentityHistogram() {
			this(Integer.MAX_VALUE);
		}

		public IdentityHistogram(int maxValue) {
			this.maxValue = maxValue;
		}

		@Override
		protected int assignBin(Integer value) {
			if (value > maxValue)
				return Integer.MAX_VALUE;
			return (int) value;
		}

		@Override
		protected Integer getLowerLimit(int binIdx) {
			if (binIdx == Integer.MAX_VALUE)
				return maxValue + 1;
			return binIdx;
		}

		@Override
		protected Integer getUpperLimit(int binIdx) {
			if (binIdx == Integer.MAX_VALUE)
				return Integer.MAX_VALUE;
			return binIdx;
		}
	}

	public static class QuantizedProbabilityHistogram extends Histogram<Double> {
		private final double quantizationFactor;

		public QuantizedProbabilityHistogram(double quantizationFactor) {
			this.quantizationFactor = quantizationFactor;
		}

		@Override
		protected int assignBin(Double value) {
			if (value < 0.0)
				return -1;
			if (value > 1.0)
				value = 1.0;
			return (int) Math.ceil(value * quantizationFactor);
		}

		@Override
		protected Double getLowerLimit(int binIdx) {
			if (binIdx == -1)
				return Double.NEGATIVE_INFINITY;
			if (binIdx == 0)
				return 0.0;
			return ((double) binIdx - 1) / quantizationFactor;
		}

		@Override
		protected Double getUpperLimit(int binIdx) {
			if (binIdx == -1)
				return 0.0;
			if (binIdx == 0.0)
				return 0.0;
			return (double) binIdx / quantizationFactor;
		}
	}

	public static class Log10Histogram extends Histogram<Double> {
		@Override
		protected int assignBin(Double value) {
			if (value == 0.0)
				return Integer.MAX_VALUE;
			return (int) Math.ceil(Math.log10(value));
		}

		@Override
		protected Double getLowerLimit(int binIdx) {
			if (binIdx == Integer.MAX_VALUE)
				return 0.0;
			return Math.pow(10.0, binIdx - 1);
		}

		@Override
		protected Double getUpperLimit(int binIdx) {
			if (binIdx == Integer.MAX_VALUE)
				return 0.0;
			return Math.pow(10.0, binIdx);
		}
	}

}
