package edu.tum.cs.math.dist;

public class HistogramBin<T> {

	public final int index;
	public final T lowerLimit, upperLimit;
	public final long count;
	public final double p;

	public HistogramBin(int index, T lowerLimit, T upperLimit, long count, double p) {
		this.index = index;
		this.lowerLimit = lowerLimit;
		this.upperLimit = upperLimit;
		this.count = count;
		this.p = p;
	}

}
