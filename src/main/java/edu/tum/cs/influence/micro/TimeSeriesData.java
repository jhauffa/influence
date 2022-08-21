package edu.tum.cs.influence.micro;

import java.io.Serializable;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.rng.UniformRandomProvider;

public class TimeSeriesData implements Serializable {

	private static final long serialVersionUID = 6984271290780118656L;

	protected final double[][] values;
	private int numMissing;
	private int curIdx = 0;

	public TimeSeriesData(int maxValues) {
		values = new double[maxValues][];
		numMissing = maxValues;
	}

	public TimeSeriesData(TimeSeriesData other) {
		values = other.values.clone();
		numMissing = other.numMissing;
		curIdx = other.curIdx;
	}

	public int addValue(double[] v) {
		if (v != null) {
			values[curIdx] = v;
			numMissing--;
		}
		return ++curIdx;
	}

	public double[] getValue(int t) {
		return values[t];
	}

	public void setValue(int t, double[] v) {
		if ((values[t] == null) && (v != null))
			numMissing--;
		else if ((values[t] != null) && (v == null))
			numMissing++;
		values[t] = v;
	}

	public int getLength() {
		return values.length;
	}

	public boolean isMissing(int t) {
		return (values[t] == null);
	}

	public int getNumMissing() {
		return numMissing;
	}

	public void setMissingValues(double[] v) {
		for (int i = 0; i < values.length; i++) {
			if (values[i] == null)
				values[i] = v.clone();
		}
		numMissing = 0;
	}

	public double[][] toRowMajorArray() {
		if (numMissing > 0)
			throw new RuntimeException("cannot flatten time series with missing values");

		// copy and transpose values
		double[][] a = new double[values[0].length][values.length];
		for (int i = 0; i < values.length; i++) {
			for (int j = 0; j < values[i].length; j++)
				a[j][i] = values[i][j];
		}
		return a;
	}

	public int[] getIndices() {
		int[] indices = new int[values.length - numMissing];
		int pos = 0;
		for (int i = 0; i < values.length; i++) {
			if (values[i] != null)
				indices[pos++] = i;
			if (pos >= indices.length)
				break;
		}
		return indices;
	}

	public void permute(int[] indices, UniformRandomProvider rand) {
		for (int i = indices.length - 1; i >= 1; i--) {
			int j = rand.nextInt(i + 1);
			double[] tmp = values[indices[i]];
			values[indices[i]] = values[indices[j]];
			values[indices[j]] = tmp;
		}
	}

	public void permute(int[] indices, int permIdx, TimeSeriesData dst) {
		if (indices.length < 2)
			return;
		List<double[]> remainingValues = new LinkedList<double[]>();
		for (int idx : indices)
			remainingValues.add(values[idx]);

		int f = (int) CombinatoricsUtils.factorial(indices.length - 1);
		for (int i = 0; i < indices.length; i++) {
			int q = permIdx / f;
			permIdx %= f;
			dst.values[indices[i]] = remainingValues.remove(q);
			if (i < (indices.length - 1))
				f /= indices.length - (i + 1);
		}
	}

}
