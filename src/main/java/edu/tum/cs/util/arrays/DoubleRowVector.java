package edu.tum.cs.util.arrays;

public class DoubleRowVector {

	protected final double[] v;
	protected final Double2DArray matrix;
	protected final int row;

	public DoubleRowVector(double[] v) {
		this.v = v;
		this.matrix = null;
		this.row = -1;
	}

	public DoubleRowVector(Double2DArray matrix, int row) {
		this.v = null;
		this.matrix = matrix;
		this.row = row;
	}

	public int size() {
		if (v != null)
			return v.length;
		return matrix.getNumColumns();
	}

	public double get(int i) {
		if (v != null)
			return v[i];
		return matrix.get(row, i);
	}

	public double[] toArray() {
		if (v != null)
			return v;

		int n = matrix.getNumColumns();
		double[] array = new double[n];
		for (int i = 0; i < n; i++)
			array[i] = matrix.get(row, i);
		return array;
	}

}
