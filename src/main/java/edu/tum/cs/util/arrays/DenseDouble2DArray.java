package edu.tum.cs.util.arrays;

/**
 * Boxing class for double[][]
 */
public class DenseDouble2DArray implements Double2DArray {

	private static final long serialVersionUID = -7028390418328112680L;

	protected final double[][] data;

	public DenseDouble2DArray(int numRows, int numColumns) {
		data = new double[numRows][numColumns];
	}

	/**
	 * Creates a deep copy of another double 2D array.
	 */
	public DenseDouble2DArray(Double2DArray other) {
		this(other.getNumRows(), other.getNumColumns());
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
				data[i][j] = other.get(i, j);
	}

	@Override
	public int getNumRows() {
		return data.length;
	}

	@Override
	public int getNumColumns() {
		return data[0].length;
	}

	@Override
	public double get(int row, int column) {
		return data[row][column];
	}

	@Override
	public void set(int row, int column, double value) {
		data[row][column] = value;
	}

	@Override
	public double adjust(int row, int column, double offset) {
		return data[row][column] += offset;
	}

}
