package edu.tum.cs.util.arrays;

public class FixedDouble2DArray implements Double2DArray {

	private static final long serialVersionUID = 891205681653687594L;
	private final int numRows;
	private final int numColumns;
	private final double value;

	public FixedDouble2DArray(int numRows, int numColumns, double value) {
		this.numRows = numRows;
		this.numColumns = numColumns;
		this.value = value;
	}

	@Override
	public int getNumRows() {
		return numRows;
	}

	@Override
	public int getNumColumns() {
		return numColumns;
	}

	@Override
	public double get(int row, int column) {
		return value;
	}

	@Override
	public void set(int row, int column, double value) {
		throw new UnsupportedOperationException("modifying a fixed array");
	}

	@Override
	public double adjust(int row, int column, double offset) {
		throw new UnsupportedOperationException("modifying a fixed array");
	}

}
