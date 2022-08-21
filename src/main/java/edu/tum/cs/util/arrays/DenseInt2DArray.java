package edu.tum.cs.util.arrays;

/**
 * Boxing Class for int[][]
 */
public class DenseInt2DArray implements Int2DArray {

	private class Iterator implements Sparse2DIterator<Integer> {
		private int row, column = -1;

		@Override
		public boolean hasNext() {
			return ((column < (getNumColumns() - 1)) || (row < (getNumRows() - 1)));
		}

		@Override
		public void advance() {
			if (++column >= getNumColumns()) {
				column = 0;
				row++;
			}
		}

		@Override
		public Integer getValue() {
			return get(row, column);
		}

		@Override
		public int getRow() {
			return row;
		}

		@Override
		public int getColumn() {
			return column;
		}
	}

	private static final long serialVersionUID = 1760910246625236302L;

	protected final int[][] data;

	public DenseInt2DArray(int numRows, int numColumns) {
		data = new int[numRows][numColumns];
	}

	public DenseInt2DArray(int[][] data) {
		this.data = data;
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
	public int internalSize() {
		return data.length * data[0].length;
	}

	@Override
	public int get(int row, int column) {
		return data[row][column];
	}

	@Override
	public void set(int row, int column, int value) {
		data[row][column] = value;
	}

	@Override
	public int adjust(int row, int column, int offset) {
		return data[row][column] += offset;
	}

	@Override
	public Sparse2DIterator<Integer> sparseIterator() {
		return new Iterator();
	}

}
