package edu.tum.cs.util.arrays;

import java.io.Serializable;

import com.carrotsearch.hppc.cursors.IntObjectCursor;

import edu.tum.cs.util.io.SerializableIntObjectScatterMap;

public class SparseObject2DArray<T extends Object> implements Object2DArray<T>, SparseArray, Serializable {

	private static final long serialVersionUID = 3862060836455635520L;

	private class Iterator implements Sparse2DIterator<T> {
		private final java.util.Iterator<IntObjectCursor<T>> it = data.iterator();
		private IntObjectCursor<T> c;

		@Override
		public boolean hasNext() {
			return it.hasNext();
		}

		@Override
		public void advance() {
			c = it.next();
		}

		@Override
		public T getValue() {
			return c.value;
		}

		@Override
		public int getRow() {
			return (c.key >>> numColumnBits);
		}

		@Override
		public int getColumn() {
			return (c.key & ((1 << numColumnBits) - 1));
		}
	}

	protected static final int defaultInitialSize = 32;

	private static final int numColumnBits = 16;

	protected final int numRows;
	protected final int numColumns;

	protected final SerializableIntObjectScatterMap<T> data;

	public SparseObject2DArray(int numRows, int numColumns) {
		this(numRows, numColumns, defaultInitialSize);
	}

	public SparseObject2DArray(int numRows, int numColumns, int initialSize) {
		if (numColumns >= (1 << numColumnBits))
			throw new IllegalArgumentException("numColumns have to use only " + numColumnBits + " bits");
		if (numRows >= (1 << (32 - numColumnBits)))
			throw new IllegalArgumentException("numRows have to use only " + (32 - numColumnBits) + " bits");

		this.numRows = numRows;
		this.numColumns = numColumns;
		this.data = new SerializableIntObjectScatterMap<T>(initialSize, LOAD_FACTOR);
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
	public int size() {
		return data.size();
	}

	protected int getKey(int row, int column) {
		// assert row < getNumRows() && column < getNumColumns() && row >= 0 && column >= 0;
		return (row << numColumnBits) | column;
	}

	@Override
	public T get(int row, int column) {
		int index = data.indexOf(getKey(row, column));
		if (index >= 0)
			return data.indexGet(index);
		return null;
	}

	@Override
	public void set(int row, int column, T value) {
		int key = getKey(row, column);
		if (value != null)
			data.put(key, value);
		else
			data.remove(key);
	}

	@Override
	public Sparse2DIterator<T> sparseIterator() {
		return new Iterator();
	}

}
