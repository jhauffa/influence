package edu.tum.cs.util.arrays;

import com.carrotsearch.hppc.cursors.IntIntCursor;

import edu.tum.cs.util.io.SerializableIntIntScatterMap;

public class SparseInt2DArray implements Int2DArray, SparseArray {

	private class Iterator implements Sparse2DIterator<Integer> {
		private final java.util.Iterator<IntIntCursor> it = data.iterator();
		private IntIntCursor c;

		@Override
		public boolean hasNext() {
			return it.hasNext();
		}

		@Override
		public void advance() {
			c = it.next();
		}

		@Override
		public Integer getValue() {
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

	private static final long serialVersionUID = -9219570292074840628L;

	protected static final int defaultInitialSize = 32;
	protected static final int defaultNoEntryValue = 0;

	private static final int numColumnBits = 8;	// columns == topics, max. 255 topics

	protected final int numRows;
	protected final int numColumns;

	protected final int noEntryValue;

	protected final SerializableIntIntScatterMap data;

	public SparseInt2DArray(int numRows, int numColumns) {
		this(numRows, numColumns, defaultNoEntryValue, defaultInitialSize);
	}

	public SparseInt2DArray(int numRows, int numColumns, int noEntryValue, int initialSize) {
		if (numColumns >= (1 << numColumnBits))
			throw new IllegalArgumentException("numColumns have to use only " + numColumnBits + " bits");
		if (numRows >= (1 << (32 - numColumnBits)))
			throw new IllegalArgumentException("numRows have to use only " + (32 - numColumnBits) + " bits");

		this.numRows = numRows;
		this.numColumns = numColumns;
		this.noEntryValue = noEntryValue;
		this.data = new SerializableIntIntScatterMap(initialSize, LOAD_FACTOR);
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
	public int internalSize() {
		return data.size();
	}

	protected int getKey(int row, int column) {
		// assert row < getNumRows() && column < getNumColumns() && row >= 0 && column >= 0;
		return (row << numColumnBits) | column;
	}

	@Override
	public int get(int row, int column) {
		int index = data.indexOf(getKey(row, column));
		if (index >= 0)
			return data.indexGet(index);
		return 0;
	}

	@Override
	public void set(int row, int column, int value) {
		int key = getKey(row, column);
		if (value != noEntryValue)
			data.put(key, value);
		else
			data.remove(key);
	}

	@Override
	public int adjust(int row, int column, int offset) {
		int key = getKey(row, column);
		int index = data.indexOf(key);
		if (index >= 0) {
			int value = data.indexGet(index) + offset;
			if (value != noEntryValue)
				data.indexReplace(index, value);
			else
				data.remove(key);
			return value;
		} else if (offset != 0) {
			data.put(key, offset);
			return offset;
		}
		return 0;
	}

	@Override
	public Sparse2DIterator<Integer> sparseIterator() {
		return new Iterator();
	}

}
