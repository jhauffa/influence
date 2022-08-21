package edu.tum.cs.util.arrays;

public class SparseShadowInt2DArray extends SparseInt2DArray {

	private static final long serialVersionUID = 1L;

	private static final int defaultNoEntryValue = -1;

	private final Int2DArray source;

	public SparseShadowInt2DArray(Int2DArray source) {
		this(source, defaultNoEntryValue, defaultInitialSize);
	}

	public SparseShadowInt2DArray(Int2DArray source, int initialSize) {
		this(source, defaultNoEntryValue, initialSize);
	}

	public SparseShadowInt2DArray(Int2DArray source, int noEntryValue, int initialSize) {
		super(source.getNumRows(), source.getNumColumns(), noEntryValue, initialSize);
		this.source = source;
	}

	@Override
	public int get(int row, int column) {
		int key = getKey(row, column);
		int index = data.indexOf(key);
		if (index >= 0)
			return data.indexGet(index);
		return source.get(row, column);
	}

	@Override
	public int adjust(int row, int column, int offset) {
		int key = getKey(row, column);
		int index = data.indexOf(key);
		if (index >= 0) {
			int value = data.indexGet(index) + offset;
			data.indexReplace(index, value);
			return value;
		} else {
			int value = source.get(row, column) + offset;
			data.put(key, value);
			return value;
		}
	}

}
