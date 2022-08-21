package edu.tum.cs.util.arrays;

import java.io.Serializable;

/**
 * Interface for accessing 2D integer arrays
 */
public interface Int2DArray extends Serializable {

	public int getNumRows();

	public int getNumColumns();

	public int internalSize();

	public int get(int row, int column);

	public void set(int row, int column, int value);

	public int adjust(int row, int column, int offset);

	/**
	 * @return an iterator over the elements explicitly stored in the underlying data structure; iteration order is
	 * 	undefined
	 */
	public Sparse2DIterator<Integer> sparseIterator();

}
