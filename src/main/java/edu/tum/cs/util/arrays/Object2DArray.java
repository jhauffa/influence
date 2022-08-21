package edu.tum.cs.util.arrays;

/**
 * Interface for accessing 2D generic object arrays
 */
public interface Object2DArray<T> {

	public int getNumRows();

	public int getNumColumns();

	public int size();

	public T get(int row, int column);

	public void set(int row, int column, T value);

	/**
	 * @return an iterator over the elements explicitly stored in the underlying data structure; iteration order is
	 * 	undefined
	 */
	public Sparse2DIterator<T> sparseIterator();

}
