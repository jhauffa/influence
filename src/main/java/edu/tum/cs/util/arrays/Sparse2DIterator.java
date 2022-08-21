package edu.tum.cs.util.arrays;

public interface Sparse2DIterator<T> {

	public boolean hasNext();

	public void advance();

	public T getValue();

	public int getRow();

	public int getColumn();

}
