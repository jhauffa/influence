package edu.tum.cs.util.arrays;

import java.io.Serializable;

/**
 * Interface for accessing 2D double arrays
 */
public interface Double2DArray extends Serializable {

	public int getNumRows();

	public int getNumColumns();

	public double get(int row, int column);

	public void set(int row, int column, double value);

	public double adjust(int row, int column, double offset);

}
