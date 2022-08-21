package edu.tum.cs.util.arrays;

import java.util.Random;

import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class Int2DArrayTest {

	private static final long seed = 12345L;
	private static final int numRows = 500;
	private static final int numCols = 10;

	@Test
	public void testArrays() {
		Random rand = new Random(seed);

		Int2DArray denseArray = new DenseInt2DArray(numRows, numCols);
		Int2DArray sparseArray = new SparseInt2DArray(numRows, numCols);

		for (int i = 0; i < 1000; i++) {
			int row = rand.nextInt(numRows);
			int col = rand.nextInt(numCols);
			int value = rand.nextInt(10);

			denseArray.set(row, col, value);
			assertEquals(value, denseArray.get(row, col));
			sparseArray.set(row, col, value);
			assertEquals(value, sparseArray.get(row, col));

			int adjustValue = -1 * rand.nextInt(10);
			assertEquals(value + adjustValue, denseArray.adjust(row, col, adjustValue));
			assertEquals(value + adjustValue, denseArray.get(row, col));
			assertEquals(value + adjustValue, sparseArray.adjust(row, col, adjustValue));
			assertEquals(value + adjustValue, sparseArray.get(row, col));
		}

		for (int row = 0; row < numRows; row++)
			for (int col = 0; col < numCols; col++)
				assertEquals(row + "x" + col, denseArray.get(row, col), sparseArray.get(row, col));

		denseArray = new DenseInt2DArray(numRows, numCols);
		Int2DArray shadow = new SparseShadowInt2DArray(denseArray);
		Int2DArray denseEmptyArray = new DenseInt2DArray(numRows, numCols);
		for (int i = 0; i < 1000; i++) {
			int row = rand.nextInt(numRows);
			int col = rand.nextInt(numCols);
			int value = rand.nextInt(10);

			denseEmptyArray.set(row, col, value);
			assertEquals(value, denseEmptyArray.get(row, col));
			shadow.set(row, col, value);
			assertEquals(value, shadow.get(row, col));

			for (int j = 0; j < 100; j++) {
				int adjustValue = -1 * rand.nextInt(10);
				if ((value + adjustValue) == -1)
					adjustValue--;
				assertEquals("it" + j, value + adjustValue, denseEmptyArray.adjust(row, col, adjustValue));
				assertEquals("it" + j, value + adjustValue, denseEmptyArray.get(row, col));
				assertEquals("it" + j, value + adjustValue, shadow.adjust(row, col, adjustValue));
				assertEquals("it" + j, value + adjustValue, shadow.get(row, col));
				value += adjustValue;

				adjustValue = rand.nextInt(10);
				if ((value + adjustValue) == -1)
					adjustValue++;
				assertEquals("it" + j, value + adjustValue, denseEmptyArray.adjust(row, col, adjustValue));
				assertEquals("it" + j, value + adjustValue, denseEmptyArray.get(row, col));
				assertEquals("it" + j, value + adjustValue, shadow.adjust(row, col, adjustValue));
				assertEquals("it" + j, value + adjustValue, shadow.get(row, col));
				value += adjustValue;
			}
		}

		for (int row = 0; row < numRows; row++)
			for (int col = 0; col < numCols; col++)
				assertEquals(row + "x" + col, denseEmptyArray.get(row, col), shadow.get(row, col));
	}

	@Test
	public void testDenseInt2DArray() {
		DenseInt2DArray array = new DenseInt2DArray(numRows, numCols);
		testSparseIterator(array);
	}

	@Test
	public void testSparseInt2DArray() {
		SparseInt2DArray array = new SparseInt2DArray(numRows, numCols);
		testSparseIterator(array);
	}

	private static void testSparseIterator(Int2DArray array) {
		// assume that any non-zero value is stored in the array
		for (int i = 0; i < numRows - 1; i++) {
			for (int j = 0; j < numCols - 1; j++)
				array.set(i, j, 1);
			array.set(i, numCols - 1, 0);
		}
		for (int i = 0; i < numCols; i++)
			array.set(numRows - 1, i, 0);

		int sum = 0;
		Sparse2DIterator<Integer> it = array.sparseIterator();
		while (it.hasNext()) {
			it.advance();
			assertTrue(it.getColumn() >= 0);
			assertTrue(it.getColumn() < numCols);
			assertTrue(it.getRow() >= 0);
			assertTrue(it.getRow() < numRows);
			sum += it.getValue();
		}
		assertEquals((numRows - 1) * (numCols - 1), sum);
	}

}
