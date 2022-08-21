package edu.tum.cs.influence.micro;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.junit.Test;

public class TimeSeriesDataTest {

	private static class Permutation extends TimeSeriesData {
		private static final long serialVersionUID = 1L;

		public Permutation(TimeSeriesData ts) {
			super(ts);
		}

		@Override
		public int hashCode() {
			return Arrays.hashCode(values);
		}

		@Override
		public boolean equals(Object o) {
			if (!(o instanceof Permutation))
				return false;
			return Arrays.equals(values, ((Permutation) o).values);
		}

		public boolean isValid() {
			Set<double[]> unique = new HashSet<double[]>(values.length);
			for (double[] v : values)
				unique.add(v);
			return (unique.size() == values.length);
		}
	}

	private static final int permLength = 7;

	@Test
	public void testPermutation() {
		TimeSeriesData data = new TimeSeriesData(permLength);
		data.setMissingValues(new double[0]);	// clones input -> unique default hashCode/equals

		// test enumerating all permutations
		int numPerm = (int) CombinatoricsUtils.factorial(permLength);
		Set<Permutation> unique = new HashSet<Permutation>(numPerm);
		int[] indices = data.getIndices();
		assertEquals(permLength, indices.length);
		TimeSeriesData dst = new TimeSeriesData(data);
		for (int i = 0; i < numPerm; i++) {
			data.permute(indices, i, dst);
			Permutation p = new Permutation(dst);
			assertTrue(p.isValid());
			unique.add(p);
		}
		assertEquals(numPerm, unique.size());
	}

}
