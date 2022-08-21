package edu.tum.cs.math.dist;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.junit.Test;

public class TruncatedBetaDistributionTest {

	private static final double[] limits = { 0.0, 1E-10, 0.1, 0.5 - 1E-10 };
	private static final int minExp = -5;
	private static final int maxExp = 1;

	@Test
	public void testMedian() {
		for (double l : limits) {
			for (int i = minExp; i <= maxExp; i++) {
				double p = Math.pow(10.0, i);
				TruncatedBetaDistribution dist = new TruncatedBetaDistribution(p, p, l, 1.0 - l);
				assertEquals(0.5, dist.median(), 1E-7);
			}
		}
	}

	private static final long seed = 73736463L;
	private static final int numSamples = 10000;
	private static final double maxRelDiff = 0.16;


	private static void sampleDist(UniformRandomProvider rand, double alpha, double beta, double l, double u) {
		TruncatedBetaDistribution dist = new TruncatedBetaDistribution(alpha, beta, l, u, rand);
		double median = dist.median();
		// Median is computed using the inverse CDF, which suffers from numerical accuracy issues close to the
		// boundaries of the support. Skip the test if the median cannot be computed with sufficient accuracy.
		if (((l == 0.0) || (u == 1.0)) && ((median < 0.001) || (median > 0.999)))
			return;

		String id = alpha + ";" + beta + ";" + l + ";" + u + ";" + median;
		int numLeft = 0;
		for (int i = 0; i < numSamples; i++) {
			double s = dist.sample();
			assertFalse(id, Double.isNaN(s));
			assertTrue(id + ";" + s, s > (l - 1E-10));
			assertTrue(id + ";" + s, s < (u + 1E-10));
			if (s < median)
				numLeft++;
		}
		double relDiff = Math.abs(1.0 - ((2.0 * numLeft) / numSamples));
		assertTrue(id + ";" + numLeft + ";" + relDiff, relDiff < maxRelDiff);
	}

	@Test
	public void testSampling() {
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
		for (int i = minExp; i <= maxExp; i++) {
			for (int j = minExp; j <= maxExp; j++) {
				double alpha = Math.pow(10.0, i);
				double beta = Math.pow(10.0, j);
				for (double l : limits) {
					sampleDist(rand, alpha, beta, l, 1.0 - l);
					if (l > 0.0) {
						sampleDist(rand, alpha, beta, l, 1.0);
						sampleDist(rand, alpha, beta, 0.0, 1.0 - l);
					}
				}
			}
		}
	}

}
