package edu.tum.cs.math;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.junit.Test;

public class MultipleHypothesisTestTest {

	private static class BasicTest extends MultipleHypothesisTest<Double, Double> {
		public BasicTest(double alpha, int expMaxIter) {
			super(alpha, expMaxIter, false);
		}

		@Override
		protected int generateSamples(int n, Double hypo, Double[] perIterationData, Double perHypothesisData) {
			int numExceedance = 0;
			for (int i = 0; i < n; i++) {
				if (Math.random() <= hypo)
					numExceedance++;
			}
			return numExceedance;
		}
	}

	private static final long seed = 1472347823L;
	private static final int numHypo = 10;
	private static final double baseSignificanceLevel = 0.01, maxAbsDiffSmall = 0.1, maxAbsDiffLarge = 0.5;
	private static final int maxIter = 1000;

	@Test
	public void testHypothesisTest() {
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
		Double[] hypo = new Double[numHypo];
		for (int i = 0; i < numHypo; i++)
			hypo[i] = rand.nextDouble();

		BasicTest test = new BasicTest(baseSignificanceLevel, maxIter);
		AtomicLong numTestsCompleted = new AtomicLong(0);
		double[] p = test.performTests(Arrays.asList(hypo), numTestsCompleted);
		assertEquals(numHypo, numTestsCompleted.get());
		for (int i = 0; i < numHypo; i++) {
			double absDiff = Math.abs(hypo[i] - p[i]);
			if (hypo[i] < 0.1)
				assertTrue(Double.toString(absDiff), absDiff < maxAbsDiffSmall);
			else
				assertTrue(Double.toString(absDiff), absDiff < maxAbsDiffLarge);
		}
	}


	// example data taken from Benjamini and Hochberg, 1995
	private static final double[] uncorrP = new double[] {
		0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344, 0.0459, 0.324, 0.4262, 0.5719, 0.6528, 0.759,
		1.0
	};
	private static final double alphaBH = 0.05;
	private static final int maxRejectIdx = 3;

	@Test
	public void testBenjaminiHochberg() {
		Integer[] rejectIdx = MultipleHypothesisTest.testSignificanceBH(uncorrP, alphaBH);
		assertEquals(maxRejectIdx + 1, rejectIdx.length);
		for (int i = 0; i <= maxRejectIdx; i++)
			assertEquals(i, (int) rejectIdx[i]);
	}

}
