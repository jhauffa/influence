package edu.tum.cs.math.dist;

import static org.junit.Assert.assertTrue;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.junit.Test;

public class BetaDistributionTest {

	private static final long seed = 1330L;

	private static final int numFitTests = 10;
	private static final double minParam = 0.001;
	private static final double maxParam = 5.0;
	private static final int numSamples = 10000;
	private static final double tolFitTest = 0.2;
	private static final double tolAvg = 0.1;
	private static final double fitTestAlpha = 0.01;
	private static final double maxRelSignificant = 0.3;

	@Test
	public void testFit() {
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
		double avgNorm = 0.0;
		int numSignificant = 0;
		for (int i = 0; i < numFitTests; i++) {
			// generate true parameters alpha, beta
			double trueAlpha = minParam + (rand.nextDouble() * (maxParam - minParam));
			double trueBeta = minParam + (rand.nextDouble() * (maxParam - minParam));

			// generate samples from Beta distribution
			BetaDistribution trueDist = new BetaDistribution(trueAlpha, trueBeta, rand);
			double[] data = new double[numSamples];
			for (int j = 0; j < numSamples; j++)
				data[j] = trueDist.sample();

			// estimate parameters from samples and compare parameters
			BetaDistribution estDist = BetaDistribution.fit(data);
			double norm = 0.0;
			double d = estDist.getAlpha() - trueAlpha;
			norm += d * d;
			d = estDist.getBeta() - trueBeta;
			norm += d * d;
			norm = Math.sqrt(norm);
			assertTrue(Double.toString(norm), norm < tolFitTest);
			avgNorm += norm;

			// test goodness-of-fit
			if (estDist.testGoodnessOfFit(data, 250) < fitTestAlpha)
				numSignificant++;
		}
		avgNorm /= numFitTests;
		assertTrue(Double.toString(avgNorm), avgNorm < tolAvg);
		assertTrue(Integer.toString(numSignificant), ((double) (numSignificant / numFitTests)) <= maxRelSignificant);
	}

}
