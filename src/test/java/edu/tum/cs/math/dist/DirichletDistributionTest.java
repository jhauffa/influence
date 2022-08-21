package edu.tum.cs.math.dist;

import static org.junit.Assert.assertTrue;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.junit.Test;

import edu.tum.cs.nlp.topic.model.TopicModelTest;

public class DirichletDistributionTest {

	private static final long seed = 1339L;

	@Test
	public void testSampling() {
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
		DirichletDistribution dist = new DirichletDistribution(10, 5.0, rand);
		double[] p = dist.sample();
		TopicModelTest.checkProbabilityDistr(p);
	}

	private static final int numFitTests = 10;
	private static final int numDim = 20;
	private static final double minS = 0.001;
	private static final double maxS = 5.0;
	private static final int numSamples = 1000;
	private static final double tolFitTest = 0.1;

	@Test
	public void testFit() {
		DirichletDistribution flatDist = new DirichletDistribution(numDim, 1.0);
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
		for (int i = 0; i < numFitTests; i++) {
			// generate true parameter vector alpha
			double s = minS + (rand.nextDouble() * (maxS - minS));
			double[] trueAlpha = new double[numDim];
			double sum = 0.0;
			for (int j = 0; j < numDim; j++) {
				trueAlpha[j] = rand.nextDouble();
				sum += trueAlpha[j];
			}
			for (int j = 0; j < numDim; j++)
				trueAlpha[j] = (trueAlpha[j] / sum) * s;

			// generate samples from Dirichlet distribution
			DirichletDistribution trueDist = new DirichletDistribution(trueAlpha, rand);
			double[][] data = new double[numSamples][];
			for (int j = 0; j < numSamples; j++)
				data[j] = trueDist.sample();

			// estimate parameters from samples and compare parameters
			DirichletDistribution estDist = DirichletDistribution.fit(data);
			double norm = 0.0;
			double[] estAlpha = estDist.getAlpha();
			for (int j = 0; j < estAlpha.length; j++) {
				double d = estAlpha[j] - trueAlpha[j];
				norm += d * d;
			}
			norm = Math.sqrt(norm);
			assertTrue(Double.toString(norm), norm < tolFitTest);

			// LL of fitted distribution should be higher than LL of flat distribution
			assertTrue(flatDist.logLikelihood(data) < estDist.logLikelihood(data));
		}
	}

}
