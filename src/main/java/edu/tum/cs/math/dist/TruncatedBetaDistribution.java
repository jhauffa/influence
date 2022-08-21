package edu.tum.cs.math.dist;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;

import edu.tum.cs.math.FastMath;

public class TruncatedBetaDistribution {

	/**
	 * Trade-off between speed and closeness to the true distribution; 24 iterations yield reasonably unbiased samples,
	 * as evidenced by the unit tests TruncatedBetaDistributionTest and ConstrainedDirichletDistributionTest.
	 */
	private static final int minGibbsIter = 24;

	private static final double logZero = 1E-16;
	private static final UniformRandomProvider defaultRand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S);

	private final double alpha, beta;
	private final double l, u;
	private final double startX;
	private final int numIter;
	private final UniformRandomProvider rand;

	public TruncatedBetaDistribution(double alpha, double beta, double l, double u, UniformRandomProvider rand) {
		this.alpha = alpha;
		this.beta = beta;
		this.l = l;
		this.u = u;
		this.rand = (rand == null) ? defaultRand : rand;

		numIter = computeNumIter(alpha, beta, l, u);
		startX = computeStartX(alpha, beta, l, u, numIter);
	}

	public TruncatedBetaDistribution(double alpha, double beta, double l, double u) {
		this(alpha, beta, l, u, null);
	}

	private static int computeNumIter(double alpha, double beta, double l, double u) {
		if ((alpha > 1.0) && (beta < 1.0)) {
			// setting numIter to ~ 11 * minGibbsIter would also work
			double d = Math.max(1E-24, Math.min(l, 1.0 - u));
			return Math.max(minGibbsIter, (int) Math.floor(-6.0 * (7.0 + FastMath.log(d))));
		}
		return minGibbsIter;
	}

	private static double computeStartX(double alpha, double beta, double l, double u, int numIter) {
		double x = (l + u) / 2.0;
		if ((alpha > 1.0) && (beta < 1.0))
			return Math.max(x, u - FastMath.pow(numIter, -1.2));
		return x;
	}

	public double median() {
		BetaDistribution dist = new BetaDistribution(alpha, beta);
		return dist.inverseCumulativeProbability((dist.cumulativeProbability(l) + dist.cumulativeProbability(u)) / 2.0);
	}

	public double sample() {
		return sample(startX, numIter, alpha, beta, l, u, rand);
	}

	public static double sample(double alpha, double beta, double l, double u, UniformRandomProvider rand) {
		int numIter = computeNumIter(alpha, beta, l, u);
		double startX = computeStartX(alpha, beta, l, u, numIter);
		return sample(startX, numIter, alpha, beta, l, u, rand);
	}

	public static double sample(double startX, int numIter, double alpha, double beta, double l, double u,
			UniformRandomProvider rand) {
		// special cases that can be handled by simple CDF inversion
		if (alpha == 1.0) {
			if (beta == 1.0)
				return l + (rand.nextDouble() * (u - l));
			return 1.0 - sampleX(1.0 - u, 1.0 - l, beta, 1.0 / beta, rand);
		} else if (beta == 1.0) {
			return sampleX(l, u, alpha, 1.0 / alpha, rand);
		}

		double x = startX;

		if ((alpha < 1.0) && (beta < 1.0)) {
			// sample from logit-transformed Beta distribution
			x = FastMath.log(x) - FastMath.log(1.0 - x);
			double lt = Math.max(logZero, l);
			lt = FastMath.log(lt) - FastMath.log(1.0 - lt);
			double ut = Math.min(u, 1.0 - logZero);
			ut = FastMath.log(ut) - FastMath.log(1.0 - ut);

			final double invAB = 1.0 / (alpha + beta);
			for (int i = 0; i < numIter; i++) {
				double p = 1.0 - rand.nextDouble();
				double b = FastMath.log1p(FastMath.exp(x) - Math.pow(p, invAB)) - (invAB * FastMath.log(p));
				// Fewer operations, but probably less accurate:
				// double b = FastMath.log(((1.0 + FastMath.exp(x)) / Math.pow(p, invAB)) - 1.0);

				double utb = Math.min(b, ut);
				p = rand.nextDouble();
				x = lt + (FastMath.log((1.0 - p) + p * FastMath.exp(alpha * (utb - lt))) / alpha);
			}
			return 1.0 / (1.0 + FastMath.exp(-x));
		}

		final double invAlpha = 1.0 / alpha;
		final double invBeta1 = 1.0 / (beta - 1.0);
		for (int i = 0; i < numIter; i++) {
			double p = rand.nextDouble();
			double yPow = FastMath.pow(p, invBeta1) * (1.0 - x);

			double ll = l, uu = u;
			if (beta < 1.0)
				ll = Math.max(ll, 1.0 - yPow);
			else
				uu = Math.min(1.0 - yPow, uu);
			x = sampleX(ll, uu, alpha, invAlpha, rand);
		}
		return x;
	}

	private static double sampleX(double l, double u, double alpha, double invAlpha, UniformRandomProvider rand) {
		double p = rand.nextDouble();
		return u * FastMath.pow(p + (1.0 - p) * FastMath.pow(l / u, alpha), invAlpha);
	}

}
