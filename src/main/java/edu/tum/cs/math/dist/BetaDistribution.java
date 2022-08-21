package edu.tum.cs.math.dist;

import java.util.Arrays;
import java.util.logging.Logger;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.ContinuedFraction;
import org.apache.commons.math3.util.Precision;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;

import edu.tum.cs.math.FastMath;
import edu.tum.cs.util.RngAdaptor;

public class BetaDistribution extends org.apache.commons.math3.distribution.BetaDistribution {

	private static final long serialVersionUID = 6202366282043803508L;
	private static final Logger logger = Logger.getLogger(BetaDistribution.class.getName());
	private static final double eps = 1E-14;

	public BetaDistribution(double alpha, double beta) {
		this(alpha, beta, RandomSource.create(RandomSource.XOR_SHIFT_1024_S));
	}

	public BetaDistribution(double alpha, double beta, UniformRandomProvider rand) {
		super(new RngAdaptor(rand), alpha, beta);
	}

	public static class BetaSummaryStatistics {
		private final boolean computeLogSum;
		private final double minValueMl;
		private final Mean mean = new Mean();
		private final Variance variance = new Variance();
		private double sumLog, sumLogInv;
		private int n;

		public BetaSummaryStatistics(boolean computeLogSum, int expectedN) {
			this.computeLogSum = computeLogSum;
			// could be simplified to (1/x)^ln(10)
			minValueMl = computeLogSum ? Math.min(0.05, FastMath.pow(10.0, -FastMath.log(expectedN))) : 0.0;
		}

		public void addValue(double v) {
			mean.increment(v);
			variance.increment(v);
			if (computeLogSum) {
				v = Math.max(minValueMl, v);
				v = Math.min(v, 1.0 - minValueMl);
				sumLog += FastMath.log(v);
				sumLogInv += FastMath.log(1.0 - v);
			}
			n++;
		}

		public double getMean() {
			return mean.getResult();
		}

		public double getVariance() {
			return variance.getResult();
		}

		public double getAvgLog() {
			return sumLog / n;
		}

		public double getAvgLogInv() {
			return sumLogInv / n;
		}

		public int getNumSamples() {
			return n;
		}
	}

	private static final int maxIterFit = 1000;
	private static final double tolFit = Math.sqrt(eps);

	public static BetaDistribution fit(double[] data) {
		return fit(data, true);
	}

	public static BetaDistribution fit(double[] data, boolean useMl) {
		BetaSummaryStatistics stats = new BetaSummaryStatistics(useMl, data.length);
		for (double v : data)
			stats.addValue(v);
		return fit(stats, useMl);
	}

	public static BetaDistribution fit(BetaSummaryStatistics stats, boolean useMl) {
		// Return uniform distribution in case of no information.
		if (stats.n == 0)
			return new BetaDistribution(1.0, 1.0);

		// Use method of moments to get an initial estimate of parameters alpha and beta.
		double mean = stats.getMean();
		double f = ((mean * (1.0 - mean)) / stats.getVariance()) - 1.0;
		double alpha = mean * f;
		double beta = (1.0 - mean) * f;
		if ((alpha < 0.0) || (beta < 0.0)) {
			logger.warning("invalid MoM estimate " + alpha + ", " + beta);
			alpha = beta = 1.0;
		}

		if (useMl) {
			double origAlpha = alpha;
			double origBeta = beta;

			// Use maximum likelihood estimation to refine the parameter values. Less robust than MoM if data points lie
			// exactly on the boundaries of the support (i.e. 0 and 1).
			double avgLog = stats.getAvgLog();
			double avgLogInv = stats.getAvgLogInv();

			// Newton-Raphson iteration
			double delta = Double.POSITIVE_INFINITY;
			for (int i = 0; i < maxIterFit; i++) {
				double tgAlpha = Gamma.trigamma(alpha);
				double tgBeta = Gamma.trigamma(beta);
				double tgAlphaBeta = Gamma.trigamma(alpha + beta);
				double c = 1.0 / (((tgAlpha - tgAlphaBeta) * (tgBeta - tgAlphaBeta)) - (tgAlphaBeta * tgAlphaBeta));

				double dgAlpha = Gamma.digamma(alpha);
				double dgBeta = Gamma.digamma(beta);
				double dgAlphaBeta = Gamma.digamma(alpha + beta);
				double alphaUpd = (c * (tgBeta - tgAlphaBeta) * (dgAlpha - dgAlphaBeta - avgLog)) +
						(c * tgAlphaBeta * (dgBeta - dgAlphaBeta - avgLogInv));
				double betaUpd = (c * tgAlphaBeta * (dgAlpha - dgAlphaBeta - avgLog)) +
						(c * (tgAlpha - tgAlphaBeta) * (dgBeta - dgAlphaBeta - avgLogInv));

				alpha -= alphaUpd;
				beta -= betaUpd;
				if ((alpha < 0.0) || (beta < 0.0)) {
					logger.warning("ML estimation failed, using MoM estimate");
					return new BetaDistribution(origAlpha, origBeta);
				}
				delta = Math.sqrt((alphaUpd * alphaUpd) + (betaUpd * betaUpd));
				if (delta < tolFit)
					break;
			}
			if (delta >= tolFit) {
				logger.warning("L2 norm > epsilon (" + delta + " > " + tolFit + ") after " + maxIterFit +
					" Newton iterations");
			}
		}
		return new BetaDistribution(alpha, beta);
	}

	/* implementation adapted from Apache Commons Math */
	private static double regularizedBeta(double x, final double a, final double b, double abTerm, double baTerm,
			double epsilon) {
		if (Double.isNaN(x) || Double.isNaN(a) || Double.isNaN(b) || x < 0 || x > 1 || a <= 0 || b <= 0)
			return Double.NaN;
		if (x > (a + 1) / (2 + b + a) && 1 - x <= (b + 1) / (2 + b + a))
			return 1 - regularizedBeta(1 - x, b, a, baTerm, abTerm, epsilon);

		ContinuedFraction fraction = new ContinuedFraction() {
			@Override
			protected double getB(int n, double x) {
				double ret;
				double m;
				if (n % 2 == 0) { // even
					m = n / 2.0;
					ret = (m * (b - m) * x) / ((a + (2 * m) - 1) * (a + (2 * m)));
				} else {
					m = (n - 1.0) / 2.0;
					ret = -((a + m) * (a + b + m) * x) / ((a + (2 * m)) * (a + (2 * m) + 1.0));
				}
				return ret;
			}

			@Override
			protected double getA(int n, double x) {
				return 1.0;
			}
		};
		if (Double.isNaN(abTerm))
			abTerm = -FastMath.log(a) - Beta.logBeta(a, b);
		return FastMath.exp((a * FastMath.log(x)) + (b * FastMath.log1p(-x)) + abTerm) *
				1.0 / fraction.evaluate(x, epsilon, Integer.MAX_VALUE);
	}

	/* implementation adapted from Apache Commons Math */
	private double kolmogorovSmirnovStatistic(double[] data) {
		Arrays.sort(data);
		double a = getAlpha();
		double b = getBeta();
		double abTerm = -FastMath.log(a) - Beta.logBeta(a, b);
		double baTerm = -FastMath.log(b) - Beta.logBeta(b, a);
		double d = 0.0;
		for (int i = 1; i <= data.length; i++) {
			// open-coded evaluation of Beta distribution CDF for value data[i - 1]
			double yi;
			if (data[i - 1] <= 0.0)
				yi = 0.0;
			else if (data[i - 1] >= 1.0)
				yi = 1.0;
			else
				yi = regularizedBeta(data[i - 1], a, b, abTerm, baTerm, eps);

			double curD = Math.max(yi - ((double) (i - 1) / data.length), ((double) i / data.length) - yi);
			if (curD > d)
				d = curD;
		}
		return d;
	}

	/* implementation adapted from Apache Commons Math */
	private static final class ChengSampler {
		static double sample(RandomGenerator random, final double alpha, final double beta) {
			final double a = Math.min(alpha, beta);
			final double b = Math.max(alpha, beta);

			if (a > 1.)
				return algorithmBB(random, alpha, a, b);
			return algorithmBC(random, alpha, b, a);
		}

		private static double algorithmBB(RandomGenerator random, final double a0, final double a, final double b) {
			final double alpha = a + b;
			final double beta = Math.sqrt((alpha - 2.) / (2. * a * b - alpha));
			final double gamma = a + 1. / beta;

			double r;
			double w;
			double t;
			do {
				final double u1 = random.nextDouble();
				final double u2 = random.nextDouble();
				final double v = beta * (FastMath.log(u1) - FastMath.log1p(-u1));
				w = a * FastMath.exp(v);
				final double z = u1 * u1 * u2;
				r = gamma * v - 1.3862944;
				final double s = a + r - w;
				if (s + 2.609438 >= 5 * z)
					break;

				t = FastMath.log(z);
				if (s >= t)
					break;
			} while (r + alpha * (FastMath.log(alpha) - FastMath.log(b + w)) < t);

			w = Math.min(w, Double.MAX_VALUE);
			return Precision.equals(a, a0) ? w / (b + w) : b / (b + w);
		}

		private static double algorithmBC(RandomGenerator random, final double a0, final double a, final double b) {
			final double alpha = a + b;
			final double beta = 1. / b;
			final double delta = 1. + a - b;
			final double k1 = delta * (0.0138889 + 0.0416667 * b) / (a * beta - 0.777778);
			final double k2 = 0.25 + (0.5 + 0.25 / delta) * b;

			double w;
			for (;;) {
				final double u1 = random.nextDouble();
				final double u2 = random.nextDouble();
				final double y = u1 * u2;
				final double z = u1 * y;
				if (u1 < 0.5) {
					if (0.25 * u2 + z - y >= k1)
						continue;
				} else {
					if (z <= 0.25) {
						final double v = beta * (FastMath.log(u1) - FastMath.log1p(-u1));
						w = a * FastMath.exp(v);
						break;
					}

					if (z >= k2)
						continue;
				}

				final double v = beta * (FastMath.log(u1) - FastMath.log1p(-u1));
				w = a * FastMath.exp(v);
				if (alpha * (FastMath.log(alpha) - FastMath.log(b + w) + v) - 1.3862944 >= FastMath.log(z))
					break;
			}

			w = Math.min(w, Double.MAX_VALUE);
			return Precision.equals(a, a0) ? w / (b + w) : b / (b + w);
		}
	}

	@Override
	public double sample() {
		return ChengSampler.sample(random, getAlpha(), getBeta());
	}

	/**
	 * Testing the null hypothesis that the given sample was drawn from the Beta distribution represented by this
	 * instance. If the returned p-value is significant, the null hypothesis can be rejected, which is evidence against
	 * the sample having been drawn from this distribution.
	 */
	public double testGoodnessOfFit(double[] values, int numMcSamples) {
		if (values.length == 0)
			return 1.0;

		double d = kolmogorovSmirnovStatistic(values);
		int n = 0;
		for (int i = 0; i < numMcSamples; i++) {
			double[] valuesSample = new double[values.length];
			for (int j = 0; j < valuesSample.length; j++)
				valuesSample[j] = sample();
			BetaDistribution distSample = BetaDistribution.fit(valuesSample, false);
			double dSample = distSample.kolmogorovSmirnovStatistic(valuesSample);
			if (dSample > d)
				n++;
		}
		return (double) (1 + n) / (numMcSamples + 1);
	}

}
