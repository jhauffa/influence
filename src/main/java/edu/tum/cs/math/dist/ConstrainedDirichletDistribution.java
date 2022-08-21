package edu.tum.cs.math.dist;

import java.util.logging.Logger;

import org.apache.commons.math3.special.Beta;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;

import edu.tum.cs.math.FastMath;

public class ConstrainedDirichletDistribution {

	private static final Logger logger = Logger.getLogger(ConstrainedDirichletDistribution.class.getName());

	// sampling parameters
	private static final int numBurnIn = 30;
	private static final int numSampleLag = 2;
	private static final double eps = 1E-21;
	private static final double constrEps = 1E-14;
	private static final UniformRandomProvider defaultRand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S);

	// distribution parameters
	private final double[] alpha;
	private final double[] c;
	private final int n;

	// sampling state
	final double[] x;
	private final double[] xi = new double[7];
	private final double[] xm = new double[7];
	private final double[] l = new double[7];
	private final double[] u = new double[7];
	private final int[] type = new int[7];
	private final double[] cc = new double[7];
	private final double[] p = new double[7];
	private final double[] pCumul = new double[7];

	private final UniformRandomProvider rand;
	private int numGibbsIter = 0;

	public ConstrainedDirichletDistribution(double[] alpha, double[] c, double r, UniformRandomProvider rand,
			boolean strict) {
		this.alpha = alpha;
		this.c = c;
		this.rand = (rand == null) ? defaultRand : rand;
		n = alpha.length;

		if (strict) {
			// sanity check for constraints
			if (n < 3)
				throw new IllegalArgumentException("3 or more dimensions required (got " + n + ")");
			double sum = 0.0;
			double min = 1.0;
			for (int i = 0; i < n; i++) {
				if (c[i] < -constrEps)
					throw new IllegalArgumentException("component " + i + " of c is negative: " + c[i]);
				sum += c[i];
				if (c[i] < min)
					min = c[i];
			}
			if (Math.abs(1.0 - sum) > constrEps)
				throw new IllegalArgumentException("c does not sum to 1: " + sum);
			double maxR = 2.0 - (2.0 * min);
			if (r > (maxR + constrEps)) {
				throw new IllegalArgumentException("L1 n-sphere does not intersect the unit simplex: r = " + r +
						" > " + maxR);
			}
		}

		x = initialValue(r);
	}

	public ConstrainedDirichletDistribution(double[] alpha, double[] c, double r, UniformRandomProvider rand) {
		this(alpha, c, r, rand, false);
	}

	public ConstrainedDirichletDistribution(double[] alpha, double[] c, double r) {
		this(alpha, c, r, null, false);
	}

	private double[] initialValue(double r) {
		int idxN = 0;
		for (int i = 1; i < n; i++) {
			if (c[i] < c[idxN])
				idxN = i;
		}

		// generate initial point within constraint region
		double[] x = new double[n];
		double sumBelow = 0.0;
		double sum = 0.0;
		for (int i = 0; i < n; i++) {
			if (i == idxN)
				continue;
			x[i] = Math.max(0.0, c[i] - ((r / 2.0) - sumBelow));
			sumBelow += c[i] - x[i];
			sum += x[i];
		}
		x[idxN] = 1.0 - sum;
		return x;
	}

	public double[] sample() {
		int numSteps = (numGibbsIter == 0) ? numBurnIn : numSampleLag;
		for (int i = 0; i < numSteps; i++) {
			gibbsStep();
			numGibbsIter++;
		}
		return x.clone();
	}

	private void gibbsStep() {
		final int idxM = n - 2;
		final int idxN = n - 1;

		for (int i = 0; i < idxM; i++) {
			double sumAbove = 0.0, sumBelow = 0.0;
			if (x[i] > c[i])
				sumAbove += x[i] - c[i];
			else
				sumBelow += c[i] - x[i];
			if (x[idxM] > c[idxM])
				sumAbove += x[idxM] - c[idxM];
			else
				sumBelow += c[idxM] - x[idxM];
			if (x[idxN] > c[idxN])
				sumAbove += x[idxN] - c[idxN];
			else
				sumBelow += c[idxN] - x[idxN];
			double sumRemaining = 1.0 - (x[i] + x[idxM] + x[idxN]);

			int s = 0;
			if (sumAbove < constrEps) {
				type[s] = 0;
				xm[s] = x[idxM];
				l[s] = Math.max(0.0, c[i] - sumBelow + c[idxM] - xm[s]);
				u[s] = Math.min(c[i] - sumBelow + c[idxM] - xm[s] + c[idxN], c[i]);
				cc[s] = 1.0 - xm[s] - sumRemaining;
				p[s] = 1.0;
				s++;
			} else {
				if (sumBelow < constrEps) {
					type[s] = 0;
					xm[s] = x[idxM];
					l[s] = c[i];
					u[s] = c[i] + sumAbove - xm[s] + c[idxM];
					cc[s] = 1.0 - xm[s] - sumRemaining;
					p[s] = segmentWeight(l[s], u[s], alpha[i], alpha[idxN], alpha[idxM], cc[s], xm[s]);
					s++;
				}
				if ((sumBelow - constrEps) <= c[i]) {
					type[s] = 2;
					xi[s] = c[i] - sumBelow;
					l[s] = c[idxM];
					u[s] = c[idxM] + sumAbove;
					cc[s] = 1.0 - xi[s] - sumRemaining;
					p[s] = segmentWeight(l[s], u[s], alpha[idxM], alpha[idxN], alpha[i], cc[s], xi[s]);
					s++;
				}
				if ((sumBelow - constrEps) <= (c[idxN] + c[i])) {
					type[s] = 0;
					xm[s] = c[idxM] + sumAbove;
					l[s] = Math.max(0.0, c[i] - sumBelow);
					u[s] = Math.min(c[i] - sumBelow + c[idxN], c[i]);
					cc[s] = 1.0 - xm[s] - sumRemaining;
					p[s] = segmentWeight(l[s], u[s], alpha[i], alpha[idxN], alpha[idxM], cc[s], xm[s]);
					s++;
				}
				if ((sumBelow - constrEps) <= (c[idxM] + c[i])) {
					type[s] = 1;
					xm[s] = c[i] - sumBelow + c[idxM];
					l[s] = Math.max(0.0, c[i] - sumBelow);
					u[s] = Math.min(c[i] - sumBelow + c[idxM], c[i]);
					cc[s] = c[i] - sumBelow + c[idxM];
					p[s] = segmentWeight(l[s], u[s], alpha[i], alpha[idxM], alpha[idxN], cc[s],
							1.0 - xm[s] - sumRemaining);
					s++;
				}
				if ((sumBelow - constrEps) <= (c[idxM] + c[idxN])) {
					type[s] = 2;
					xi[s] = c[i] + sumAbove;
					l[s] = Math.max(0.0, c[idxM] - sumBelow);
					u[s] = Math.min(c[idxM] - sumBelow + c[idxN], c[idxM]);
					cc[s] = 1.0 - xi[s] - sumRemaining;
					p[s] = segmentWeight(l[s], u[s], alpha[idxM], alpha[idxN], alpha[i], cc[s], xi[s]);
					s++;
				}
				if ((sumBelow - constrEps) <= c[idxN]) {
					type[s] = 1;
					xm[s] = c[i] + sumAbove + c[idxM];
					l[s] = c[i];
					u[s] = c[i] + sumAbove;
					cc[s] = c[i] + sumAbove + c[idxM];
					p[s] = segmentWeight(l[s], u[s], alpha[i], alpha[idxM], alpha[idxN], cc[s],
							1.0 - xm[s] - sumRemaining);
					s++;
				}
				if ((sumBelow - constrEps) <= c[idxM]) {
					type[s] = 0;
					xm[s] = c[idxM] - sumBelow;
					l[s] = c[i];
					u[s] = c[i] + sumAbove;
					cc[s] = 1.0 - xm[s] - sumRemaining;
					p[s] = segmentWeight(l[s], u[s], alpha[i], alpha[idxN], alpha[idxM], cc[s], xm[s]);
					s++;
				}
			}

			// Due to numerical inaccuracy, there may be no matching segments. In this case, keep the current parameters
			// and try again with a different component i.
			if (s > 0) {
				int sIdx = 0;
				if (s > 1) {
					double pSum = 0.0;
					for (int j = 0; j < s; j++) {
						if (Double.isNaN(p[j]))
							throw new RuntimeException("segment " + j + ": probability is NaN");
						pSum += p[j];
					}

					for (int j = 0; j < s; j++) {
						double pj = 1.0;
						if (pSum > eps)
							pj = p[j];
						if (j == 0)
							pCumul[j] = pj;
						else
							pCumul[j] = pCumul[j - 1] + pj;
					}
					sIdx = DiscreteDistributionSampler.sample1D(rand.nextDouble(), pCumul, s);
				}

				if (type[sIdx] == 0) {
					x[i] = sampleX(l[sIdx], u[sIdx], alpha[i], alpha[idxN], cc[sIdx]);
					x[idxM] = xm[sIdx];
				} else if (type[sIdx] == 1) {
					x[i] = sampleX(l[sIdx], u[sIdx], alpha[i], alpha[idxM], cc[sIdx]);
					x[idxM] = xm[sIdx] - x[i];
				} else {
					x[i] = xi[sIdx];
					x[idxM] = sampleX(l[sIdx], u[sIdx], alpha[idxM], alpha[idxN], cc[sIdx]);
				}
				x[idxN] = 1.0 - (sumRemaining + x[i] + x[idxM]);
			} else {
				logger.warning("no matching segment, sumAbove = " + sumAbove + ", sumBelow = " + sumBelow);
			}
		}
	}

	private double sampleX(double l, double u, double alphaI, double alphaN, double c) {
		if (u <= eps)
			return 0.0;
		if (l >= u)
			return u;
		if (c < eps)
			return l + (rand.nextDouble() * (u - l));

		// sample y from TBeta(alpha_i, alpha_n; l', u'), where l' = l/c, u' = u/c, x = c*y
		double ly = Math.max(0.0, l / c);
		double uy = Math.min(u / c, 1.0);
		return c * TruncatedBetaDistribution.sample(alphaI, alphaN, ly, uy, rand);
	}

	private static double segmentWeight(double l, double u, double alphaI, double alphaN, double alphaM, double c,
			double d) {
		if ((l >= u) || (c < eps))
			return 0.0;

		double ly = Math.max(0.0, l / c);
		double uy = Math.min(u / c, 1.0);
		d = Math.pow(Math.max(constrEps, d), alphaM - 1.0);

		if (alphaI == 1.0) {
			if (alphaN == 1.0)
				return d * c * (uy - ly);
			return d * c * (FastMath.pow(1.0 - ly, alphaN) - FastMath.pow(1.0 - uy, alphaN)) / (alphaN * alphaN);
		} else if (alphaN == 1.0) {
			return d * c * (FastMath.pow(uy, alphaI) - FastMath.pow(ly, alphaI)) / (alphaI * alphaI);
		} else {
			double w = d * c * FastMath.exp(Beta.logBeta(alphaI, alphaN));
			if (uy == 1.0) {
				if (ly == 0.0)
					return w;
				return w * (1.0 - Beta.regularizedBeta(ly, alphaI, alphaN));
			} else if (ly == 0.0)
				return w * Beta.regularizedBeta(uy, alphaI, alphaN);
			return w * (Beta.regularizedBeta(uy, alphaI, alphaN) - Beta.regularizedBeta(ly, alphaI, alphaN));
		}
	}

}
