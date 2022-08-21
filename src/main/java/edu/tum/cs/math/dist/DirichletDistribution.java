package edu.tum.cs.math.dist;

import java.util.Arrays;
import java.util.logging.Logger;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.BoxMullerGaussianSampler;
import org.apache.commons.rng.simple.RandomSource;

import edu.tum.cs.util.arrays.DoubleRowVector;

/**
 * Computation of log-likelihood and parameter fitting of the Dirichlet distribution according to Minka, 2000,
 * "Estimating a Dirichlet distribution".
 */
public class DirichletDistribution {

	private static final Logger logger = Logger.getLogger(DirichletDistribution.class.getName());
	private static final UniformRandomProvider defaultRand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S);

	private final DoubleRowVector alpha;
	private final UniformRandomProvider rand;
	private final BoxMullerGaussianSampler gaussian;

	public DirichletDistribution(int n, double alphaUniform) {
		this(n, alphaUniform, null);
	}

	public DirichletDistribution(int n, double alphaUniform, UniformRandomProvider rand) {
		double[] a = new double[n];
		Arrays.fill(a, alphaUniform);
		alpha = new DoubleRowVector(a);
		this.rand = (rand == null) ? defaultRand : rand;
		gaussian = new BoxMullerGaussianSampler(rand, 0, 1);
	}

	public DirichletDistribution(double[] alpha) {
		this(alpha, null);
	}

	public DirichletDistribution(double[] alpha, UniformRandomProvider rand) {
		this(new DoubleRowVector(alpha), rand);
	}

	public DirichletDistribution(DoubleRowVector alpha) {
		this(alpha, null);
	}

	public DirichletDistribution(DoubleRowVector alpha, UniformRandomProvider rand) {
		this.alpha = alpha;
		this.rand = (rand == null) ? defaultRand : rand;
		gaussian = new BoxMullerGaussianSampler(rand, 0, 1);
	}

	public double[] getAlpha() {
		return alpha.toArray();
	}

	private double sampleGammaAhrensDieter(double alpha, double theta) {
		double oneOverAlpha = 1 / alpha;
		double bGSOptim = 1 + alpha / Math.E;

		while (true) {
			double u = rand.nextDouble();
			double p = bGSOptim * u;
			if (p <= 1) {
				double x = Math.pow(p, oneOverAlpha);
				double u2 = rand.nextDouble();
				if (u2 > Math.exp(-x))
					continue;
				return theta * x;
			}

			double x = -Math.log((bGSOptim - p) * oneOverAlpha);
			double u2 = rand.nextDouble();
			if (u2 <= Math.pow(x, alpha - 1))
				return theta * x;
		}
	}

	private double sampleGammaMarsagliaTsang(double alpha, double theta) {
		double dOptim = alpha - (1d / 3);
		double cOptim = (1d / 3) / Math.sqrt(dOptim);

		while (true) {
			double x = gaussian.sample();
			double oPcTx = 1 + cOptim * x;
			double v = oPcTx * oPcTx * oPcTx;
			if (v <= 0)
				continue;

			double x2 = x * x;
			double u = rand.nextDouble();
			if (u < 1 - 0.0331 * x2 * x2)
				return theta * dOptim * v;
			if (Math.log(u) < 0.5 * x2 + dOptim * (1 - v + Math.log(v)))
				return theta * dOptim * v;
		}
	}

	private double sampleGamma(double alpha, double theta) {
		// This code and the two functions it calls have been adapted from Apache Commons RNG, so we don't have to
		// create a separate object for each combination of parameters.
		if (alpha < 1)
			return sampleGammaAhrensDieter(alpha, theta);
		return sampleGammaMarsagliaTsang(alpha, theta);
	}

	public double[] sample(double[] v) {
		int n = alpha.size();
		double sum = 0.0;
		for (int i = 0; i < n; i++) {
			v[i] = sampleGamma(alpha.get(i), 1.0);
			sum += v[i];
		}
		for (int i = 0; i < n; i++)
			v[i] /= sum;
		return v;
	}

	public double[] sample() {
		double[] v = new double[alpha.size()];
		return sample(v);
	}

	public double logLikelihood(double[][] data) {
		int m = data.length;
		double s1 = 0.0, s2 = 0.0, s3 = 0.0;
		for (int i = 0; i < alpha.size(); i++) {
			double ai = alpha.get(i);
			s1 += ai;
			s2 += Gamma.logGamma(ai);
			double sd = 0.0;
			for (int j = 0; j < m; j++)
				sd += Math.log(data[j][i] > 0.0 ? data[j][i] : Double.MIN_NORMAL);
			s3 += (ai - 1.0) * (sd / m);
		}
		return m * (Gamma.logGamma(s1) - s2 + s3);
	}

	private static final int maxIterFit = 1000;
	private static final double tolFit = 1E-7;

	public static DirichletDistribution fit(double[][] data) {
		int m = data.length;
		int n = data[0].length;

		double[] e = new double[n];
		double e2 = 0.0;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++)
				e[j] += data[i][j] / n;
			e2 += (data[i][0] * data[i][0]) / n;
		}
		double s = (e[0] - e2) / (e2 - (e[0] * e[0]));
		double[] alpha = new double[n];
		for (int i = 0; i < n; i++)
			alpha[i] = s * e[i];

		double norm = 0.0;
		for (int i = 0; i < maxIterFit; i++) {
			double s1 = 0.0;
			for (int j = 0; j < n; j++)
				s1 += alpha[j];
			s1 = Gamma.digamma(s1);

			norm = 0.0;
			for (int j = 0; j < n; j++) {
				double s2 = 0.0;
				for (int k = 0; k < m; k++)
					s2 += Math.log(data[k][j] > 0.0 ? data[k][j] : Double.MIN_NORMAL);
				double alphaNew = invDigamma(s1 + (s2 / m));
				double d = alpha[j] - alphaNew;
				norm += d * d;
				alpha[j] = alphaNew;
			}
			norm = Math.sqrt(norm);
			if (norm < tolFit)
				break;
		}
		if (norm >= tolFit) {
			logger.warning("L2 norm > epsilon (" + norm + " > " + tolFit + ") after " + maxIterFit +
					" fixed-point iterations");
		}
		return new DirichletDistribution(alpha);
	}

	private static final int maxIterInvDigamma = 10;
	private static final double tolInvDigamma = 1E-12;

	private static double invDigamma(double y) {
		double x = (y < -2.22) ? -1.0 / (y + Gamma.GAMMA) : Math.exp(y) + 0.5;
		double norm = 0.0;
		for (int i = 0; i < maxIterInvDigamma; i++) {
			double xUpd = (Gamma.digamma(x) - y) / Gamma.trigamma(x);
			x = x - xUpd;
			norm = Math.abs(xUpd);
			if (norm < tolInvDigamma)
				break;
		}
		if (norm >= tolInvDigamma) {
			logger.warning("L1 norm > epsilon (" + norm + " > " + tolInvDigamma + ") after " + maxIterInvDigamma +
					" fixed-point iterations");
		}
		return x;
	}

}
