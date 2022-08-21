package edu.tum.cs.math.dist;

import java.util.Arrays;

import edu.tum.cs.math.FastMath;

public class DiscreteDistribution {

	public static final double LOG2 = Math.log(2.0);

	/**
	 * Jensen-Shannon divergence, base 2; 0 if p and q are equal, 1 if they are maximally different.
	 */
	public static double distJS2(double[] p, double[] q) {
		return distJSe(p, q) / LOG2;
	}

	/** Jensen-Shannon divergence, base e */
	public static double distJSe(double[] p, double[] q) {
		double dist = 0.0;
		for (int i = 0; i < p.length; i++) {
			double m2 = p[i] + q[i];
			if (m2 > 0.0) {
				// Exclude components equal to 0 from the computation.
				if ((p[i] == 0.0) || (q[i] == 0.0)) {
					dist += m2 * LOG2;
				} else {
					dist += (p[i] * (LOG2 + FastMath.log(p[i] / m2))) +
						(q[i] * (LOG2 + FastMath.log(q[i] / m2)));
				}
			}
		}
		return dist / 2;
	}

	public static double simCosine(double[] p, double[] q) {
		double dotProductPQ = 0.0, magnitudeP = 0.0, magnitudeQ = 0.0;
		for (int i = 0; i < p.length; i++) {
			dotProductPQ += p[i] * q[i];
			magnitudeP += p[i] * p[i];
			magnitudeQ += q[i] * q[i];
		}
		return dotProductPQ / (Math.sqrt(magnitudeP) * Math.sqrt(magnitudeQ));
	}

	/**
	 * @return 0 if index of the largest component is the same in p and q, 1 otherwise
	 */
	public static double distZeroOne(double[] p, double[] q) {
		int maxPIdx = 0, maxQIdx = 0;
		double maxP = 0.0, maxQ = 0.0;
		for (int i = 0; i < p.length; i++) {
			if (p[i] > maxP) {
				maxP = p[i];
				maxPIdx = i;
			}
			if (q[i] > maxQ) {
				maxQ = q[i];
				maxQIdx = i;
			}
		}
		return (maxPIdx == maxQIdx) ? 0.0 : 1.0;
	}

	/**
	 * Normalizes a vector so that its entries sum up to 1.
	 * If the sum of all entries is equal to 0.0, the vector is set to the uniform distribution.
	 */
	public static void normalize(double[] vector) {
		double sum = 0.0;
		for (int i = 0; i < vector.length; i++)
			sum += vector[i];
		if (sum != 0.0) {
			for (int i = 0; i < vector.length; i++)
				vector[i] /= sum;
		} else
			Arrays.fill(vector, 1.0 / vector.length);
	}

	public static double[] cumulative(double[] p) {
		double[] cumul = new double[p.length];
		cumul[0] = p[0];
		for (int i = 1; i < p.length; i++)
			cumul[i] = cumul[i - 1] + p[i];
		return cumul;
	}

	public static double entropy(double[] p) {
		double entropy = 0.0;
		for (double value : p)
			if (value != 0.0)
				entropy += value * (Math.log(value) / LOG2);
		return -entropy;
	}

	public static double perplexity(double logLikelihood, int numSamples) {
		return Math.exp(-(logLikelihood / numSamples));
	}

}
