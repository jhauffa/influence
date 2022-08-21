package edu.tum.cs.math;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;

import edu.tum.cs.math.dist.BetaDistribution;

public abstract class MultipleHypothesisTest<H, D> {

	private static final int initDelta = 10;
	private static final double a = 1.25;	// scale factor for delta
	private static final double eps = 0.01; 	// limit of spending sequence for Clopper-Pearson conf.int.
	private static final int r = 10000;	// inverse scale factor of spending sequence

	private final double alpha;
	private final int expMaxIter;
	private final boolean perHypoLimit;

	protected MultipleHypothesisTest(double alpha, int expMaxIter, boolean perHypoLimit) {
		this.alpha = alpha;
		this.expMaxIter = expMaxIter;
		this.perHypoLimit = perHypoLimit;
	}

	protected D[] preparePerIterationData(int n, int baseIter, List<H> hypos) {
		return null;
	}

	protected Map<H, D> preparePerHypothesisData(List<H> hypos, int maxIter) {
		return null;
	}

	/** generate n MC samples for the hypothesis and count exceedance of reference statistic */
	protected abstract int generateSamples(int n, H hypo, D[] perIterationData, D perHypothesisData);

	public double[] performTests(List<H> hypos, AtomicLong numTestsCompleted) {
		return performTests(hypos, numTestsCompleted, null);
	}

	public double[] performTests(List<H> hypos, AtomicLong numTestsCompleted, boolean[] insufficientData) {
		int n = hypos.size();
		int maxIter = expMaxIter * n;
		int delta = initDelta;

		boolean[] isActive = new boolean[n];
		Arrays.fill(isActive, true);
		boolean[] isRejected = new boolean[n];
		double[][] confInt = new double[2][n];
		Arrays.fill(confInt[1], 1.0);
		int[] numExceedance = new int[n];
		int[] numSamples = new int[n];
		Map<H, D> perHypothesisData = preparePerHypothesisData(hypos, maxIter);

		int numIter = 0, numOverallIter = 0;
		int numUnclassified = n;
		int numUnclassifiedPrev = n;
		while ((numUnclassified > 0) && (numOverallIter < maxIter)) {
			D[] perIterationData = preparePerIterationData(delta, numIter, hypos);
			for (int j = 0; j < n; j++) {
				if (!isActive[j] || isRejected[j])
					continue;

				// perform delta iterations of MC test
				H hypo = hypos.get(j);
				D curHypoData = null;
				if (perHypothesisData != null)
					curHypoData = perHypothesisData.get(hypo);
				int curNumExceedance = generateSamples(delta, hypo, perIterationData, curHypoData);
				if (curNumExceedance < 0) {
					confInt[0][j] = confInt[1][j] = 0.0;
					if (insufficientData != null)
						insufficientData[j] = true;
					continue;
				}
				numExceedance[j] += curNumExceedance;
				int numSamplesPrev = numSamples[j];
				numSamples[j] += delta;
				numOverallIter += delta;
				if (perHypoLimit && (numSamples[j] >= expMaxIter)) {
					confInt[0][j] = confInt[1][j] = 0.0;
					continue;
				}

				// update confidence intervals
				double l, u;
				double rho = (((numSamples[j] * eps) / (numSamples[j] + r)) -
						((numSamplesPrev * eps) / (numSamplesPrev + r))) / (2.0 * n);
				if (numExceedance[j] == 0) {
					l = 0.0;
					u = 1.0 - FastMath.pow(rho, 1.0 / numSamples[j]);
				} else if (numExceedance[j] == numSamples[j]) {
					l = FastMath.pow(rho, 1.0 / numSamples[j]);
					u = 1.0;
				} else {
					BetaDistribution beta = new BetaDistribution(numSamples[j] - numExceedance[j],
							numExceedance[j] + 1);
					l = 1.0 - beta.inverseCumulativeProbability(rho);
					beta = new BetaDistribution(numSamples[j] + 1 - numExceedance[j], numExceedance[j]);
					u = 1.0 - beta.inverseCumulativeProbability(1.0 - rho);
				}
				confInt[0][j] = Math.max(confInt[0][j], l);
				confInt[1][j] = Math.min(u, confInt[1][j]);
			}

			Arrays.fill(isActive, false);
			for (Integer idx : testSignificanceBH(confInt[0], alpha))
				isActive[idx] = true;
			Arrays.fill(isRejected, false);
			for (Integer idx : testSignificanceBH(confInt[1], alpha))
				isRejected[idx] = true;

			numUnclassifiedPrev = numUnclassified;
			numUnclassified = 0;
			for (int i = 0; i < n; i++) {
				if (isActive[i] && !isRejected[i])
					numUnclassified++;
			}
			numTestsCompleted.addAndGet(numUnclassifiedPrev - numUnclassified);

			numIter += delta;
			delta = Math.min((int) Math.floor(a * delta), 100000);
		}

		numTestsCompleted.addAndGet(numUnclassified);
		double[] p = new double[n];
		for (int i = 0; i < n; i++)
			p[i] = (double) (1 + numExceedance[i]) / (numSamples[i] + 1);
		return p;
	}

	public static Integer[] testSignificanceBH(final double[] p, double alpha) {
		Integer[] idx = new Integer[p.length];
		for (int i = 0; i < p.length; i++)
			idx[i] = i;
		Arrays.sort(idx, new Comparator<Integer>() {
			@Override
			public int compare(Integer i1, Integer i2) {
				return Double.compare(p[i1], p[i2]);
			}
		});

		int i;
		for (i = p.length - 1; i >= 0; i--) {
			if (p[idx[i]] <= (((double) (i + 1) / p.length) * alpha))
				return Arrays.copyOfRange(idx, 0, i + 1);
		}
		return new Integer[0];
	}

}
