package edu.tum.cs.math;

import java.util.logging.Logger;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.rank.Median;

public class SparseVectorAutoregressiveModel {

	private static final Logger logger = Logger.getLogger(SparseVectorAutoregressiveModel.class.getName());
	private final RealMatrix coeff, t, p;

	protected SparseVectorAutoregressiveModel(RealMatrix coeff, RealMatrix t, RealMatrix p) {
		this.coeff = coeff;
		this.t = t;
		this.p = p;
	}

	static class DebiasedLassoHypothesisTest {
		private static final double debiasResol = 1.3;
		private static final int debiasMaxIter = 32;
		private static final double debiasEps = 1E-2;

		private final Array2DRowRealMatrix x;
		private final RealMatrix mxt;
		private final RealVector xSdInv, a;
		private final int n, d;
		private final double lambda;
		private final SqrtLasso sqrtLasso;
		private final NormalDistribution stdNormal = new NormalDistribution();

		private ArrayRealVector prevBeta, beta, vs, accDelta;

		public DebiasedLassoHypothesisTest(double[][] xData, boolean useFixedMu) {
			n = xData.length;
			d = xData[0].length;

			// The reference implementation of Javanmard and Montanari, 2014, deviates from the choice of lambda as
			// sqrt(log(d)/n), proposed by Li et al., 2018, and earlier papers on the SQRT-Lasso. This might be in
			// reference to Sun and Zhang, 2013 or Belloni et al., 2014, who derive a tighter upper bound on the optimal
			// lambda. Similar reasoning probably applies to the initial value of mu, computed below.
			lambda = Math.sqrt(stdNormal.inverseCumulativeProbability(1.0 - (0.1 / d)) / n);

			// column-wise mean centering (opt. coeff. are the same for centered and uncentered data)
			for (int i = 0; i < d; i++) {
				double mean = 0.0;
				for (int j = 0; j < n; j++)
					mean += xData[j][i];
				mean /= n;
				for (int j = 0; j < n; j++)
					xData[j][i] -= mean;
			}

			// transform to unit covariance
			// Tibshirani et al., 2012: "we typically would not standardize if the features were measured in the same
			// units", but the reference implementation of Javanmard and Montanari does it unconditionally, and it is
			// unclear to what extent their algorithm depends on it.
			xSdInv = new ArrayRealVector(d);
			for (int i = 0; i < d; i++) {
				double sdInv = 0.0;
				for (int j = 0; j < n; j++) {
					double v = xData[j][i];
					sdInv += v * v;
				}
				if (sdInv > 0.0) {
					sdInv = 1.0 / Math.sqrt(sdInv / n);
					for (int j = 0; j < n; j++)
						xData[j][i] *= sdInv;
				}
				xSdInv.setEntry(i, sdInv);
			}

			this.x = new Array2DRowRealMatrix(xData, false);
			sqrtLasso = new SqrtLasso(x);

			// optimize de-biasing matrix
			Array2DRowRealMatrix covEst = InPlaceLinAlg.multiply((Array2DRowRealMatrix) sqrtLasso.xt.multiply(x),+
					1.0 / n);
			RealVector covDiag = InPlaceLinAlg.getDiag(covEst, new ArrayRealVector(d));
			Array2DRowRealMatrix m;
			if (covDiag.getMinValue() > 0.0) {
				Array2DRowRealMatrix covOffDiag = InPlaceLinAlg.setDiag((Array2DRowRealMatrix) covEst.copy(), 0.0);
				double initMu = stdNormal.inverseCumulativeProbability(1.0 - (0.1 / (d * d))) * sqrtLasso.f;
				prevBeta = new ArrayRealVector(d);
				beta = new ArrayRealVector(d);
				vs = new ArrayRealVector(d);
				accDelta = new ArrayRealVector(d);
				m = optimizeDebiasingMatrix(covOffDiag, covDiag, initMu, useFixedMu);
			} else {
				// Could set m = I_d as per algorithm 1, step 5 of Javanmard and Montanari, 2014, but then a = covEst
				// and the result of the hypothesis test is undefined.
				throw new RuntimeException("got 0 on diagonal of estimated covariance matrix");
			}
			mxt = m.multiply(sqrtLasso.xt);

			// prepare for estimation of noise level
			Array2DRowRealMatrix mCov = (Array2DRowRealMatrix) m.multiply(covEst);
			a = new ArrayRealVector(d);	// a = diag(mCov * m^T)
			for (int i = 0; i < d; i++)
				a.setEntry(i, InPlaceLinAlg.dotProductOfRows(mCov, m, i));
		}

		private int optimizeDebiasingMatrixRow(Array2DRowRealMatrix covOffDiag, RealVector covDiag, int rowIdx,
				double mu, double mu0, RealVector beta) {
			beta.set(0.0);
			beta.setEntry(rowIdx, (1.0 - mu0) / covDiag.getEntry(rowIdx));
			if (mu >= mu0)
				return 0;

			int d = covDiag.getDimension();
			vs = (ArrayRealVector) InPlaceLinAlg.operate(covOffDiag, beta, vs).mapMultiplyToSelf(-1.0);
			double deltaNorm2 = 1.0;
			double prevNorm2 = 1.0;
			int numIter = 1;
			int prevNumIter = 1;
			accDelta.set(0.0);
			while ((numIter <= debiasMaxIter) && (deltaNorm2 >= (debiasEps * prevNorm2))) {
				for (int i = 0; i < d; i++) {
					double prevVal = beta.getEntry(i);
					double v = vs.getEntry(i);
					if (i == rowIdx)
						v += 1.0;
					double newVal = SqrtLasso.softThreshold(v, mu) / covDiag.getEntry(i);
					beta.setEntry(i, newVal);

					double delta = prevVal - newVal;
					vs = (ArrayRealVector) InPlaceLinAlg.addScaledCol(vs, covOffDiag, i, delta, vs);
					accDelta.addToEntry(i, -delta);
				}

				if (++numIter >= (2 * prevNumIter)) {
					deltaNorm2 = accDelta.getNorm();
					accDelta.set(0.0);
					prevNorm2 = beta.getNorm();
					prevNumIter = numIter;
					if (numIter > 10)
						vs = (ArrayRealVector) InPlaceLinAlg.operate(covOffDiag, beta, vs).mapMultiplyToSelf(-1.0);
				}
			}
			return numIter;
		}

		/**
		 * Straightforward re-implementation of the undocumented R routine of Javadi et al. (reference implementation of
		 * Javanmard and Montanari, 2014), which finds a suitable value of mu and solves the constrained convex
		 * optimization problem that determines the de-biasing matrix.
		 */
		private Array2DRowRealMatrix optimizeDebiasingMatrix(Array2DRowRealMatrix covOffDiag, RealVector covDiag,
				double initMu, boolean useFixedMu) {
			int d = covDiag.getDimension();
			double[] mu0 = new double[d];
			for (int i = 0; i < d; i++) {
				mu0[i] = InPlaceLinAlg.normLInfOfRow(covOffDiag, i);
				mu0[i] /= covDiag.getEntry(i);
				mu0[i] /= 1.0 + mu0[i];
			}

			Array2DRowRealMatrix m = new Array2DRowRealMatrix(d, d);
			for (int i = 0; i < d; i++) {
				double mu = initMu;
				boolean muStop = false;
				int numAttempt = 1;
				boolean incr = false;
				while (!muStop && (numAttempt < 10)) {
					ArrayRealVector tmp = beta;
					beta = prevBeta;
					prevBeta = tmp;
					int numIter = optimizeDebiasingMatrixRow(covOffDiag, covDiag, i, mu, mu0[i], beta);
					if (useFixedMu)
						break;

					if (numAttempt == 1) {
						if (numIter > debiasMaxIter) {
							incr = true;
							mu *= debiasResol;
						} else {
							incr = false;
							mu /= debiasResol;
						}
					} else {
						if (incr) {
							if (numIter > debiasMaxIter)
								mu *= debiasResol;
							else
								muStop = true;
						} else {
							if (numIter > debiasMaxIter) {
								mu *= debiasResol;
								tmp = beta;
								beta = prevBeta;
								prevBeta = tmp;
								muStop = true;
							} else
								mu /= debiasResol;
						}
					}
					numAttempt++;
				}
				m.setRow(i, beta.getDataRef());	// faster than setRowVector, which performs a per-element bounds check
			}
			return m;
		}

		private static double medianAbsoluteDeviation(double[] data) {
			Median m = new Median();
			double median = m.evaluate(data);
			double[] adjData = new double[data.length];
			for (int i = 0; i < data.length; i++)
				adjData[i] = Math.abs(data[i] - median);
			median = m.evaluate(adjData);
			return 1.4826 * median;
		}

		// see Javanmard and Montanari, 2014, section 3 (formula 17)
		private static double estimateNoiseSd(RealVector coeff, RealVector a, int n, double f) {
			double[] yNorm = new double[coeff.getDimension()];
			for (int i = 0; i < yNorm.length; i++)
				yNorm[i] = (coeff.getEntry(i) / Math.sqrt(a.getEntry(i))) / f;
			double sd0 = medianAbsoluteDeviation(yNorm);

			double y2Norm = 0.0;
			double aTrace = 0.0;
			for (int i = 0; i < yNorm.length; i++) {
				if (Math.abs(yNorm[i]) <= (3 * sd0)) {
					y2Norm += coeff.getEntry(i) * coeff.getEntry(i);
					aTrace += a.getEntry(i);
				}
			}
			if (aTrace == 0.0) {
				logger.fine("estimation of noise SD failed, falling back to " + sd0);
				return sd0;
			}
			return Math.sqrt(n * (y2Norm / aTrace));
		}

		public RealVector[] performTest(RealVector y) {
			// mean-centering of response vector
			double mean = 0.0;
			for (int i = 0; i < y.getDimension(); i++)
				mean += y.getEntry(i);
			mean /= y.getDimension();
			for (int i = 0; i < y.getDimension(); i++)
				y.setEntry(i, y.getEntry(i) - mean);

			// estimate coefficients, apply de-biasing matrix, estimate noise level
			RealVector coeff = sqrtLasso.minimize(y, lambda);
			RealVector coeffUnbiased = coeff.add(mxt.operate(y.subtract(x.operate(coeff))).mapDivideToSelf(n));
			double noiseSd = estimateNoiseSd(coeffUnbiased, a, n, sqrtLasso.f);
			if (noiseSd == 0.0)
				throw new RuntimeException("error estimating noise SD");

			// undo unit covariance transformation
			coeff = coeff.ebeMultiply(xSdInv);
			coeffUnbiased = coeffUnbiased.ebeMultiply(xSdInv);

			// compute test statistic Z and p-values
			RealVector z = new ArrayRealVector(d);
			RealVector p = new ArrayRealVector(d);
			for (int i = 0; i < d; i++) {
				double t = 0.0;
				if ((coeffUnbiased.getEntry(i) != 0.0) &&	// implies xSdInv > 0
					(a.getEntry(i) != 0.0)) {
					t = Math.abs(coeffUnbiased.getEntry(i)) /
						(sqrtLasso.f * noiseSd * xSdInv.getEntry(i) * Math.sqrt(a.getEntry(i)));
				}
				z.setEntry(i, t);
				p.setEntry(i, 2.0 * (1.0 - stdNormal.cumulativeProbability(t)));
			}
			return new RealVector[] { coeff, z, p };
		}
	}

	public static SparseVectorAutoregressiveModel fit(double[][] obs, boolean useFixedMu) {
		return fit(obs, obs.length, useFixedMu);
	}

	/**
	 * Express sparse VAR as Lasso optimization problem and solve with debiased Lasso of Javanmard and Montanari,
	 * 2014. Assumes mean-centered observations.
	 */
	public static SparseVectorAutoregressiveModel fit(double[][] obs, int maxDim, boolean useFixedMu) {
		int d = obs.length;
		int n = obs[0].length;
		if (n > (2 * d))	// check n << d to see if Sqrt-Lasso is still appropriate
			logger.warning("d = " + d + ", n = " + n);

		// Express VAR as a series of multiple linear regression problems, one for each component of the response
		// variable: Given observation vectors o0 .. oN, build matrix X = (o0 .. oN-1)^T. Perform Lasso regression
		// on X and each vector yi = (o1i .. oNi)^T.
		double[][] xData = new double[n - 1][d];
		for (int i = 0; i < (n - 1); i++)
			for (int j = 0; j < d; j++)
				xData[i][j] = obs[j][i];
		DebiasedLassoHypothesisTest test = new DebiasedLassoHypothesisTest(xData, useFixedMu);

		int rd = Math.min(d, maxDim);
		RealMatrix coeff = new Array2DRowRealMatrix(rd, d);
		RealMatrix z = new Array2DRowRealMatrix(rd, d);
		RealMatrix p = new Array2DRowRealMatrix(rd, d);
		for (int i = 0; i < rd; i++) {
			double[] yData = new double[n - 1];
			System.arraycopy(obs[i], 1, yData, 0, n - 1);
			RealVector y = new ArrayRealVector(yData, false);

			RealVector[] res = test.performTest(y);
			coeff.setRowVector(i, res[0]);
			z.setRowVector(i, res[1]);
			p.setRowVector(i, res[2]);
		}
		return new SparseVectorAutoregressiveModel(coeff, z, p);
	}

	/**
	 * Generate a VAR with random coefficients; partial re-implementation of the R routine sugm.generator by Li et al.
	 * (package "flare") with graph pattern "random".
	 */
	public static SparseVectorAutoregressiveModel random(RandomGenerator rand, int d, double pNonZero, double rho) {
		if (pNonZero < 0.5)
			pNonZero = Math.sqrt(0.5 * pNonZero);
		else
			pNonZero = 1.0 - Math.sqrt(0.5 - (0.5 * pNonZero));

		RealMatrix omega = new Array2DRowRealMatrix(d, d);
		for (int i = 0; i < d; i++) {
			for (int j = 0; j < d; j++)
				omega.setEntry(i, j, rand.nextDouble() * 0.5);
		}
		omega = omega.add(omega.transpose());
		for (int i = 0; i < d; i++) {
			for (int j = 0; j < d; j++) {
				if ((i == j) || (omega.getEntry(i, j) >= pNonZero))
					omega.setEntry(i, j, 0.0);
				else
					omega.setEntry(i, j, 1.0);
			}
		}
		omega = omega.scalarMultiply(0.3);

		EigenDecomposition eigen = new EigenDecomposition(omega);
		double diag = Math.abs(StatUtils.min(eigen.getRealEigenvalues())) + 0.2;
		for (int i = 0; i < d; i++)
			omega.setEntry(i, i, diag);
		LUDecomposition lu = new LUDecomposition(omega);
		PearsonsCorrelation cor = new PearsonsCorrelation();
		RealMatrix sigma = cor.covarianceToCorrelation(lu.getSolver().getInverse());
		lu = new LUDecomposition(sigma);
		omega = lu.getSolver().getInverse();

		eigen = new EigenDecomposition(omega);
		double maxEigen = 0.0;
		for (int i = 0; i < d; i++) {
			double e = Math.abs(eigen.getRealEigenvalue(i));
			if (e > maxEigen)
				maxEigen = e;
		}
		omega = omega.scalarMultiply(rho / maxEigen);

		return new SparseVectorAutoregressiveModel(omega, null, null);
	}

	/**
	 * @return coefficient matrix B in concise notation
	 * see https://en.wikipedia.org/wiki/Vector_autoregression#Concise_matrix_notation
	 */
	public RealMatrix getCoefficients() {
		return coeff;
	}

	public RealMatrix getTestStatistics() {
		return t;
	}

	/** @return p-values associated with the coefficients */
	public RealMatrix getPValues() {
		return p;
	}

	public RealMatrix generate(RandomGenerator rand, int n, int numBurnIn, double noiseCov) {
		int d = coeff.getRowDimension();

		// generate Gaussian noise
		double[][] cov = new double[d][d];
		for (int i = 0; i < d; i++)
			cov[i][i] = noiseCov;
		MultivariateNormalDistribution noiseDist = new MultivariateNormalDistribution(rand, new double[d], cov);
		RealMatrix noise = new Array2DRowRealMatrix(d, numBurnIn + n);
		for (int i = 0; i < noise.getColumnDimension(); i++)
			noise.setColumn(i, noiseDist.sample());

		// simulate VAR
		RealMatrix samples = new Array2DRowRealMatrix(d, n);
		RealVector v = new ArrayRealVector(d);
		for (int i = 0; i < (numBurnIn + n); i++) {
			v = coeff.operate(v).add(noise.getColumnVector(i));
			if (i >= numBurnIn)
				samples.setColumnVector(i - numBurnIn, v);
		}
		return samples;
	}

}
