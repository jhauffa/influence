package edu.tum.cs.math;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class SqrtLasso {

	private static final int numSteps = 5;
	private static final double gdEps = 1E-4;
	private static final int gdMaxIter = (int) 1E5;
	private static final int gdMaxLineSearchIter = gdMaxIter / 10;
	private static final double gdLineSearchEps = 1E-14;
	private static final double gdRate = 1.0;

	private final Array2DRowRealMatrix x;
	public final Array2DRowRealMatrix xt;
	public final double f;
	private RealVector gradPrev, paramDelta, pNew, err;

	public SqrtLasso(Array2DRowRealMatrix x) {
		this.x = x;
		xt = (Array2DRowRealMatrix) x.transpose();
		int n = x.getRowDimension();
		f = 1.0 / Math.sqrt(n);

		err = new ArrayRealVector(n);
		int d = x.getColumnDimension();
		gradPrev = new ArrayRealVector(d);
		paramDelta = new ArrayRealVector(d);
		pNew = new ArrayRealVector(d);
	}

	public static double softThreshold(double v, double l) {
		if (v >= l)
			return v - l;
		if (v <= -l)
			return v + l;
		return 0.0;
	}

	private static RealVector softThreshold(RealVector x, double l) {
		for (int i = 0; i < x.getDimension(); i++)
			x.setEntry(i, softThreshold(x.getEntry(i), l));
		return x;
	}

	/**
	 * Use proximal gradient descent with line search to solve the SQRT-Lasso optimization problem. Assume that x
	 * and y have a mean of zero and x has unit variance.
	 */
	private RealVector minimizeStep(RealVector y, double lambda, RealVector p) {
		double rate = gdRate;
		double paramDeltaNorm;
		int numIter = -1;
		do {
			if (++numIter >= gdMaxIter)
				throw new RuntimeException("SQRT-Lasso did not converge after maximum number of iterations");

			RealVector pPrev = p;
			InPlaceLinAlg.subtract(y, InPlaceLinAlg.operate(x, pPrev, err), err);
			double evalPrev = err.getNorm();
			gradPrev = InPlaceLinAlg.operate(xt, err, gradPrev).mapMultiplyToSelf(-f / evalPrev);

			int numLineSearchIter = -1;
			double approxDelta = 0.0;
			while (true) {
				if (++numLineSearchIter >= gdMaxLineSearchIter)
					throw new RuntimeException("stuck in line search: delta = " + approxDelta);

				p = softThreshold(InPlaceLinAlg.subtract(pPrev, InPlaceLinAlg.multiply(gradPrev, rate, pNew), pNew),
						lambda * rate);
				paramDelta = InPlaceLinAlg.subtract(p, pPrev, paramDelta);
				paramDeltaNorm = paramDelta.getNorm();

				double eval = InPlaceLinAlg.subtract(y, InPlaceLinAlg.operate(x, p, err), err).getNorm();
				approxDelta = (f * (eval - evalPrev)) - gradPrev.dotProduct(paramDelta) -
						(0.5 * rate * paramDeltaNorm * paramDeltaNorm);
				if (approxDelta <= gdLineSearchEps)
					break;
				rate *= 0.5;
			}
			pNew = pPrev;
		} while (paramDeltaNorm >= gdEps);
		return p;
	}

	/** Perform pathwise optimization according to Li et al., 2018 (algorithm 3). */
	public RealVector minimize(RealVector y, double lambda) {
		RealVector p = new ArrayRealVector(x.getColumnDimension());	// all parameters zero: solution at init. lambda
		double yNorm = y.getNorm();
		if (yNorm > 0.0) {
			double curLambda = InPlaceLinAlg.operate(xt, y, pNew).mapMultiplyToSelf(f / yNorm).getLInfNorm();
			double eta = Math.pow(lambda / curLambda, 1.0 / numSteps);
			for (int i = 0; i < numSteps; i++) {
				curLambda *= eta;
				p = minimizeStep(y, curLambda, p);
			}
		}
		return p;
	}

}
