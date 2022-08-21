package edu.tum.cs.math;

import java.util.logging.Logger;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class MultiDimensionalScaling {

	private static final Logger logger = Logger.getLogger(MultiDimensionalScaling.class.getName());
	private static final double tol = 1E-6;

	private final int n;
	private final SymmetricEigenDecomposition eigen;

	public MultiDimensionalScaling(double[][] d) {
		this(new Array2DRowRealMatrix(d, false));
	}

	public MultiDimensionalScaling(RealMatrix d) {
		n = d.getRowDimension();

		// build matrix of squared dissimilarities
		RealMatrix d2 = new Array2DRowRealMatrix(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				double dd = d.getEntry(i, j);
				d2.setEntry(i, j, dd * dd);
			}
		}

		// double centering
		RealMatrix j = (new Array2DRowRealMatrix(n, n)).scalarAdd(1.0 / n);
		RealMatrix c = MatrixUtils.createRealIdentityMatrix(n).subtract(j);
		d2 = c.multiply(d2).multiply(c).scalarMultiply(-0.5);

		// perform eigendecomposition and store eigenvectors/values for later use
		eigen = new SymmetricEigenDecomposition(d2, 50);

		// diagnose deviations from metricity in the dissimilarity matrix
		double sum = 0.0, sumPositive = 0.0;
		for (double v : eigen.getRealEigenvalues()) {
			if (v > 0.0)
				sumPositive += v;
			sum += Math.abs(v);
		}
		double r = sumPositive / sum;
		if (Math.abs(1.0 - r) > tol) {
			logger.warning("deviation from metricity: ratio of positive eigenvalues = " + r);
			logger.warning("largest positive EV = " + eigen.getRealEigenvalues()[0] + ", largest negative EV = " +
					eigen.getRealEigenvalues()[n - 1] + ", ratio of abs. EV = " +
					(Math.abs(eigen.getRealEigenvalues()[n - 1]) / eigen.getRealEigenvalues()[0]));
		}
	}

	public RealMatrix getEmbedding(int dim) {
		RealMatrix e = eigen.getEigenvectors().getSubMatrix(0, n - 1, 0, dim - 1);
		RealMatrix v = new Array2DRowRealMatrix(dim, dim);
		for (int i = 0; i < dim; i++)
			v.setEntry(i, i, Math.sqrt(eigen.getRealEigenvalues()[i]));
		return e.multiply(v);
	}

}
