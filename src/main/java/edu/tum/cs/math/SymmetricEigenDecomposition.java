package edu.tum.cs.math;

import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.Precision;

/**
 * Lifted from Apache Commons Math with minor modifications: assume that input matrix is symmetric, make maximum number
 * of iterations configurable.
 */
public class SymmetricEigenDecomposition {

	private static final int defaultMaxIter = 30;

	private final int maxIter;
	private double[] main;
	private double[] secondary;
	private double[] realEigenvalues;
	private RealMatrix eigenvectors;

	public SymmetricEigenDecomposition(RealMatrix matrix) {
		this(matrix, defaultMaxIter);
	}

	public SymmetricEigenDecomposition(RealMatrix matrix, int maxIter) {
		this.maxIter = maxIter;
		double[][] q = transformToTridiagonal(matrix);
		findEigenVectors(q);
	}

	private double[][] transformToTridiagonal(RealMatrix matrix) {
		final double[][] householderVectors = matrix.getData();
		final int m = householderVectors.length;
		main = new double[m];
		secondary = new double[m - 1];

		final double[] z = new double[m];
		for (int k = 0; k < m - 1; k++) {

			//zero-out a row and a column simultaneously
			final double[] hK = householderVectors[k];
			main[k] = hK[k];
			double xNormSqr = 0;
			for (int j = k + 1; j < m; ++j) {
				final double c = hK[j];
				xNormSqr += c * c;
			}
			final double a = (hK[k + 1] > 0) ? -Math.sqrt(xNormSqr) : Math.sqrt(xNormSqr);
			secondary[k] = a;
			if (a != 0.0) {
				// apply Householder transform from left and right simultaneously

				hK[k + 1] -= a;
				final double beta = -1 / (a * hK[k + 1]);

				// compute a = beta A v, where v is the Householder vector
				// this loop is written in such a way
				//   1) only the upper triangular part of the matrix is accessed
				//   2) access is cache-friendly for a matrix stored in rows
				Arrays.fill(z, k + 1, m, 0);
				for (int i = k + 1; i < m; ++i) {
					final double[] hI = householderVectors[i];
					final double hKI = hK[i];
					double zI = hI[i] * hKI;
					for (int j = i + 1; j < m; ++j) {
						final double hIJ = hI[j];
						zI   += hIJ * hK[j];
						z[j] += hIJ * hKI;
					}
					z[i] = beta * (z[i] + zI);
				}

				// compute gamma = beta vT z / 2
						double gamma = 0;
				for (int i = k + 1; i < m; ++i) {
					gamma += z[i] * hK[i];
				}
				gamma *= beta / 2;

				// compute z = z - gamma v
				for (int i = k + 1; i < m; ++i) {
					z[i] -= gamma * hK[i];
				}

				// update matrix: A = A - v zT - z vT
				// only the upper triangular part of the matrix is updated
				for (int i = k + 1; i < m; ++i) {
					final double[] hI = householderVectors[i];
					for (int j = i; j < m; ++j) {
						hI[j] -= hK[i] * z[j] + z[i] * hK[j];
					}
				}
			}
		}
		main[m - 1] = householderVectors[m - 1][m - 1];

		// generate Q: build up first part of the matrix by applying Householder transforms
		double[][] q = new double[m][m];
		for (int k = m - 1; k >= 1; --k) {
			final double[] hK = householderVectors[k - 1];
			q[k][k] = 1;
			if (hK[k] != 0.0) {
				final double inv = 1.0 / (secondary[k - 1] * hK[k]);
				double beta = 1.0 / secondary[k - 1];
				q[k][k] = 1 + beta * hK[k];
				for (int i = k + 1; i < m; ++i) {
					q[i][k] = beta * hK[i];
				}
				for (int j = k + 1; j < m; ++j) {
					beta = 0;
					for (int i = k + 1; i < m; ++i) {
						beta += q[i][j] * hK[i];
					}
					beta *= inv;
					q[k][j] = beta * hK[k];
					for (int i = k + 1; i < m; ++i) {
						q[i][j] += beta * hK[i];
					}
				}
			}
		}
		q[0][0] = 1;
		return q;
	}

	private void findEigenVectors(final double[][] z) {
		final int n = main.length;
		realEigenvalues = new double[n];
		final double[] e = new double[n];
		for (int i = 0; i < n - 1; i++) {
			realEigenvalues[i] = main[i];
			e[i] = secondary[i];
		}
		realEigenvalues[n - 1] = main[n - 1];
		e[n - 1] = 0;

		// Determine the largest main and secondary value in absolute term.
		double maxAbsoluteValue = 0;
		for (int i = 0; i < n; i++) {
			if (Math.abs(realEigenvalues[i]) > maxAbsoluteValue) {
				maxAbsoluteValue = Math.abs(realEigenvalues[i]);
			}
			if (Math.abs(e[i]) > maxAbsoluteValue) {
				maxAbsoluteValue = Math.abs(e[i]);
			}
		}
		// Make null any main and secondary value too small to be significant
		if (maxAbsoluteValue != 0) {
			for (int i=0; i < n; i++) {
				if (Math.abs(realEigenvalues[i]) <= Precision.EPSILON * maxAbsoluteValue) {
					realEigenvalues[i] = 0;
				}
				if (Math.abs(e[i]) <= Precision.EPSILON * maxAbsoluteValue) {
					e[i]=0;
				}
			}
		}

		for (int j = 0; j < n; j++) {
			int its = 0;
			int m;
			do {
				for (m = j; m < n - 1; m++) {
					double delta = Math.abs(realEigenvalues[m]) +
							Math.abs(realEigenvalues[m + 1]);
					if (Math.abs(e[m]) + delta == delta) {
						break;
					}
				}
				if (m != j) {
					if (its == maxIter) {
						throw new RuntimeException("maximum number of iterations");
					}
					its++;
					double q = (realEigenvalues[j + 1] - realEigenvalues[j]) / (2 * e[j]);
					double t = Math.sqrt(1 + q * q);
					if (q < 0.0) {
						q = realEigenvalues[m] - realEigenvalues[j] + e[j] / (q - t);
					} else {
						q = realEigenvalues[m] - realEigenvalues[j] + e[j] / (q + t);
					}
					double u = 0.0;
					double s = 1.0;
					double c = 1.0;
					int i;
					for (i = m - 1; i >= j; i--) {
						double p = s * e[i];
						double h = c * e[i];
						if (Math.abs(p) >= Math.abs(q)) {
							c = q / p;
							t = Math.sqrt(c * c + 1.0);
							e[i + 1] = p * t;
							s = 1.0 / t;
							c *= s;
						} else {
							s = p / q;
							t = Math.sqrt(s * s + 1.0);
							e[i + 1] = q * t;
							c = 1.0 / t;
							s *= c;
						}
						if (e[i + 1] == 0.0) {
							realEigenvalues[i + 1] -= u;
							e[m] = 0.0;
							break;
						}
						q = realEigenvalues[i + 1] - u;
						t = (realEigenvalues[i] - q) * s + 2.0 * c * h;
						u = s * t;
						realEigenvalues[i + 1] = q + u;
						q = c * t - h;
						for (int ia = 0; ia < n; ia++) {
							p = z[ia][i + 1];
							z[ia][i + 1] = s * z[ia][i] + c * p;
							z[ia][i] = c * z[ia][i] - s * p;
						}
					}
					if (t == 0.0 && i >= j) {
						continue;
					}
					realEigenvalues[j] -= u;
					e[j] = q;
					e[m] = 0.0;
				}
			} while (m != j);
		}

		//Sort the eigen values (and vectors) in increase order
		for (int i = 0; i < n; i++) {
			int k = i;
			double p = realEigenvalues[i];
			for (int j = i + 1; j < n; j++) {
				if (realEigenvalues[j] > p) {
					k = j;
					p = realEigenvalues[j];
				}
			}
			if (k != i) {
				realEigenvalues[k] = realEigenvalues[i];
				realEigenvalues[i] = p;
				for (int j = 0; j < n; j++) {
					p = z[j][i];
					z[j][i] = z[j][k];
					z[j][k] = p;
				}
			}
		}

		// Determine the largest eigen value in absolute term.
		maxAbsoluteValue = 0;
		for (int i = 0; i < n; i++) {
			if (Math.abs(realEigenvalues[i]) > maxAbsoluteValue) {
				maxAbsoluteValue = Math.abs(realEigenvalues[i]);
			}
		}
		// Make null any eigen value too small to be significant
		if (maxAbsoluteValue != 0.0) {
			for (int i=0; i < n; i++) {
				if (Math.abs(realEigenvalues[i]) < Precision.EPSILON * maxAbsoluteValue) {
					realEigenvalues[i] = 0;
				}
			}
		}

		eigenvectors = new Array2DRowRealMatrix(z);
	}

	public RealMatrix getEigenvectors() {
		return eigenvectors;
	}

	public double[] getRealEigenvalues() {
		return realEigenvalues;
	}

}
