package edu.tum.cs.math;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class InPlaceLinAlg {

	/** dest := a - b (dest may be one of a, b) */
	public static RealVector subtract(RealVector a, RealVector b, RealVector dest) {
		for (int i = 0; i < dest.getDimension(); i++)
			dest.setEntry(i, a.getEntry(i) - b.getEntry(i));
		return dest;
	}

	/** dest := c * a (dest may be a) */
	public static RealVector multiply(RealVector a, double c, RealVector dest) {
		for (int i = 0; i < dest.getDimension(); i++)
			dest.setEntry(i, c * a.getEntry(i));
		return dest;
	}

	/** m := c * m */
	public static Array2DRowRealMatrix multiply(Array2DRowRealMatrix m, double c) {
		double[][] mData = m.getDataRef();
		for (int i = 0; i < mData.length; i++) {
			for (int j = 0; j < mData[i].length; j++)
				mData[i][j] *= c;
		}
		return m;
	}

	/** dest := M * a (dest must not be a) */
	public static RealVector operate(Array2DRowRealMatrix m, RealVector a, RealVector dest) {
		double[][] mData = m.getDataRef();
		for (int i = 0; i < dest.getDimension(); i++) {
			double d = 0.0;
			for (int j = 0; j < a.getDimension(); j++)
				d += mData[i][j] * a.getEntry(j);
			dest.setEntry(i, d);
		}
		return dest;
	}

	public static RealVector getDiag(Array2DRowRealMatrix m, RealVector dest) {
		double[][] mData = m.getDataRef();
		for (int i = 0; i < mData.length; i++)
			dest.setEntry(i, mData[i][i]);
		return dest;
	}

	public static Array2DRowRealMatrix setDiag(Array2DRowRealMatrix m, double c) {
		double[][] mData = m.getDataRef();
		for (int i = 0; i < mData.length; i++)
			mData[i][i] = c;
		return m;
	}

	public static double dotProductOfRows(Array2DRowRealMatrix m1, Array2DRowRealMatrix m2, int row) {
		double[][] m1Data = m1.getDataRef();
		double[][] m2Data = m2.getDataRef();
		double v = 0.0;
		for (int i = 0; i < m1Data[row].length; i++)
			v += m1Data[row][i] * m2Data[row][i];
		return v;
	}

	public static double normLInfOfRow(Array2DRowRealMatrix m, int row) {
		double[][] mData = m.getDataRef();
		double norm = 0.0;
		for (int i = 0; i < mData[row].length; i++) {
			double v = Math.abs(mData[row][i]);
			if (v > norm)
				norm = v;
		}
		return norm;
	}

	/** dest := a + c * m_(*,col) (dest may be a) */
	public static RealVector addScaledCol(RealVector a, Array2DRowRealMatrix m, int col, double c, RealVector dest) {
		double[][] mData = m.getDataRef();
		for (int i = 0; i < mData.length; i++)
			dest.setEntry(i, a.getEntry(i) + (c * mData[i][col]));
		return dest;
	}

}
