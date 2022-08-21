package edu.tum.cs.influence.micro;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.moment.Mean;

import com.carrotsearch.hppc.sorting.IndirectComparator;
import com.carrotsearch.hppc.sorting.IndirectSort;

import edu.tum.cs.math.FastMath;
import edu.tum.cs.math.SparseVectorAutoregressiveModel;

public class GrangerCausality {

	private static final double zeroThreshold = 1E-7;

	public static class PreprocessingStatistics {
		public int numZeroRows, totalRows;
		public int numZeroesReplaced, totalValues;
		public boolean allRowsEliminated;
	}

	private double[] testStatistics;
	private final double coeffRatio;
	private final boolean isSignificant;
	private final double[] srcAbsCoeff, dstAbsCoeff;
	private final PreprocessingStatistics preprocStatistics;

	public GrangerCausality() {
		this(new PreprocessingStatistics());
	}

	private GrangerCausality(PreprocessingStatistics stats) {
		this(new double[0], 0.0, 1.0, null, null, stats);
	}

	protected GrangerCausality(double[] testStatistics, double coeffRatio, double alpha, double[] srcAbsCoeff,
			double[] dstAbsCoeff, PreprocessingStatistics preprocStatistics) {
		this.testStatistics = testStatistics;
		this.coeffRatio = coeffRatio;
		isSignificant = isSignificantJJ(alpha);
		this.srcAbsCoeff = srcAbsCoeff;
		this.dstAbsCoeff = dstAbsCoeff;
		this.preprocStatistics = preprocStatistics;
	}

	private static double[][] removeTimeConstantComponents(double[][] data, int[] map, PreprocessingStatistics stats) {
		int numRemaining = data.length;
		stats.totalRows += numRemaining;
		for (int i = 0; i < data.length; i++) {
			boolean isConst = true;
			for (int j = 0; j < data[i].length; j++) {
				if (Math.abs(data[i][j]) > zeroThreshold) {
					isConst = false;
					break;
				}
			}
			if (isConst) {
				data[i] = null;
				numRemaining--;
				stats.numZeroRows++;
			}
		}

		double[][] compactData = new double[numRemaining][];
		int idx = 0;
		for (int i = 0; i < data.length; i++) {
			if (data[i] != null) {
				compactData[idx] = data[i];
				if (map != null)
					map[idx] = i;
				idx++;
			}
		}
		return compactData;
	}

	private static double[][] additiveLogratioTransform(double[][] data, PreprocessingStatistics stats) {
		int d = data.length;
		int n = data[0].length;
		double[][] tfData = new double[d - 1][n];
		for (int i = 0; i < n; i++) {
			// zero replacement according to Martin-Fernandez et al., 2000
			int numZeroes = 0;
			for (int j = 0; j < d; j++) {
				if (data[j][i] == 0.0)
					numZeroes++;
			}
			double repl = (zeroThreshold * (numZeroes + 1) * (d - numZeroes)) / (d * d);

			for (int j = 0; j < d; j++) {
				if (data[j][i] == 0.0) {
					data[j][i] = repl;
					stats.numZeroesReplaced++;
				} else
					data[j][i] *= (1.0 - (numZeroes * repl));
				stats.totalValues++;
			}

			// actual transformation
			double denom = FastMath.log(data[d - 1][i]);
			for (int j = 0; j < (d - 1); j++)
				tfData[j][i] = FastMath.log(data[j][i]) - denom;
		}
		return tfData;
	}

	private static void meanCenter(double[] v) {
		Mean m = new Mean();
		double mean = m.evaluate(v);
		for (int i = 0; i < v.length; i++)
			v[i] -= mean;
	}

	private static double[][] concatColumns(double[][] m1, double[][] m2) {
		double[][] m = new double[m1.length + m2.length][];
		int idx = 0;
		for (double[] r : m1)
			m[idx++] = r;
		for (double[] r : m2)
			m[idx++] = r;
		return m;
	}

	private static double[][] preprocessData(double[][] data, boolean transform, int[] map,
			PreprocessingStatistics stats) {
		if (transform) {
			// Remove all time-constant zero components, as they do not contribute to the VAR, and the logratio tf.
			// does not work in the presence of zeroes.
			data = removeTimeConstantComponents(data, map, stats);
			data = additiveLogratioTransform(data, stats);
		}

		// Mean-center each component, then remove all remaining time-constant components that do not contribute to
		// the VAR.
		for (double[] row : data)
			meanCenter(row);
		data = removeTimeConstantComponents(data, map, stats);	// just checks for zeroes, so repeat after centering
		return data;
	}

	public static GrangerCausality compute(double[][] srcData, double[][] dstData, double alpha, int numDim,
			boolean transform, boolean recoverCoeff) {
		int[] srcMap = null, dstMap = null;
		if (recoverCoeff) {
			srcMap = new int[numDim];
			dstMap = new int[numDim];
		}
		PreprocessingStatistics stats = new PreprocessingStatistics();
		srcData = preprocessData(srcData, transform, srcMap, stats);
		dstData = preprocessData(dstData, transform, dstMap, stats);

		// Fit a VAR(1) model to the column-wise concatenation of influencee's and influencer's time series data,
		// but only compute the coefficients corresponding to the influencee's data.
		int srcLen = srcData.length;
		int dstLen = dstData.length;
		if ((srcLen == 0) || (dstLen == 0)) {	// not enough information to determine GC
			stats.allRowsEliminated = true;
			return new GrangerCausality(stats);
		}

		double[][] dstSrcDataArray = concatColumns(dstData, srcData);
		SparseVectorAutoregressiveModel var = SparseVectorAutoregressiveModel.fit(dstSrcDataArray, dstLen, true);

		// Save test statistics of all relevant coefficients to allow for significance testing at different
		// thresholds via the Javanmard-Javadi procedure.
		RealMatrix coeff = var.getCoefficients();
		RealMatrix testStat = var.getTestStatistics();
		int idx = 0;
		double[] aggTestStatistics = new double[dstLen * srcLen];
		double[] srcAbsCoeff = null, dstAbsCoeff = null;
		if (srcMap != null)
			srcAbsCoeff = new double[srcMap.length];
		if (dstMap != null)
			dstAbsCoeff = new double[dstMap.length];
		double coeffSrc = 0.0, coeffDst = 0.0;
		for (int i = 0; i < dstLen; i++) {
			for (int j = 0; j < (dstLen + srcLen); j++) {
				double absCoeff = Math.abs(coeff.getEntry(i, j));
				if (j >= dstLen) {
					coeffSrc += absCoeff;
					if (srcAbsCoeff != null)
						srcAbsCoeff[srcMap[j - dstLen]] += absCoeff;
					aggTestStatistics[idx++] = testStat.getEntry(i, j);
				} else {
					coeffDst += absCoeff;
					if (dstAbsCoeff != null)
						dstAbsCoeff[dstMap[j]] += absCoeff;
				}
			}
		}

		double coeffRatio = 0.0;
		if ((coeffSrc > 0.0) || (coeffDst > 0.0))
			coeffRatio = coeffSrc / (coeffSrc + coeffDst);
		return new GrangerCausality(aggTestStatistics, coeffRatio, alpha, srcAbsCoeff, dstAbsCoeff, stats);
	}

	public boolean isSignificantJJ() {
		return isSignificant;
	}

	/**
	 * Returns true if at least one significant test statistic remains after controlling FDR using the procedure of
	 * Javanmard and Javadi, 2018.
	 */
	public boolean isSignificantJJ(double alpha) {
		if (isEmpty())
			return false;
		Arrays.sort(testStatistics);
		return testStatistics[testStatistics.length - 1] > computeThresholdJJ(testStatistics, alpha);
	}

	private static double computeThresholdJJ(double[] z, double alpha) {
		int d = z.length;
		NormalDistribution norm = new NormalDistribution();
		double maxThreshold = 2.0 * FastMath.log(d);
		double maxY = Math.sqrt(maxThreshold - 2.0 * FastMath.log(FastMath.log(d)));
		int i = 0;
		while (i < d) {
			while (((i + 1) < d) && (z[i] == z[i + 1]))
				i++;
			double y = norm.inverseCumulativeProbability(1.0 - ((alpha * (d - i)) / (2.0 * d)));
			if (y > maxY)
				break;
			if (y <= z[i])
				return y;
			i++;
		}
		return maxThreshold;
	}

	public static <K> void testGroupSignificanceJJ(Map<K, GrangerCausality> edges, double alpha) {
		int d = 0;
		for (GrangerCausality gc : edges.values())
			d += gc.testStatistics.length;

		double[] testStatistics = new double[d];
		List<K> assignment = new ArrayList<K>(d);
		int offset = 0;
		for (Map.Entry<K, GrangerCausality> e : edges.entrySet()) {
			double[] z = e.getValue().testStatistics;
			System.arraycopy(z, 0, testStatistics, offset, z.length);
			for (int i = 0; i < z.length; i++)
				assignment.add(e.getKey());
			offset += z.length;
		}

		int[] perm = IndirectSort.mergesort(0, d, new IndirectComparator.AscendingDoubleComparator(testStatistics));
		double[] testStatisticsSorted = new double[d];
		for (int i = 0; i < perm.length; i++)
			testStatisticsSorted[i] = testStatistics[perm[i]];
		testStatistics = testStatisticsSorted;

		double threshold = GrangerCausality.computeThresholdJJ(testStatistics, alpha);
		Set<K> sig = new HashSet<K>();
		for (int i = testStatistics.length - 1; i >= 0; i--) {
			if (testStatistics[i] < threshold)
				break;
			sig.add(assignment.get(perm[i]));
		}
		edges.keySet().retainAll(sig);
	}

	public boolean isEmpty() {
		return (testStatistics.length == 0);
	}

	public double getMagnitude() {
		return coeffRatio;
	}

	public double[] getSrcAbsCoeff() {
		return srcAbsCoeff;
	}

	public double[] getDstAbsCoeff() {
		return dstAbsCoeff;
	}

	public PreprocessingStatistics getPreprocessingStatistics() {
		return preprocStatistics;
	}

}
