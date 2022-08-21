package edu.tum.cs.time.hmm;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.rng.RestorableUniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationInteger;
import be.ac.ulg.montefiore.run.jahmm.OpdfInteger;
import edu.tum.cs.math.dist.DiscreteDistribution;

public class HmmDivergence {

	private static final RestorableUniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S);

	// values taken from NetworkX
	private static final int maxPowerIter = 100;
	private static final double tolPowerIter = 1E-6;

	// find the largest eigenvalue via power iteration
	public static double findLargestEigenvalue(RealMatrix m) {
		// randomly initialize candidate eigenvector b
		RealVector b = new ArrayRealVector(m.getColumnDimension());
		for (int i = 0; i < b.getDimension(); i++)
			b.setEntry(i, rand.nextDouble());

		// obtain the largest eigenvector by iterative approximation
		RealVector bPrev;
		int numIter = 0;
		do {
			bPrev = b;
			b = m.operate(b);
			b.mapDivideToSelf(b.getL1Norm());
		} while ((b.getL1Distance(bPrev) > (b.getDimension() * tolPowerIter)) && (++numIter < maxPowerIter));

		// compute Rayleigh coefficient to get the associated eigenvalue
		return b.dotProduct(m.operate(b)) / b.dotProduct(b);
	}

	private static double[][] getTransitionProbabilities(Hmm<?> hmm) {
		double[][] pTrans = new double[hmm.nbStates()][hmm.nbStates()];
		for (int i = 0; i < hmm.nbStates(); i++) {
			for (int j = 0; j < hmm.nbStates(); j++)
				pTrans[i][j] = hmm.getAij(i, j);
		}
		return pTrans;
	}

	private static double[][] getOutputProbabilities(Hmm<? extends ObservationInteger> hmm) {
		double[][] pOut = new double[hmm.nbStates()][];
		for (int i = 0; i < hmm.nbStates(); i++) {
			OpdfInteger pdf = (OpdfInteger) hmm.getOpdf(i);
			pOut[i] = new double[pdf.nbEntries()];
			for (int j = 0; j < pdf.nbEntries(); j++)
				pOut[i][j] = pdf.probability(new ObservationInteger(j));
		}
		return pOut;
	}

	public static double divergenceRate(Hmm<? extends ObservationInteger> hmm1,
			Hmm<? extends ObservationInteger> hmm2) {
		return divergenceRate(getTransitionProbabilities(hmm1), getOutputProbabilities(hmm1),
				getTransitionProbabilities(hmm2), getOutputProbabilities(hmm2));
	}

	private static RealMatrix buildTransferMatrix(double[][] pTrans1, double[][] pOut1, double[][] pTrans2,
			double[][] pOut2) {
		int n1 = pTrans1.length;
		int n2 = pTrans2.length;
		int nOut = pOut1[0].length;

		double[][] mm = new double[n1 * n2][n1 * n2];
		for (int i = 0; i < n1; i++) {
			for (int j = 0; j < n1; j++) {
				for (int m = 0; m < n2; m++) {
					for (int n = 0; n < n2; n++) {
						int r = ((i + 1) * (m + 1)) - 1;
						int c = ((j + 1) * (n + 1)) - 1;
						for (int o = 0; o < nOut; o++)
							mm[r][c] += pTrans1[i][j] * pOut1[i][o] * pTrans2[m][n] * pOut2[m][o];
					}
				}
			}
		}
		return new Array2DRowRealMatrix(mm, false);
	}

	// implementation of Yang et al., 2020, "Measures of distinguishability between stochastic processes"
	public static double divergenceRate(double[][] pTrans1, double[][] pOut1, double[][] pTrans2, double[][] pOut2) {
		if (pOut1[0].length != pOut2[0].length)
			throw new RuntimeException("HMMs must have the same vocabulary");

		double muPq = findLargestEigenvalue(buildTransferMatrix(pTrans1, pOut1, pTrans2, pOut2));
		double muPp = findLargestEigenvalue(buildTransferMatrix(pTrans1, pOut1, pTrans1, pOut1));
		double muQq = findLargestEigenvalue(buildTransferMatrix(pTrans2, pOut2, pTrans2, pOut2));
		return -0.5 * Math.log(muPq / (Math.sqrt(muPp) * Math.sqrt(muQq))) / Math.log(2.0);
	}

	public static double distTransitionMatrix(double[][] pTrans1, double[][] pTrans2) {
		double sumDist = 0.0;
		for (int i = 0; i < pTrans1.length; i++)
			sumDist += Math.sqrt(DiscreteDistribution.distJS2(pTrans1[i], pTrans2[i]));
		return sumDist / pTrans1.length;
	}

}
