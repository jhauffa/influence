package edu.tum.cs.time.hmm;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Test;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.HmmBuilder;
import be.ac.ulg.montefiore.run.jahmm.ObservationInteger;
import be.ac.ulg.montefiore.run.jahmm.OpdfInteger;

public class HmmDivergenceTest {

	private static final double eps = 1E-6;

	@Test
	public void testLargestEigenvalue() {
		RealMatrix m = new Array2DRowRealMatrix(new double[][] { { 2.0, 1.0 }, { 1.0, 2.0 } }, false);
		double v = HmmDivergence.findLargestEigenvalue(m);
		assertEquals(3.0, v, eps);
	}

	private static Hmm<ObservationInteger> buildPerturbedCoinProcess(double p) {
		double[][] a = new double[][] { { 1.0 - p, p }, { p, 1.0 - p } };
		List<OpdfInteger> pdfs = new ArrayList<OpdfInteger>(2);
		pdfs.add(new OpdfInteger(new double[] { 1.0, 0.0 }));
		pdfs.add(new OpdfInteger(new double[] { 0.0, 1.0 }));
		return (new HmmBuilder<ObservationInteger>(2)).withUniformPi().withA(a).withOpdfs(pdfs).build();
	}

	private static double analyticDivergenceRate(double p, double q) {
		return -0.5 * Math.log((p * q + (1.0 - p) * (1.0 - q)) /
				Math.sqrt((p * p + (1.0 - p) * (1.0 - p)) * (q * q + (1.0 - q) * (1.0 - q)))) / Math.log(2.0);
	}

	@Test
	public void testDivergenceRate() {
		Hmm<ObservationInteger> hmm1 = buildPerturbedCoinProcess(0.2);
		Hmm<ObservationInteger> hmm2 = buildPerturbedCoinProcess(0.5);
		assertEquals(0.0, HmmDivergence.divergenceRate(hmm1, hmm1), eps);
		double rFwd = HmmDivergence.divergenceRate(hmm1, hmm2);
		assertTrue(rFwd > 0.0);
		double rRev = HmmDivergence.divergenceRate(hmm2, hmm1);
		assertTrue(rRev > 0.0);
		assertEquals(rFwd, rRev, eps);
		assertEquals(analyticDivergenceRate(0.2, 0.5), rFwd, eps);
	}

}
