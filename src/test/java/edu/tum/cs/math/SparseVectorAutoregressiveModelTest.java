package edu.tum.cs.math;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;

public class SparseVectorAutoregressiveModelTest {

	private static final double[] coeffHypo = {
		1.78897214299606, 0, 0, 0.385917198573085, 0, 0, 0, 0, 0, -0.530822612617507, 0, 0, 0, 0, -0.056872212678753, 0, 0, 0, 0, 0
	};

	private static final double[] p = {
		0, 0.0122120125112268, 0.74273960062942, 1.05982802389715e-05, 0.443676994645771, 0.839791916728969, 0.618659634879648, 0.55259507982187, 0.0730375641819363, 0.00033718099823199, 0.777980892203379, 0.846668717857303, 0.951740930282193, 0.593280787157083, 0.0537165121262881, 0.392424191551892, 0.534478916131054, 0.397493670462872, 0.709798561653612, 0.487915820409302
	};

	private static final double eps = 1E-2;

	@Test
	public void testHypothesisTest() {
		RealMatrix x = new Array2DRowRealMatrix(SqrtLassoTest.xData, true);
		RealVector y = new ArrayRealVector(SqrtLassoTest.yData, true);
		SparseVectorAutoregressiveModel.DebiasedLassoHypothesisTest test =
				new SparseVectorAutoregressiveModel.DebiasedLassoHypothesisTest(x.getData(), false);
		RealVector[] res = test.performTest(y);
		assertEquals(3, res.length);
		assertEquals(coeffHypo.length, res[0].getDimension());
		for (int i = 0; i < coeffHypo.length; i++)
			assertEquals(coeffHypo[i], res[0].getEntry(i), eps);
		assertEquals(p.length, res[2].getDimension());
		for (int i = 0; i < p.length; i++)
			assertEquals(p[i], res[2].getEntry(i), eps);
	}

	private static final long seed = 4546546L;
	private static final int numDim = 60;
	private static final int numSamples = 20;
	private static final double maxError = 0.0025;

	@Test
	public void testParameterFitting() {
		RandomGenerator rand = new MersenneTwister(seed);
		double pNonZero = Math.sqrt(numSamples - 1) / (2.0 * numDim * Math.log(numDim));

		// generate samples from VAR given a random starting point
		SparseVectorAutoregressiveModel trueModel = SparseVectorAutoregressiveModel.random(rand, numDim, pNonZero, 0.2);
		RealMatrix samples = trueModel.generate(rand, numSamples * 2, 150, 0.25);

		// attempt to recover VAR coefficients
		SparseVectorAutoregressiveModel estModel1 = SparseVectorAutoregressiveModel.fit(
				samples.getSubMatrix(0, numDim - 1, 0, numSamples - 1).getData(), true);
		RealMatrix delta = estModel1.getCoefficients().subtract(trueModel.getCoefficients());
		double mse1 = delta.multiply(delta.transpose()).getTrace() / (numDim * numDim);
		assertTrue(Double.toString(mse1), mse1 < maxError);
		SparseVectorAutoregressiveModel estModel2 = SparseVectorAutoregressiveModel.fit(samples.getData(), true);
		delta = estModel2.getCoefficients().subtract(trueModel.getCoefficients());
		double mse2 = delta.multiply(delta.transpose()).getTrace() / (numDim * numDim);
		assertTrue(mse1 + ";" + mse2, mse2 < mse1);
	}

}
