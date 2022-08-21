package edu.tum.cs.math.dist;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.Test;

public class DiscreteDistributionTest {

	@Test
	public void testJSD() {
		double[] pUnif = new double[10];
		Arrays.fill(pUnif, 1.0 / pUnif.length);
		assertEquals(0.0, DiscreteDistribution.distJS2(pUnif, pUnif), 0.0);
		assertEquals(1.0, DiscreteDistribution.distJS2(new double[] { 1.0, 0.0, 0.0, 0.0 },
				new double[] { 0.0, 1.0, 0.0, 0.0 }), 0.0);
	}

}
