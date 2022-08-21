package edu.tum.cs.math.dist;

import static org.junit.Assert.assertEquals;

import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

public class DiscreteDistributionSamplerTest {

	@Test
	public void testSamplingWithoutReplacement() {
		// when sampling all values from a finite range, the result should be a permutation
		DiscreteDistributionSampler sampler = new DiscreteDistributionSampler();
		int[] perm = sampler.sample1DUniformWithoutReplacement(100, 100);
		Set<Integer> unique = new HashSet<Integer>(perm.length);
		for (int e : perm)
			unique.add(e);
		assertEquals(perm.length, unique.size());
	}

}
