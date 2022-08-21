package edu.tum.cs.util;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.rng.UniformRandomProvider;

public class RngAdaptor implements RandomGenerator {

	private final UniformRandomProvider rng;

	public RngAdaptor(UniformRandomProvider rng) {
		this.rng = rng;
	}

	@Override
	public boolean nextBoolean() {
		return rng.nextBoolean();
	}

	@Override
	public void nextBytes(byte[] b) {
		rng.nextBytes(b);
	}

	@Override
	public double nextDouble() {
		return rng.nextDouble();
	}

	@Override
	public float nextFloat() {
		return rng.nextFloat();
	}

	@Override
	public double nextGaussian() {
		throw new UnsupportedOperationException();	// not implemented
	}

	@Override
	public int nextInt() {
		return rng.nextInt();
	}

	@Override
	public int nextInt(int n) {
		return rng.nextInt(n);
	}

	@Override
	public long nextLong() {
		return rng.nextLong();
	}

	@Override
	public void setSeed(int seed) {
		throw new UnsupportedOperationException();	// would need to recreate UniformRandomProvider from source
	}

	@Override
	public void setSeed(int[] seed) {
		throw new UnsupportedOperationException();	// would need to recreate UniformRandomProvider from source
	}

	@Override
	public void setSeed(long seed) {
		throw new UnsupportedOperationException();	// would need to recreate UniformRandomProvider from source
	}

}
