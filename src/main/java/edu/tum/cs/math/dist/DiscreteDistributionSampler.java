package edu.tum.cs.math.dist;

import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import org.apache.commons.rng.RestorableUniformRandomProvider;
import org.apache.commons.rng.core.RandomProviderDefaultState;
import org.apache.commons.rng.simple.RandomSource;

public class DiscreteDistributionSampler implements Serializable {

	private static final long serialVersionUID = -6779194641750140596L;

	public transient RestorableUniformRandomProvider prng;

	public DiscreteDistributionSampler() {
		prng = RandomSource.create(RandomSource.XOR_SHIFT_1024_S);
	}

	public DiscreteDistributionSampler(DiscreteDistributionSampler other) {
		this(other.prng.nextLong());
	}

	public DiscreteDistributionSampler(long seed) {
		prng = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
	}

	public int sample1DUniform(int numEvents) {
		if (numEvents == 1)
			return 0;
		return prng.nextInt(numEvents);
	}

	public int sample1D(double[] pCumul) {
		return sample1D(prng.nextDouble(), pCumul, pCumul.length);
	}

	public int sample1D(double[] pCumul, int n) {
		return sample1D(prng.nextDouble(), pCumul, n);
	}

	public static int sample1D(double p, double[] pCumul) {
		return sample1D(p, pCumul, pCumul.length);
	}

	public static int sample1D(double p, double[] pCumul, int n) {
		p *= pCumul[n - 1];
		int left = 0;
		int right = n - 1;
		int v = right >>> 1;
		// find v so that p >= pCumul[v - 1] && p < pCumul[v]
		while (left < right) {
			if (p >= pCumul[v])
				left = v + 1;
			else
				right = v;
			v = (left + right) >>> 1;
		}
		return v;
	}

	public int[] sample1DUniformWithoutReplacement(int numEvents, int numSamples) {
		numSamples = Math.min(numEvents, numSamples);
		int[] samples = new int[numSamples];
		int n = 0, t = 0;
		while (n < numSamples) {
			double u = prng.nextDouble();
			if ((u * (numEvents - t)) < (numSamples - n))
				samples[n++] = t;
			t++;
		}
		return samples;
	}

	private void writeObject(ObjectOutputStream out) throws IOException {
		out.defaultWriteObject();
		out.writeObject(((RandomProviderDefaultState) prng.saveState()).getState());
	}

	private void readObject(java.io.ObjectInputStream in) throws IOException, ClassNotFoundException {
		in.defaultReadObject();
		prng = RandomSource.create(RandomSource.XOR_SHIFT_1024_S);
		prng.restoreState(new RandomProviderDefaultState((byte[]) in.readObject()));
	}

}
