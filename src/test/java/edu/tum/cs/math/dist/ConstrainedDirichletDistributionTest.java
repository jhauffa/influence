package edu.tum.cs.math.dist;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.junit.Test;

public class ConstrainedDirichletDistributionTest {

	private static final double eps = 1E-10;

	private static void checkDistribution(String id, double[] x, double[] c, double r) {
		double sum = 0.0;
		double sumAbsDiff = 0.0;
		for (int i = 0; i < x.length; i++) {
			assertTrue(id + ": x[" + i + "] = " + x[i], x[i] > -eps);
			sum += x[i];
			sumAbsDiff += Math.abs(x[i] - c[i]);
		}
		assertEquals(id, 1.0, sum, eps);
		assertEquals(id, r, sumAbsDiff, eps);
	}

	private static final long seed = 13371338L;
	private static final int numTestDistrSamples = 100;
	private static final int[] dimensions = { 3, 10, 200 };
	private static enum CenterVariants { BARYCENTER, OFF_CENTER, CORNER_1, CORNER_N };
	private static enum RadiusVariants { INTERSECT_INNER, MAXIMUM_INNER, INTERSECT_OUTER, MAXIMUM_OUTER };
	private static enum AlphaVariants { UNIFORM, SYMMETRIC_SPARSE, SYMMETRIC_DENSE, ASYMMETRIC, ASYMMETRIC_N };

	@Test
	public void testSampling() {
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
		for (int n : dimensions) {
			for (CenterVariants cv : CenterVariants.values()) {
				double[] c = new double[n];
				double cMin = 0.0;
				switch (cv) {
				case BARYCENTER:
					Arrays.fill(c, 1.0 / n);
					cMin = 1.0 / n;
					break;
				case OFF_CENTER:
					for (int i = 0; i < (n / 2); i++) {
						c[2 * i] = (1.0 / n) + (1.0 / (2 * n));
						c[2 * i + 1] = (1.0 / n) - (1.0 / (2 * n));
					}
					if ((n % 2) == 1)
						c[n - 1] = 1.0 / n;
					cMin = (1.0 / n) - (1.0 / (2 * n));
					break;
				case CORNER_1:
					c[0] = 1.0;
					break;
				case CORNER_N:
					c[n - 1] = 1.0;
					break;
				}

				for (RadiusVariants rv : RadiusVariants.values()) {
					double r = 0.0;
					switch (rv) {
					case INTERSECT_INNER:
						r = (1.0 - (1.0 / n)) / 2.0;
						break;
					case MAXIMUM_INNER:
						r = 1.0 - (1.0 / n);
						break;
					case INTERSECT_OUTER:
						r = ((1.0 - (1.0 / n)) + (2.0 - (2.0 * cMin))) / 2.0;
						break;
					case MAXIMUM_OUTER:
						r = 2.0 - (2.0 * cMin);
						break;
					}

					for (AlphaVariants av : AlphaVariants.values()) {
						double[] alpha = new double[n];
						switch (av) {
						case UNIFORM:
							Arrays.fill(alpha, 1.0);
							break;
						case SYMMETRIC_SPARSE:
							Arrays.fill(alpha, 0.1);
							break;
						case SYMMETRIC_DENSE:
							Arrays.fill(alpha, 10.0);
							break;
						case ASYMMETRIC:
							for (int i = 0; i < n; i++)
								alpha[i] = (i + 1) * 0.3;
							break;
						case ASYMMETRIC_N:
							for (int i = 0; i < (n - 1); i++)
								alpha[i] = (i + 1) * 0.3;
							alpha[n - 1] = 1.0;
							break;
						}

						String id = n + "-" + cv + "-" + rv + "-" + av;
						ConstrainedDirichletDistribution dist = new ConstrainedDirichletDistribution(alpha, c, r, rand);
						checkDistribution("init-" + id, dist.x, c, r);	// initial state
						checkDistribution("burnin-" + id, dist.sample(), c, r);	// state after burn-in
						for (int i = 0; i < numTestDistrSamples; i++)
							checkDistribution("sample-" + i + "-" + id, dist.sample(), c, r);	// state after lag
					}
				}
			}
		}
	}


	private static final double[] typicalAlpha = {
		0.0008104277205577801, 0.002399302458610846, 0.0012983374920475725, 0.0007206217250330912, 0.001264593802616527,
		0.0031008761502928964, 0.012847727941347491, 0.0006495064959995071, 0.0017244359408497654,
		0.0036243036912556093, 0.0006177435011521368, 0.013735195103145043, 0.013355890190236315, 0.0021096782841007926,
		0.002161777300417673, 0.0028235076402647855, 0.002478394515613635, 0.0010558349851108322, 0.001493902247725142,
		0.002056064681495052, 0.010881168896780105, 0.0016877044586594355, 0.02739233025614923, 0.0024983597554522215,
		0.008524513777939748, 0.0024084824394314712, 0.0009922486718967542, 0.0008417966128777247,
		0.0017297158979513662, 0.00048096516588401483, 0.0028378938293745653, 0.004461835495589483,
		0.0024281218673292437, 0.001481351772888192, 0.0003257931478476927, 0.0011928096296026259, 0.009828055480426583,
		0.0005177418691815221, 0.0014493442597134918, 0.000575590570353794, 0.0009370714498753967, 0.000907783825384802,
		0.0012880682660478004, 0.0009740697052998162, 0.026005186767064956, 0.0012195443717341451, 0.01586029303578373,
		0.07003872710703464, 0.002964575757378324, 0.00354371066768484, 0.0004388262958637487, 0.0033572869639097867,
		0.0047171482335151, 0.0010636175528425898, 0.02045524255581931, 0.0024117991986736947, 0.0025607057419328684,
		0.046692034844982314, 0.002166247427981139, 0.11299395218054076, 0.010783048686776726, 0.0011344878822794214,
		0.0020298184493623282, 0.0022433996877128685, 0.019853953979434555, 0.0008524770775692356,
		0.0009823981402538545, 0.0018053055343887634, 0.001472437390884757, 0.003381766614127223, 0.0015503002508639458,
		0.0006020099494867585, 0.005246580591604407, 0.0014067173582097464, 0.05242434898824254, 0.0006283573731112969,
		0.0017874715394979247, 0.027111493559285025, 0.001451459902491616, 0.003061187903629909, 0.0015316208440770973,
		0.0009052870205427502, 0.009293350098694567, 0.0010633803735597119, 0.00910080646494858, 0.016665057757745953,
		0.0004914565987789147, 0.0010293319861682437, 0.035964241094363324, 0.008348928441112522, 0.0034016112267820703,
		0.0008578194308639995, 0.00047845125790671537, 0.0018711870618725243, 0.0007260355429333649,
		0.00044673152310053784, 0.0004493478106586614, 0.0016116131985771888, 0.0008576155330873244,
		0.0005282333004744064, 0.006022181778766256, 0.0004888849417407031, 0.003688756311334507, 0.0026914192154777854,
		0.002500026712324937, 0.0007152791560293635, 0.0008287702204444501, 0.002779183368986716, 0.0203824359641677,
		0.0018830146291848205, 0.005005481006910321, 0.0017046999151801739, 0.0067082822588582156,
		0.0015956495190410627, 0.21493351572724176, 0.0009184702699272461, 0.0016111155650691301, 0.004028296557061569,
		0.0017040089887328978, 0.005674023267589124, 0.0011185742911991943, 0.0012303884377255967,
		0.0031370267991522585, 0.05817412422220666, 0.0012360448726458874, 0.02887029498164493, 0.0033143738721637464,
		0.008273726222718372, 0.0013893084160678173, 0.0006812553082574134, 0.0001470148382277636,
		0.0012595154674669476, 0.00262242840578705, 0.0019328793326867707, 0.0007575106776943911, 0.0025893836909075363,
		0.0004336211910761012, 0.0024451135534129488, 0.0011451947736801236, 0.0007705477690282891,
		0.003680738986112686, 0.0020415974192589705, 0.017961998877061073, 0.0020199418118046607,
		0.00044932914354601047, 0.0049369488627468185, 0.0013042474424475173, 0.0005941889596499517,
		0.0065110318319204314, 0.006349303972798182
	};

	@Test
	public void testTypicalLDAAlpha() {
		int n = typicalAlpha.length;
		double[] c = new double[n];
		Arrays.fill(c, 1.0 / n);
		double r = 1.0 - (1.0 / n);

		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
		ConstrainedDirichletDistribution dist = new ConstrainedDirichletDistribution(typicalAlpha, c, r, rand);
		checkDistribution("init", dist.x, c, r);
		checkDistribution("burnin", dist.sample(), c, r);
		for (int i = 0; i < numTestDistrSamples; i++)
			checkDistribution("sample-" + i, dist.sample(), c, r);
	}


	private static final int[] numDim = { 3, 5, 20 };
	private static final int numSamplesPerDim = 1000;
	private static final double[] alphaValues = { 0.01, 1.0, 10.0 };
	private static final double maxRelDiff = 0.1;

	@Test
	public void testDistributionProperties() {
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
		for (int d : numDim) {
			for (double a : alphaValues) {
				double[] alpha = new double[d];
				Arrays.fill(alpha, a);
				double[] c = new double[d];
				Arrays.fill(c, 1.0 / d);
				double r;
				if (a == 1.0)
					r = 2.0 - (2.0 / d);	// maximum outer
				else
					r = 1.0 - (1.0 / d);	// maximum inner
				ConstrainedDirichletDistribution dist = new ConstrainedDirichletDistribution(alpha, c, r, rand);

				int[] histo = new int[d];
				for (int i = 0; i < (d * numSamplesPerDim); i++) {
					double[] x = dist.sample();

					double ext = (a >= 1.0) ? 0.0 : 1.0;
					int[] extIdx = new int[d];
					int lastExt = 0;
					for (int j = 0; j < x.length; j++) {
						if (((a >= 1.0) && (x[j] > ext)) ||
							((a < 1.0) && (x[j] < ext))) {
							ext = x[j];
							lastExt = 0;
							extIdx[lastExt++] = j;
						} else if (x[j] == ext)
							extIdx[lastExt++] = j;
					}
					if (a == 1.0)
						assertTrue(Math.abs(1.0 - ext) <= eps);
					int idx = (int) (rand.nextDouble() * lastExt);
					histo[extIdx[idx]]++;
				}

				double max = 0.0;
				for (int n : histo) {
					assertTrue("alpha = " + a + ", " + d + " dim.", n > 0);
					double diff = Math.abs(1.0 - ((double) n / numSamplesPerDim));
					if (diff > max)
						max = diff;
				}
				assertTrue("alpha = " + a + ", " + d + " dim.: " + max + " >= " + maxRelDiff, max < maxRelDiff);
			}
		}
	}

}
