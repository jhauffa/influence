package edu.tum.cs.time.hmm;

import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardCalculatorLogSpace;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.HmmLogSpace;
import be.ac.ulg.montefiore.run.jahmm.LogSpace;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ForwardBackwardCalculator.Computation;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchLearnerLogSpace;

public class BaumWelchLearnerIncomplete extends BaumWelchLearnerLogSpace {

	private static class ForwardBackwardCalculatorIncomplete extends ForwardBackwardCalculatorLogSpace {
		private final boolean isIncomplete;

		public <O extends Observation> ForwardBackwardCalculatorIncomplete(List<? extends O> oseq, boolean isIncomplete,
				HmmLogSpace<O> hmm, EnumSet<Computation> flags) {
			super(oseq, hmm, EnumSet.noneOf(Computation.class));
			this.isIncomplete = isIncomplete;
			compute(oseq, hmm, flags);
		}

		@Override
		protected <O extends Observation> void computeAlphaInit(Hmm<? super O> hmm, O o, int i) {
			double pi = isIncomplete ? -Math.log(hmm.nbStates()) : hmm.getPi(i);
			alpha[0][i] = LogSpace.product(pi, hmm.getOpdf(i).logProbability(o));
		}
	}

	public <O extends Observation> HmmLogSpace<O> iterate(HmmLogSpace<O> hmm,
			List<? extends List<? extends O>> sequences, boolean[] isIncomplete) {
		HmmLogSpace<O> nhmm;
		try {
			nhmm = hmm.clone();
		} catch(CloneNotSupportedException e) {
			throw new InternalError();
		}

		double allGamma[][][] = new double[sequences.size()][][];

		/* update transition probabilities Aij */
		double aijDen[] = new double[hmm.nbStates()];
		double aijNum[][] = new double[hmm.nbStates()][hmm.nbStates()];

		Arrays.fill(aijDen, LogSpace.ZERO);
		for (int i = 0; i < hmm.nbStates(); i++)
			Arrays.fill(aijNum[i], LogSpace.ZERO);

		int g = 0;
		for (List<? extends O> obsSeq : sequences) {
			ForwardBackwardCalculatorLogSpace fbc = new ForwardBackwardCalculatorIncomplete(obsSeq, isIncomplete[g],
					nhmm, EnumSet.of(Computation.ALPHA, Computation.BETA));

			double xi[][][] = estimateXi(obsSeq, fbc, nhmm);
			double gamma[][] = allGamma[g++] = estimateGamma(xi, fbc);

			double[] os = new double[obsSeq.size()];
			for (int i = 0; i < hmm.nbStates(); i++) {
				os[0] = aijDen[i];
				for (int t = 0; t < os.length - 1; t++)
					os[t + 1] = gamma[t][i];
				aijDen[i] = LogSpace.INST.sum(os);

				for (int j = 0; j < hmm.nbStates(); j++) {
					os[0] = aijNum[i][j];
					for (int t = 0; t < os.length - 1; t++)
						os[t + 1] = xi[t][i][j];
					aijNum[i][j] = LogSpace.INST.sum(os);
				}
			}
		}

		for (int i = 0; i < hmm.nbStates(); i++)
			for (int j = 0; j < hmm.nbStates(); j++)
				nhmm.setAij(i, j, LogSpace.quotient(aijNum[i][j], aijDen[i]));

		/* update initial probabilities Pi */
		int numComplete = 0;
		for (boolean b : isIncomplete)
			if (!b)
				numComplete++;
		double piDen = LogSpace.INST.log(numComplete);

		double[] s = new double[numComplete];
		for (int i = 0; i < hmm.nbStates(); i++) {
			int idx = 0;
			for (int o = 0; o < sequences.size(); o++)
				if (!isIncomplete[o])
					s[idx++] = allGamma[o][0][i];
			nhmm.setPi(i, LogSpace.quotient(LogSpace.INST.sum(s), piDen));
		}

		/* update output PDFs */
		List<O> observations = Observation.flat(sequences);
		double[] weights = new double[observations.size()];
		for (int i = 0; i < hmm.nbStates(); i++) {
			int obsIdx = 0;
			for (double[][] gamma : allGamma) {
				for (double[] gammaT : gamma)
					weights[obsIdx++] = gammaT[i];
			}
			double sum = LogSpace.INST.sum(weights);

			for (int j = 0; j < obsIdx; j++)
				weights[j] = LogSpace.INST.exp(LogSpace.quotient(weights[j], sum));

			nhmm.getOpdf(i).fit(observations, weights);
		}

		return nhmm;
	}

	public static <O extends Observation> double lnProbability(HmmLogSpace<O> hmm, List<? extends O> oseq,
			boolean isIncomplete) {
		ForwardBackwardCalculatorIncomplete fbc = new ForwardBackwardCalculatorIncomplete(oseq, isIncomplete,
				hmm, EnumSet.of(Computation.ALPHA, Computation.PROBABILITY));
		return fbc.lnProbability();
	}

}
