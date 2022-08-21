package edu.tum.cs.time.hmm;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.HmmBuilder;
import be.ac.ulg.montefiore.run.jahmm.HmmLogSpace;
import be.ac.ulg.montefiore.run.jahmm.KMeansCalculator;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfFactory;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianFactory;
import be.ac.ulg.montefiore.run.jahmm.OpdfLogNormal;
import be.ac.ulg.montefiore.run.jahmm.OpdfLogNormalFactory;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryReader;
import be.ac.ulg.montefiore.run.jahmm.io.HmmBinaryWriter;

public class InteractionHmm extends Hmm<ObservationReal> {

	private static final long serialVersionUID = 4984913823217904424L;

	private static class CustomHmmBuilder extends HmmBuilder<ObservationReal> {
		public CustomHmmBuilder(int nbStates) {
			super(nbStates);
		}

		public CustomHmmBuilder withOpdfClustering(OpdfFactory<? extends Opdf<ObservationReal>> opdfFactory,
				List<? extends List<? extends ObservationReal>> sequences, boolean logTransform) {
			List<ObservationReal> observations = new ArrayList<ObservationReal>();
			for (List<? extends ObservationReal> s : sequences) {
				if (logTransform) {
					for (ObservationReal obs : s)
						observations.add(new ObservationReal(Math.log(obs.value)));
				} else
					observations.addAll(s);
			}

			KMeansCalculator<ObservationReal> kmc = new KMeansCalculator<ObservationReal>(nbStates, observations);
			opdfs = new ArrayList<Opdf<ObservationReal>>(nbStates);
			for (int i = 0; i < nbStates; i++) {
				Collection<ObservationReal> clusterObservations = kmc.cluster(i);

				Opdf<ObservationReal> opdf = opdfFactory.factor();
				if (!clusterObservations.isEmpty()) {
					if (logTransform) {
						List<ObservationReal> tfObs = new ArrayList<ObservationReal>(clusterObservations.size());
						for (ObservationReal obs : clusterObservations)
							tfObs.add(new ObservationReal(Math.exp(obs.value)));
						clusterObservations = tfObs;
					}
					opdf.fit(clusterObservations);
				} else
					System.err.println("empty cluster, using default parameters");
				opdfs.add(opdf);
			}
			return this;
		}
	}

	private final boolean logNormalEmissions;

	public InteractionHmm(Hmm<ObservationReal> hmm, boolean logNormalEmissions) {
		super(hmm.nbStates());
		for (int i = 0; i < nbStates(); i++) {
			setPi(i, hmm.getPi(i));
			setOpdf(i, hmm.getOpdf(i));
			for (int j = 0; j < nbStates(); j++)
				setAij(i, j, hmm.getAij(i, j));
		}
		this.logNormalEmissions = logNormalEmissions;
	}

	public boolean hasLogNormalEmissions() {
		return logNormalEmissions;
	}

	public double computeLogLikelihood(List<TimeSequence> sequences) {
		double logLikelihood = 0.0;
		HmmLogSpace<ObservationReal> hmmLog = new HmmLogSpace<ObservationReal>(this);
		for (TimeSequence sequence : sequences)
			logLikelihood += BaumWelchLearnerIncomplete.lnProbability(hmmLog, sequence.data, sequence.isIncomplete);
		return logLikelihood;
	}

	public double computeBic(double logLikelihood, int numObs) {
		int numStates = nbStates();
		int numParam = (numStates * numStates) + (2 * numStates) - 1;	// same for normal and log-normal emissions
		return (-2.0 * logLikelihood) + (numParam * Math.log(numObs));
	}

	public static void write(File f, Hmm<ObservationReal> hmm) throws IOException {
		OutputStream os = new FileOutputStream(f);
		try {
			HmmBinaryWriter.write(os, hmm);
		} finally {
			os.close();
		}
	}

	@SuppressWarnings("unchecked")
	public static InteractionHmm read(File f) throws IOException {
		InputStream is = new FileInputStream(f);
		try {
			Hmm<ObservationReal> hmm = (Hmm<ObservationReal>) HmmBinaryReader.read(is);
			boolean logNormalEmissions = hmm.getOpdf(0) instanceof OpdfLogNormal;
			return new InteractionHmm(hmm, logNormalEmissions);
		} finally {
			is.close();
		}
	}

	public static InteractionHmm fit(int numStates, boolean logNormalEmissions, int maxIter, double minLikelihoodDelta,
			List<TimeSequence> sequences) {
		return fit(numStates, logNormalEmissions, maxIter, minLikelihoodDelta, null, sequences, 0, null);
	}

	public static InteractionHmm fit(int numStates, boolean logNormalEmissions, int maxIter, double minLikelihoodDelta,
			List<TimeSequence> sequences, int reportInterval, File intermediateStateFile) {
		return fit(numStates, logNormalEmissions, maxIter, minLikelihoodDelta, null, sequences, reportInterval,
				intermediateStateFile);
	}

	public static InteractionHmm resume(InteractionHmm curHmm, int maxIter, double minLikelihoodDelta,
			List<TimeSequence> sequences, int reportInterval, File intermediateStateFile) {
		return fit(curHmm.nbStates(), curHmm.logNormalEmissions, maxIter, minLikelihoodDelta, curHmm, sequences,
				reportInterval, intermediateStateFile);
	}

	private static InteractionHmm fit(int numStates, boolean logNormalEmissions, int maxIter, double minLikelihoodDelta,
			Hmm<ObservationReal> curHmm, List<TimeSequence> sequences, int reportInterval, File intermediateStateFile) {
		List<List<ObservationReal>> sequenceData = TimeSequence.getObservations(sequences);
		boolean[] isIncomplete = TimeSequence.getCompleteness(sequences);

		Hmm<ObservationReal> initialHmm;
		if (curHmm == null) {
			// use k-means to generate initial output parameters
			OpdfFactory<? extends Opdf<ObservationReal>> factory;
			if (logNormalEmissions)
				factory = new OpdfLogNormalFactory();
			else
				factory = new OpdfGaussianFactory();
			initialHmm = (new CustomHmmBuilder(numStates))
					.withOpdfClustering(factory, sequenceData, logNormalEmissions)
					.withUniformPi().withUniformA().build();
		} else {
			initialHmm = curHmm;
		}

		HmmLogSpace<ObservationReal> hmm = new HmmLogSpace<ObservationReal>(initialHmm);
		HmmLogSpace<ObservationReal> prevHmm;
		try {
			prevHmm = hmm.clone();
		} catch (CloneNotSupportedException ex) {
			throw new RuntimeException(ex);
		}

		// ML estimation of parameters
		BaumWelchLearnerIncomplete bwLearner = new BaumWelchLearnerIncomplete();
		double prevLl = 0.0, ll = Double.NEGATIVE_INFINITY;
		int i;
		for (i = 0; i < maxIter; i++) {
			hmm = bwLearner.iterate(prevHmm, sequenceData, isIncomplete);

			prevLl = ll;
			ll = 0.0;
			for (TimeSequence sequence : sequences)
				ll += BaumWelchLearnerIncomplete.lnProbability(hmm, sequence.data, sequence.isIncomplete);

			if (ll < prevLl)
				break;
			prevHmm = hmm;
			if ((ll - prevLl) < minLikelihoodDelta)
				break;

			if ((reportInterval > 0) && (((i + 1) % reportInterval) == 0)) {
				System.err.println((i + 1) + " iterations, LL = " + ll);
				try {
					InteractionHmm.write(intermediateStateFile, hmm.toHmm());
				} catch (IOException ex) {
					System.err.println("warning: could not save intermediate HMM state");
					ex.printStackTrace(System.err);
				}
			}
		}

		System.err.println("finished training HMM with " + numStates + " states after " + (i - 1) +
				" iterations, LL = " + ll + ", delta = " + (ll - prevLl));
		return new InteractionHmm(prevHmm.toHmm(), logNormalEmissions);
	}

}
