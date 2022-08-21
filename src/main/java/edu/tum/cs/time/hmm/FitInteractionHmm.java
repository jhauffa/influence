package edu.tum.cs.time.hmm;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import be.ac.ulg.montefiore.run.jahmm.LogSpace;
import edu.tum.cs.math.FastMath;
import edu.tum.cs.math.dist.DiscreteDistribution;
import edu.tum.cs.util.ExperimentConfiguration;

public class FitInteractionHmm {

	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(FitInteractionHmm.class);

	private static final int numCvFolds = 10;
	private static final int minNumStates = 1;
	private static final int maxNumStates = 21;
	private static final int maxIterDevel = 100;
	private static final double minDeltaDevel = 0.1;
	private static final int maxIterFit = 20000;
	private static final double minDeltaFit = 1E-6;
	private static String outputDir = ".";

	private static class CvTask implements Callable<CvTask> {
		private final List<List<TimeSequence>> folds;
		public final int foldIdx;
		public final int numStates;
		private final boolean logNormalEmissions;
		public double bic;

		public CvTask(List<List<TimeSequence>> folds, int foldIdx, int numStates, boolean logNormalEmissions) {
			this.folds = folds;
			this.foldIdx = foldIdx;
			this.numStates = numStates;
			this.logNormalEmissions = logNormalEmissions;
		}

		public CvTask call() throws Exception {
			List<TimeSequence> foldEval = folds.get(foldIdx);
			List<TimeSequence> foldFit = new ArrayList<TimeSequence>();
			for (int i = 0; i < numCvFolds; i++) {
				if (i != foldIdx)
					foldFit.addAll(folds.get(i));
			}

			try {
				InteractionHmm hmm = InteractionHmm.fit(numStates, logNormalEmissions, maxIterDevel, minDeltaDevel,
						foldFit);
				double ll = hmm.computeLogLikelihood(foldEval);
				int numObs = TimeSequence.countObservations(foldEval);
				bic = hmm.computeBic(ll, numObs);
			} catch (IllegalArgumentException ex) {
				System.err.println("warning: CV HMM converged to degenerate solution, ignoring");
				bic = Double.POSITIVE_INFINITY;
			}
			return this;
		}
	}

	private static int findOptimalNumStates(List<TimeSequence> sequences, boolean logNormalEmissions) {
		List<List<TimeSequence>> folds = new ArrayList<List<TimeSequence>>(numCvFolds);
		int foldSize = sequences.size() / numCvFolds;
		if (foldSize == 0)
			throw new RuntimeException("insufficient training data for CV");
		int startIdx = 0, endIdx = foldSize;
		for (int i = 0; i < numCvFolds; i++) {
			if (i < (numCvFolds - 1))
				folds.add(sequences.subList(startIdx, endIdx));
			else
				folds.add(sequences.subList(startIdx, sequences.size()));
			startIdx = endIdx;
			endIdx += foldSize;
		}

		double[][] bic = new double[maxNumStates - minNumStates][numCvFolds + 1];
		int numThreads = cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_THREADS,
				Runtime.getRuntime().availableProcessors());
		ExecutorService pool = Executors.newFixedThreadPool(numThreads);
		try {
			CompletionService<CvTask> tasks = new ExecutorCompletionService<CvTask>(pool);
			for (int n = minNumStates; n < maxNumStates; n++)
				for (int i = 0; i < numCvFolds; i++)
					tasks.submit(new CvTask(folds, i, n, logNormalEmissions));

			for (int i = 0; i < ((maxNumStates - minNumStates) * numCvFolds); i++) {
				CvTask result;
				try {
					Future<CvTask> finishedTask = tasks.take();
					result = finishedTask.get();
				} catch (Exception ex) {
					throw new RuntimeException("error in CV task", ex);
				}
				bic[result.numStates - minNumStates][result.foldIdx] = result.bic;
			}
		} finally {
			pool.shutdownNow();
		}

		double lowestBic = Double.POSITIVE_INFINITY;
		int bestNumStates = 0;
		for (int n = 0; n < (maxNumStates - minNumStates); n++) {
			int effectiveNumCvFolds = 0;
			for (int i = 0; i < numCvFolds; i++) {
				if (!Double.isInfinite(bic[n][i])) {
					bic[n][numCvFolds] += bic[n][i];
					effectiveNumCvFolds++;
				}
			}

			bic[n][numCvFolds] /= effectiveNumCvFolds;
			if (bic[n][numCvFolds] < lowestBic) {
				lowestBic = bic[n][numCvFolds];
				bestNumStates = minNumStates + n;
			}
		}

		for (int i = 0; i < bic.length; i++) {
			System.out.print((minNumStates + i));
			for (int j = 0; j < bic[i].length; j++)
				System.out.print(";" + bic[i][j]);
			System.out.println();
		}
		return bestNumStates;
	}

	private static void computePerplexity(InteractionHmm hmm, List<TimeSequence> sequencesFit,
			List<TimeSequence> sequencesEval) {
		// compare perplexity on held-out and training data to detect overfitting
		double ll = hmm.computeLogLikelihood(sequencesFit);
		int numObs = TimeSequence.countObservations(sequencesFit);
		double bic = hmm.computeBic(ll, numObs);
		double pp = DiscreteDistribution.perplexity(ll, numObs);
		System.out.println("training set: LL = " + ll + ", BIC = " + bic + ", perplexity = " + pp);
		ll = hmm.computeLogLikelihood(sequencesEval);
		numObs = TimeSequence.countObservations(sequencesEval);
		bic = hmm.computeBic(ll, numObs);
		pp = DiscreteDistribution.perplexity(ll, numObs);
		System.out.println("test set: LL = " + ll + ", BIC = " + bic + ", perplexity = " + pp);
	}

	private static void startHmmTraining(String[] args) throws Exception {
		if (args.length > 2)
			outputDir = args[2];
		boolean logNormalEmissions = false;
		if (args.length > 3)
			logNormalEmissions = Boolean.parseBoolean(args[3]);
		double subsamplingFactor = 1.0;
		if (args.length > 4)
			subsamplingFactor = Double.parseDouble(args[4]) / 100.0;

		// read sequences and split randomly into 10% development, 80% training, 10% test
		List<TimeSequence> sequences = TimeSequence.readSequences(new File(args[1]));
		Collections.shuffle(sequences);
		if (subsamplingFactor < 1.0)
			sequences = sequences.subList(0, (int) (sequences.size() * subsamplingFactor));
		int numDevelSequences = sequences.size() / 10;
		List<TimeSequence> sequencesDevel = sequences.subList(0, numDevelSequences);
		List<TimeSequence> sequencesFit = sequences.subList(numDevelSequences, sequences.size() - numDevelSequences);
		List<TimeSequence> sequencesEval = sequences.subList(sequences.size() - numDevelSequences, sequences.size());
		System.out.println("#seq = " + sequences.size() + ", #devel = " + sequencesDevel.size() +
				", #fit = " + sequencesFit.size() + ", #eval = " + sequencesEval.size());

		// save the individual partitions to make the experiment reproducible
		TimeSequence.writeSequences(new File(outputDir, "chains-devel.txt.gz"), sequencesDevel);
		TimeSequence.writeSequences(new File(outputDir, "chains-fit.txt.gz"), sequencesFit);
		TimeSequence.writeSequences(new File(outputDir, "chains-eval.txt.gz"), sequencesEval);

		// use development set to determine optimal number of states according to BIC and fit HMM to training data
		int numStates = findOptimalNumStates(sequencesDevel, logNormalEmissions);
		System.out.println("fitting HMM with " + numStates + " states and " + (logNormalEmissions ? "log-" : "") +
				"normal emissions to training set...");
		File paramFile = new File(outputDir, "hmm-param.bin");
		InteractionHmm hmm = InteractionHmm.fit(numStates, logNormalEmissions, maxIterFit, minDeltaFit, sequencesFit,
				100, paramFile);
		InteractionHmm.write(paramFile, hmm);
		computePerplexity(hmm, sequencesFit, sequencesEval);
	}

	private static void resumeHmmTraining(String[] args) throws Exception {
		InteractionHmm hmm = InteractionHmm.read(new File(args[1]));
		outputDir = args[2];
		List<TimeSequence> sequencesFit = TimeSequence.readSequences(new File(outputDir, "chains-fit.txt.gz"));
		List<TimeSequence> sequencesEval = TimeSequence.readSequences(new File(outputDir, "chains-eval.txt.gz"));

		int maxIter = maxIterFit;
		if (args.length > 3)
			maxIter = Integer.parseInt(args[3]);
		if (maxIter > 0) {
			File paramFile = new File(outputDir, "hmm-param.bin");
			paramFile.renameTo(new File(outputDir, "hmm-param.old"));

			hmm = InteractionHmm.resume(hmm, maxIter, minDeltaFit, sequencesFit, 100, paramFile);
			InteractionHmm.write(paramFile, hmm);
		}
		computePerplexity(hmm, sequencesFit, sequencesEval);
	}

	private static class LogSpaceJafama extends LogSpace {
		@Override
		public double expImpl(double x) {
			return FastMath.exp(x);
		}

		@Override
		public double logImpl(double x) {
			return FastMath.log(x);
		}
	}

	public static void main(String[] args) throws Exception {
		LogSpace.INST = new LogSpaceJafama();

		if (args.length < 1) {
			System.err.println("usage: " + FitInteractionHmm.class.getSimpleName() + "fit/resume");
			System.err.println("arguments for 'fit': chains [out] [log-normal?] [%subsampling]");
			System.err.println("arguments for 'resume': hmm out [max. iter.]");
			return;
		}
		if (args[0].equals("fit"))
			startHmmTraining(args);
		else if (args[0].equals("resume"))
			resumeHmmTraining(args);
		else
			System.err.println("unknown command '" + args[0] + "'");
	}

}
