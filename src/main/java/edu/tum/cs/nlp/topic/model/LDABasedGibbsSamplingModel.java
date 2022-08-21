package edu.tum.cs.nlp.topic.model;

import java.io.File;
import java.util.Arrays;
import java.util.Comparator;

import edu.tum.cs.math.dist.DiscreteDistributionSampler;
import edu.tum.cs.nlp.corpus.Corpus;
import edu.tum.cs.nlp.corpus.ProcessedDocument;
import edu.tum.cs.util.arrays.DenseDouble2DArray;
import edu.tum.cs.util.arrays.DenseInt2DArray;
import edu.tum.cs.util.arrays.Double2DArray;
import edu.tum.cs.util.arrays.FixedDouble2DArray;
import edu.tum.cs.util.arrays.Int2DArray;

public abstract class LDABasedGibbsSamplingModel<T extends ProcessedDocument> implements TopicModel<T> {

	private static final long serialVersionUID = 4643635705853597312L;

	public static enum SamplingType { REGULAR, FASTLDA };

	public abstract class SamplingState implements TopicModel.SamplingState {
		private static final long serialVersionUID = -468391495909452990L;

		protected final boolean isQuery;
		protected final SamplingType type;

		protected Int2DArray numWordTopic;
		protected int[] sumNumTopic;

		// FastLDA
		protected Int2DArray sortedWordTopic;
		protected int[] sortedTopic;
		protected transient Int2DArray invSortedWordTopic;
		protected transient int[] invSortedTopic;
		protected double[] bNormSquared;

		protected SamplingState(boolean isQuery, SamplingType type) {
			this.isQuery = isQuery;
			this.type = type;

			if (!isQuery) {
				numWordTopic = new DenseInt2DArray(numWords, numTopics);
				sumNumTopic = new int[numTopics];

				if (type == SamplingType.FASTLDA) {
					sortedWordTopic = new DenseInt2DArray(numWords, numTopics);
					invSortedWordTopic = new DenseInt2DArray(numWords, numTopics);
					sortedTopic = new int[numTopics];
					invSortedTopic = new int[numTopics];
					bNormSquared = new double[numWords];
				}
			}
		}

		public Int2DArray getNumWordTopic() {
			return numWordTopic;
		}

		protected void initFastLDA() {
			Integer[] curSortedTopic = new Integer[numTopics];

			// initialize sortedTopic
			for (int i = 0; i < numTopics; i++)
				curSortedTopic[i] = i;
			Arrays.sort(curSortedTopic, new Comparator<Integer>() {
				@Override
				public int compare(Integer idx1, Integer idx2) {
					return Integer.compare(sumNumTopic[idx2], sumNumTopic[idx1]);
				} });
			for (int i = 0; i < numTopics; i++) {
				sortedTopic[i] = curSortedTopic[i];
				invSortedTopic[curSortedTopic[i]] = i;
			}

			// initialize sortedWordTopic
			for (int i = 0; i < numWords; i++) {
				for (int j = 0; j < numTopics; j++)
					curSortedTopic[j] = j;
				final int word = i;
				Arrays.sort(curSortedTopic, new Comparator<Integer>() {
					@Override
					public int compare(Integer idx1, Integer idx2) {
						return Integer.compare(numWordTopic.get(word, idx2), numWordTopic.get(word, idx1));
					} });
				for (int j = 0; j < numTopics; j++) {
					sortedWordTopic.set(i, j, curSortedTopic[j]);
					invSortedWordTopic.set(i, curSortedTopic[j], j);
				}
			}

			// initialize bNormSquared
			for (int i = 0; i < numWords; i++) {
				for (int j = 0; j < numTopics; j++) {
					double v = numWordTopic.get(i, j) + beta.get(j, i);
					bNormSquared[i] += v * v;
				}
			}
		}

		protected void updateSortInc(int incTopic, int incWord) {
			for (int k = invSortedTopic[incTopic]; k > 0; k--) {
				int curIndex = sortedTopic[k];
				int prevIndex = sortedTopic[k - 1];
				if (sumNumTopic[curIndex] <= sumNumTopic[prevIndex])
					break;
				sortedTopic[k] = prevIndex;
				sortedTopic[k - 1] = curIndex;
				invSortedTopic[prevIndex] = k;
				invSortedTopic[curIndex] = k - 1;
			}

			for (int k = invSortedWordTopic.get(incWord, incTopic); k > 0; k--) {
				int curIndex = sortedWordTopic.get(incWord, k);
				int prevIndex = sortedWordTopic.get(incWord, k - 1);
				if (numWordTopic.get(incWord, curIndex) <= numWordTopic.get(incWord, prevIndex))
					break;
				sortedWordTopic.set(incWord, k, prevIndex);
				sortedWordTopic.set(incWord, k - 1, curIndex);
				invSortedWordTopic.set(incWord, prevIndex, k);
				invSortedWordTopic.set(incWord, curIndex, k - 1);
			}
		}

		protected void updateSortDec(int decTopic, int decWord) {
			for (int k = invSortedTopic[decTopic]; k < (numTopics - 1); k++) {
				int curIndex = sortedTopic[k];
				int nextIndex = sortedTopic[k + 1];
				if (sumNumTopic[curIndex] >= sumNumTopic[nextIndex])
					break;
				sortedTopic[k] = nextIndex;
				sortedTopic[k + 1] = curIndex;
				invSortedTopic[nextIndex] = k;
				invSortedTopic[curIndex] = k + 1;
			}

			for (int k = invSortedWordTopic.get(decWord, decTopic); k < (numTopics - 1); k++) {
				int curIndex = sortedWordTopic.get(decWord, k);
				int nextIndex = sortedWordTopic.get(decWord, k + 1);
				if (numWordTopic.get(decWord, curIndex) >= numWordTopic.get(decWord, nextIndex))
					break;
				sortedWordTopic.set(decWord, k, nextIndex);
				sortedWordTopic.set(decWord, k + 1, curIndex);
				invSortedWordTopic.set(decWord, nextIndex, k);
				invSortedWordTopic.set(decWord, curIndex, k + 1);
			}
		}

		protected void trackWordTopicChange(int wordIndex, int topic) {
		}

		protected void compact(boolean keepNumWordTopic) {
			if (!keepNumWordTopic)
				numWordTopic = null;
			if (samplingType == SamplingType.FASTLDA) {
				invSortedWordTopic = null;
				invSortedTopic = null;
			} else
				sumNumTopic = null;
		}

		@Override
		public double[][] estimatePhi() {
			double[][] localPhi = new double[numTopics][numWords];
			for (int i = 0; i < numTopics; i++)
				for (int j = 0; j < numWords; j++)
					localPhi[i][j] = (numWordTopic.get(j, i) + beta.get(i, j)) / (sumNumTopic[i] + sumBeta[i]);
			return localPhi;
		}
	}

	// model parameters
	protected int numWords;
	protected final int numTopics;
	protected double[] alpha;
	protected Double2DArray beta;
	protected double sumAlpha;
	protected double[] sumBeta;

	// sampling parameters

	protected SamplingType samplingType = SamplingType.FASTLDA;

	/** number of iterations for the burn in */
	protected int burnIn = 2000;

	/** after how many iterations a log message is written */
	protected int burnInReportRate = 500;

	/** after how many iterations the model is saved to disk */
	protected int burnInSaveRate = 100;

	protected transient File intermediateStateFile = new File("topic-model-state.ser.gz");

	protected int likelihoodSamples = 100;

	private static final int defaultQuerySamples = 20;
	private static final int defaultQueryBurnIn = 200;

	protected int querySamples = defaultQuerySamples;
	protected int queryLag = 5;

	protected boolean keepNumWordTopic = false;
	protected DiscreteDistributionSampler sampler = new DiscreteDistributionSampler();

	/** uniform priors */
	public LDABasedGibbsSamplingModel(int numWords, int numTopics, double alpha, double beta) {
		this.numWords = numWords;
		this.numTopics = numTopics;
		setUniformAlpha(alpha);
		setUniformBeta(beta);
	}

	/** non-uniform priors */
	public LDABasedGibbsSamplingModel(int numWords, int numTopics, double[] alpha, double[] beta) {
		this.numWords = numWords;
		this.numTopics = numTopics;
		setAlpha(alpha);
		setAllTopicsBeta(beta);
	}

	/** non-uniform priors, separate priors for each topic */
	public LDABasedGibbsSamplingModel(int numWords, int numTopics, double[] alpha, Double2DArray beta) {
		this.numWords = numWords;
		this.numTopics = numTopics;
		setAlpha(alpha);
		setBeta(beta);
	}

	/** copy constructor that copies model and sampling parameters only - called by {@link #createSimilar()} */
	protected LDABasedGibbsSamplingModel(LDABasedGibbsSamplingModel<T> other) {
		this(other.numWords, other.numTopics, other.alpha, other.beta);
		burnIn = other.burnIn;
		burnInReportRate = other.burnInReportRate;
		burnInSaveRate = other.burnInSaveRate;
		likelihoodSamples = other.likelihoodSamples;
		keepNumWordTopic = other.keepNumWordTopic;
		sampler = new DiscreteDistributionSampler(other.sampler);
	}

	public abstract LDABasedGibbsSamplingModel<T> createSimilar();

	public double[] getAlpha() {
		return alpha;
	}

	public void setUniformAlpha(double alpha) {
		this.alpha = new double[numTopics];
		Arrays.fill(this.alpha, alpha);
		sumAlpha = numTopics * alpha;
	}

	public void setAlpha(double[] alpha) {
		this.alpha = alpha;
		sumAlpha = 0.0;
		for (double a : alpha)
			sumAlpha += a;
	}

	public Double2DArray getBeta() {
		return beta;
	}

	public void setUniformBeta(double beta) {
		this.beta = new FixedDouble2DArray(numTopics, numWords, beta);
		sumBeta = new double[numTopics];
		Arrays.fill(sumBeta, numWords * beta);
	}

	public void setAllTopicsBeta(double[] beta) {
		this.beta = new DenseDouble2DArray(numTopics, numWords);
		double sum = 0.0;
		for (int i = 0; i < numWords; i++) {
			sum += beta[i];
			for (int j = 0; j < numTopics; j++)
				this.beta.set(j, i, beta[i]);
		}
		sumBeta = new double[numTopics];
		Arrays.fill(sumBeta, sum);
	}

	public void setBeta(Double2DArray beta) {
		this.beta = new DenseDouble2DArray(numTopics, numWords);
		sumBeta = new double[numTopics];
		for (int i = 0; i < numTopics; i++) {
			sumBeta[i] = 0.0;
			for (int j = 0; j < numWords; j++) {
				double v = beta.get(i, j);
				this.beta.set(i, j, v);
				sumBeta[i] += v;
			}
		}
	}

	public void setSamplingType(SamplingType samplingType) {
		this.samplingType = samplingType;
	}

	public int getBurnIn() {
		return burnIn;
	}

	public void setBurnIn(int burnIn) {
		this.burnIn = burnIn;
	}

	public int getBurnInReportRate() {
		return burnInReportRate;
	}

	public void setBurnInReportRate(int burnInReportRate) {
		this.burnInReportRate = burnInReportRate;
	}

	public int getBurnInSaveRate() {
		return burnInSaveRate;
	}

	public void setBurnInSaveRate(int burnInSaveRate) {
		this.burnInSaveRate = burnInSaveRate;
	}

	public File getIntermediateStateFile() {
		return intermediateStateFile;
	}

	public void setIntermediateStateFile(File intermediateStateFile) {
		this.intermediateStateFile = intermediateStateFile;
	}

	public int getLikelihoodSamples() {
		return likelihoodSamples;
	}

	public void setLikelihoodSamples(int likelihoodSamples) {
		this.likelihoodSamples = likelihoodSamples;
	}

	public int getQuerySamples() {
		return querySamples;
	}

	public void setQuerySamples(int querySamples) {
		this.querySamples = querySamples;
	}

	public int getQueryLag() {
		return queryLag;
	}

	public void setQueryLag(int queryLag) {
		this.queryLag = queryLag;
	}

	public enum ResampleMode { ALL, LOG, NONE };

	public abstract double queryLogLikelihoodLR(Corpus<T> corpus, ResampleMode mode);

	@Override
	public double queryLogLikelihood(Corpus<T> corpus) {
		return queryLogLikelihoodLR(corpus, ResampleMode.ALL);
	}

	/**
	 * @param keepNumWordTopic	If true, keep statistics of word-topic assignment around after training and allow them
	 * 	to be retrieved via {@link SamplingState#getNumWordTopic()}. Needed for Online LDA.
	 */
	public void setKeepNumWordTopic(boolean keepNumWordTopic) {
		this.keepNumWordTopic = keepNumWordTopic;
	}

	public void setSampler(DiscreteDistributionSampler sampler) {
		this.sampler = sampler;
	}

	public abstract SamplingState train(Corpus<T> corpus);
	public abstract SamplingState resumeTraining(Corpus<T> corpus, int numRemainingIter);

	public void configureForQuerying() {
		setQuerySamples(defaultQuerySamples);
		setBurnIn(defaultQueryBurnIn);
		setBurnInReportRate(defaultQueryBurnIn / 2);

		// workaround for loading old models
		if (sampler == null)
			sampler = new DiscreteDistributionSampler();
		if (samplingType == null)
			samplingType = SamplingType.FASTLDA;
	}

}
