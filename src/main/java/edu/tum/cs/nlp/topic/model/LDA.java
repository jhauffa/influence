package edu.tum.cs.nlp.topic.model;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import edu.tum.cs.math.dist.DiscreteDistributionSampler;
import edu.tum.cs.nlp.corpus.Corpus;
import edu.tum.cs.nlp.corpus.Document;
import edu.tum.cs.nlp.corpus.ProcessedDocument;
import edu.tum.cs.util.arrays.Double2DArray;
import edu.tum.cs.util.arrays.Int2DArray;

public class LDA extends LDABasedGibbsSamplingModel<ProcessedDocument> {

	private static final long serialVersionUID = 7204074497992525027L;
	private static final Logger logger = Logger.getLogger(LDA.class.getName());

	public interface QueryResult {
		public double[][] estimateTheta();
	}

	public class MeanQueryResult implements QueryResult {
		public final List<QueryResult> results = new ArrayList<QueryResult>();

		public void add(QueryResult result) {
			results.add(result);
		}

		@Override
		public double[][] estimateTheta() {
			double[][] meanThetas = null;
			double z = 1.0 / results.size();
			for (QueryResult result : results) {
				double[][] thetas = result.estimateTheta();
				if (meanThetas == null)
					meanThetas = new double[thetas.length][thetas[0].length];
				for (int i = 0; i < thetas.length; i++)
					for (int j = 0; j < thetas[i].length; j++)
						meanThetas[i][j] += z * thetas[i][j];
			}
			return meanThetas;
		}
	}

	public class SamplingState extends LDABasedGibbsSamplingModel<ProcessedDocument>.SamplingState
			implements QueryResult {
		private static final long serialVersionUID = -7165937751352963568L;

		protected int[][] z;
		protected transient int[][] numDocumentTopic;

		// FastLDA
		protected transient double[] aNormSquared;

		protected SamplingState(boolean isQuery, SamplingType type) {
			super(isQuery, type);
		}

		public SamplingState(int numDocuments, boolean isQuery, SamplingType type) {
			super(isQuery, type);
			z = new int[numDocuments][];
			numDocumentTopic = new int[numDocuments][numTopics];

			if (type == SamplingType.FASTLDA)
				aNormSquared = new double[numDocuments];
		}

		public QueryResult cloneResults() {
			SamplingState resultState = new SamplingState(true, type);
			resultState.numDocumentTopic = new int[numDocumentTopic.length][];
			for (int i = 0; i < numDocumentTopic.length; i++)
				resultState.numDocumentTopic[i] = numDocumentTopic[i].clone();
			return resultState;
		}

		protected void compact(boolean keepNumWordTopic) {
			super.compact(keepNumWordTopic);
			if (samplingType == SamplingType.FASTLDA)
				aNormSquared = null;
		}

		@Override
		public double[][] estimateTheta() {
			int numDocuments = numDocumentTopic.length;
			double[][] localTheta = new double[numDocuments][numTopics];
			for (int posDocument = 0; posDocument < numDocuments; posDocument++) {
				int sumNumDocumentTopic = 0;
				for (int i = 0; i < numTopics; i++)
					sumNumDocumentTopic += numDocumentTopic[posDocument][i];
				for (int i = 0; i < numTopics; i++)
					localTheta[posDocument][i] = (numDocumentTopic[posDocument][i] + alpha[i]) /
						(sumNumDocumentTopic + sumAlpha);
			}
			return localTheta;
		}

		@Override
		public int[][] getTopicAssignment() {
			return z;
		}
	}

	/** model state after training (for querying) */
	protected SamplingState modelState;

	/** topic-word distributions */
	protected double[][] phi;

	/** uniform priors */
	public LDA(int numWords, int numTopics, double alpha, double beta) {
		super(numWords, numTopics, alpha, beta);
	}

	/** non-uniform priors */
	public LDA(int numWords, int numTopics, double[] alpha, double[] beta) {
		super(numWords, numTopics, alpha, beta);
	}

	/** non-uniform priors, separate priors for each topic */
	public LDA(int numWords, int numTopics, double[] alpha, Double2DArray beta) {
		super(numWords, numTopics, alpha, beta);
	}

	protected LDA(LDA other) {
		super(other);
	}

	@Override
	public LDA createSimilar() {
		return new LDA(this);
	}

	@Override
	public SamplingState train(Corpus<ProcessedDocument> corpus) {
		modelState = burnIn(sampler, corpus, false);
		phi = modelState.estimatePhi();	// estimate the parameters of the multinomial topic-word distributions
		modelState.compact(keepNumWordTopic);
		return modelState;
	}

	@Override
	public SamplingState resumeTraining(Corpus<ProcessedDocument> corpus, int numRemainingIter) {
		throw new UnsupportedOperationException("not implemented yet");
	}

	protected SamplingState burnIn(DiscreteDistributionSampler sampler, Corpus<ProcessedDocument> corpus,
			boolean isQuery) {
		logger.info("Starting sampling with type " + samplingType);
		// random assignment of topics to words
		SamplingState state = new SamplingState(corpus.size(), isQuery, samplingType);
		assignTopicsRandom(sampler, corpus, state);

		// sample until required number of iterations for burn in is reached
		int reportIdx = 0;
		for (int i = 0; i < burnIn; i++) {
			sample(sampler, corpus, state);
			if (++reportIdx == burnInReportRate) {
				logger.info("Burn in step " + (i + 1) + " of " + burnIn);
				reportIdx = 0;
			}
		}
		logger.info("Burn in finished");
		return state;
	}

	protected void assignTopicsRandom(DiscreteDistributionSampler sampler, Corpus<ProcessedDocument> corpus,
			SamplingState state) {
		for (ProcessedDocument document : corpus) {
			int posDocument = document.getIndex();
			state.z[posDocument] = new int[document.size()];

			int posWord = 0;
			for (int wordIndex : document.getWordIds()) {
				int topic = sampler.sample1DUniform(numTopics);
				state.z[posDocument][posWord] = topic;

				state.numDocumentTopic[posDocument][topic]++;
				if (!state.isQuery) {
					state.numWordTopic.adjust(wordIndex, topic, 1);
					state.sumNumTopic[topic]++;
				}

				posWord++;
			}
		}

		// initialize data structures for FastLDA
		int numDocuments = corpus.size();
		if (state.type == SamplingType.FASTLDA) {
			if (!state.isQuery)
				state.initFastLDA();

			// initialize aNormSquared
			for (int i = 0; i < numDocuments; i++) {
				for (int j = 0; j < numTopics; j++) {
					double v = state.numDocumentTopic[i][j] + alpha[j];
					state.aNormSquared[i] += v * v;
				}
			}
		}
	}

	protected void sample(DiscreteDistributionSampler sampler, Corpus<ProcessedDocument> corpus, SamplingState state) {
		for (ProcessedDocument document : corpus)
			sampleDocument(sampler, document, state);
	}

	private void sampleDocument(DiscreteDistributionSampler sampler, ProcessedDocument document, SamplingState state) {
		switch (state.type) {
		case REGULAR:
			sampleDocumentRegular(sampler, document, state);
			break;
		case FASTLDA:
			sampleDocumentFastLda(sampler, document, state);
			break;
		}
	}

	private int sampleTopicAssignment(int wordIndex, int[] numTopic, Int2DArray numWordTopic, int[] sumNumTopic,
			double[] cumulTopicProbs, boolean isQuery) {
		double prevCumulProb = 0.0;
		for (int i = 0; i < numTopics; i++) {
			if (!isQuery) {
				cumulTopicProbs[i] = prevCumulProb + ((numTopic[i] + alpha[i]) *
						(numWordTopic.get(wordIndex, i) + beta.get(i, wordIndex))) / (sumNumTopic[i] + sumBeta[i]);
			} else {
				cumulTopicProbs[i] = prevCumulProb + (numTopic[i] + alpha[i]) * phi[i][wordIndex];
			}

			prevCumulProb = cumulTopicProbs[i];
		}
		return sampler.sample1D(cumulTopicProbs);
	}

	private int sampleTopicAssignmentQuery(int wordIndex, int[] numTopic, double[] cumulTopicProbs) {
		return sampleTopicAssignment(wordIndex, numTopic, null, null, cumulTopicProbs, true);
	}

	private void sampleDocumentRegular(DiscreteDistributionSampler sampler, ProcessedDocument document,
			SamplingState state) {
		double[] cumulTopicProbs = new double[numTopics];
		int posDocument = document.getIndex();
		int posWord = 0;
		for (int wordIndex : document.getWordIds()) {
			// do not count currently assigned topic
			int topic = state.z[posDocument][posWord];
			state.numDocumentTopic[posDocument][topic]--;
			if (!state.isQuery) {
				state.numWordTopic.adjust(wordIndex, topic, -1);
				state.sumNumTopic[topic]--;
			}

			// sample new topic
			topic = sampleTopicAssignment(wordIndex, state.numDocumentTopic[posDocument], state.numWordTopic,
					state.sumNumTopic, cumulTopicProbs, state.isQuery);

			// assign new topic and update counts
			state.z[posDocument][posWord] = topic;
			state.numDocumentTopic[posDocument][topic]++;
			if (!state.isQuery) {
				state.numWordTopic.adjust(wordIndex, topic, 1);
				state.sumNumTopic[topic]++;
			}
			posWord++;
		}
	}

	private void sampleDocumentFastLda(DiscreteDistributionSampler sampler, ProcessedDocument document,
			SamplingState state) {
		double[] cumulTopicProbs = new double[numTopics];
		int posDocument = document.getIndex();
		int posWord = 0;
		for (int wordIndex : document.getWordIds()) {
			// do not count currently assigned topic
			int topic = state.z[posDocument][posWord];
			state.numDocumentTopic[posDocument][topic]--;
			state.aNormSquared[posDocument] -=
					(2.0 * (state.numDocumentTopic[posDocument][topic] + alpha[topic])) + 1.0;
			SamplingState topicState = state;
			if (!state.isQuery) {
				state.numWordTopic.adjust(wordIndex, topic, -1);
				state.sumNumTopic[topic]--;
				state.bNormSquared[wordIndex] -=
						(2.0 * (state.numWordTopic.get(wordIndex, topic) + beta.get(topic, wordIndex))) + 1.0;
				state.updateSortDec(topic, wordIndex);
				state.trackWordTopicChange(wordIndex, topic);
			} else
				topicState = modelState;

			// Actual FastLDA algorithm: compute series of upper bounds on the cumulative topic probability to
			// efficiently determine which topic a sample from a uniform distribution corresponds to.
			// The code corresponds to algorithm 4.1 in Porteous et al., 2008.
			double u = sampler.prng.nextDouble();
			double prevCumulProb = 0.0;
			double upperBound = 0.0, prevUpperBound;
			double aNormSquared = state.aNormSquared[posDocument];
			double bNormSquared = topicState.bNormSquared[wordIndex];
	sampleTopic:
			for (int i = 0; i < numTopics; i++) {
				int t = topicState.sortedWordTopic.get(wordIndex, i);

				// compute cumulative probabilities up to current topic
				double docTerm = state.numDocumentTopic[posDocument][t] + alpha[t];
				double topicTerm;
				if (!state.isQuery) {
					topicTerm = state.numWordTopic.get(wordIndex, t) + beta.get(t, wordIndex);
					cumulTopicProbs[i] = prevCumulProb + (docTerm *
							(topicTerm / (state.sumNumTopic[t] + sumBeta[t])));
				} else {
					topicTerm = phi[t][wordIndex] * (topicState.sumNumTopic[t] + sumBeta[t]);
					cumulTopicProbs[i] = prevCumulProb + (docTerm * phi[t][wordIndex]);
				}
				prevCumulProb = cumulTopicProbs[i];

				// compute new upper bound
				prevUpperBound = upperBound;
				aNormSquared -= docTerm * docTerm;
				if (aNormSquared < 0.0)
					aNormSquared = 0.0;
				bNormSquared -= topicTerm * topicTerm;
				if (bNormSquared < 0.0)
					bNormSquared = 0.0;

				double abNorm = Math.sqrt(aNormSquared * bNormSquared);
				if ((i == (numTopics - 1)) && (abNorm > 0.01))
					throw new RuntimeException("bad upper bound, a norm = " + aNormSquared +
							", b norm = " + bNormSquared);
				int minTopic = topicState.sortedTopic[numTopics - 1];
				double cNorm = 1.0 / (topicState.sumNumTopic[minTopic] + sumBeta[minTopic]);
				upperBound = cumulTopicProbs[i] + (abNorm * cNorm);

				// sample topic
				double p = u * upperBound;
				if (p > cumulTopicProbs[i])
					continue;
				if ((i == 0) || (p > cumulTopicProbs[i - 1])) {
					topic = t;
					break;
				}
				p = (((u * prevUpperBound) - cumulTopicProbs[i - 1]) * upperBound) / (prevUpperBound - upperBound);
				for (int j = 0; j < i; j++) {
					if (p <= cumulTopicProbs[j]) {
						topic = topicState.sortedWordTopic.get(wordIndex, j);
						break sampleTopic;
					}
				}
			}

			// assign new topic and update counts
			state.z[posDocument][posWord] = topic;
			state.numDocumentTopic[posDocument][topic]++;
			state.aNormSquared[posDocument] +=
					(2.0 * (state.numDocumentTopic[posDocument][topic] + alpha[topic])) - 1.0;
			if (!state.isQuery) {
				state.numWordTopic.adjust(wordIndex, topic, 1);
				state.sumNumTopic[topic]++;
				state.bNormSquared[wordIndex] +=
						(2.0 * (state.numWordTopic.get(wordIndex, topic) + beta.get(topic, wordIndex))) - 1.0;
				state.updateSortInc(topic, wordIndex);
				state.trackWordTopicChange(wordIndex, topic);
			}

			posWord++;
		}
	}

	@Override
	public double[][] getTopicWordDistr() {
		return phi;
	}

	public QueryResult query(Corpus<ProcessedDocument> corpus) {
		SamplingState state = burnIn(sampler, corpus, true);
		if (querySamples == 1)
			return state;

		MeanQueryResult result = new MeanQueryResult();
		result.add(state);
		for (int i = 1; i < querySamples; i++) {
			for (int j = 0; j < queryLag; j++)
				sample(sampler, corpus, state);
			result.add(state.cloneResults());
		}
		return result;
	}

	private void sampleSingleParticleLR(ProcessedDocument document, double[] tokenProb, ResampleMode mode) {
		int n = document.size();
		int[] y = document.getWordIds();
		int[] z = new int[n];
		int[] numTopic = new int[numTopics];
		double[] cumulTopicProbs = new double[numTopics];

		for (int i = 0; i < n; i++) {
			// resample topic assignments of previous words
			if (mode == ResampleMode.LOG) {
				if (i > 0) {
					int numResample = Math.max(1, (int) (Math.log(n) * i / n));
					for (int j = 0; j < numResample; j++) {
						int idx = sampler.prng.nextInt(i);
						numTopic[z[idx]]--;
						z[idx] = sampleTopicAssignmentQuery(y[idx], numTopic, cumulTopicProbs);
						numTopic[z[idx]]++;
					}
				}
			} else if (mode == ResampleMode.ALL) {
				for (int j = 0; j < i; j++) {
					numTopic[z[j]]--;
					z[j] = sampleTopicAssignmentQuery(y[j], numTopic, cumulTopicProbs);
					numTopic[z[j]]++;
				}
			}	// skip resampling in ResampleMode.NONE

			// compute predictive probability of current word, given the topic assignments of previous words
			for (int k = 0; k < numTopics; k++)
				tokenProb[i] += (numTopic[k] + alpha[k]) / (i + sumAlpha) * phi[k][y[i]];
			if (i == (n - 1))
				break;

			// sample topic assignment of current word
			z[i] = sampleTopicAssignmentQuery(y[i], numTopic, cumulTopicProbs);
			numTopic[z[i]]++;
		}
	}

	@Override
	public double queryLogLikelihoodLR(Corpus<ProcessedDocument> corpus, ResampleMode mode) {
		int numParticles = likelihoodSamples;
		double logLikelihood = 0.0;
		for (ProcessedDocument document : corpus) {
			double llDoc = 0.0;
			double[][] tokenProb = new double[numParticles][document.size()];
			for (int i = 0; i < numParticles; i++)
				sampleSingleParticleLR(document, tokenProb[i], mode);
			for (int i = 0; i < document.size(); i++) {
				double sum = 0.0;
				for (int j = 0; j < numParticles; j++)
					sum += tokenProb[j][i];
				llDoc += Math.log(sum) - Math.log(numParticles);
			}
			logLikelihood += llDoc;
		}
		return logLikelihood;
	}

	public double queryApproxLogLikelihood(Corpus<ProcessedDocument> corpus) {
		double logLikelihood = 0.0;
		double[] w = new double[numTopics];

		for (ProcessedDocument document : corpus) {
			double[] z = new double[numTopics];
			double sumZ = sumAlpha;

			for (int wordIndex : document.getWordIds()) {
				double sumW = 0.0;
				for (int k = 0; k < numTopics; k++) {
					w[k] = ((alpha[k] + z[k]) / sumZ) * phi[k][wordIndex];
					sumW += w[k];
				}
				logLikelihood += Math.log(sumW);

				for (int k = 0; k < numTopics; k++) {
					double wNorm = w[k] / sumW;
					z[k] += wNorm;
					sumZ += wNorm;
				}
			}
		}
		return logLikelihood;
	}

	public static <T extends Document> double estimateLogLikelihood(Corpus<ProcessedDocument> corpus, double[][] phi,
			double[][] theta) {
		double logLikelihood = 0.0;
		for (ProcessedDocument document : corpus)
			logLikelihood += estimateLogLikelihood(document, phi, theta[document.getIndex()]);
		return logLikelihood;
	}

	public static <T extends Document> double estimateLogLikelihood(ProcessedDocument document, double[][] phi,
			double[] theta) {
		double logLikelihood = 0.0;
		for (int wordIndex : document.getWordIds()) {
			double sumProb = 0.0;
			for (int i = 0; i < phi.length; i++)
				sumProb += theta[i] * phi[i][wordIndex];
			logLikelihood += Math.log(sumProb);
		}
		return logLikelihood;
	}

}
