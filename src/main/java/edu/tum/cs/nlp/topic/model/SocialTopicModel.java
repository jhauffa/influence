package edu.tum.cs.nlp.topic.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import be.ac.ulg.montefiore.run.jahmm.LogSpace;
import edu.tum.cs.math.dist.DiscreteDistributionSampler;
import edu.tum.cs.nlp.corpus.Corpus;
import edu.tum.cs.nlp.corpus.MessageSequenceCorpus;
import edu.tum.cs.nlp.corpus.ProcessedMessage;
import edu.tum.cs.util.arrays.Double2DArray;
import edu.tum.cs.util.arrays.Int2DArray;
import edu.tum.cs.util.arrays.Object2DArray;
import edu.tum.cs.util.arrays.Sparse2DIterator;
import edu.tum.cs.util.arrays.SparseObject2DArray;

public class SocialTopicModel extends LDABasedGibbsSamplingModel<ProcessedMessage> {

	private static final long serialVersionUID = 1848390119860180988L;
	private static final Logger logger = Logger.getLogger(SocialTopicModel.class.getName());

	protected static class PerAuthorRecipientData {
		public final int[][] numTopic;
		public int[] sumNumTopic;
		public int[][] numTransition;
		public int[] sumNumTransition;
		public double[] aNormSquared;	// FastLDA

		public PerAuthorRecipientData(int numStates, int numTopics) {
			numTopic = new int[numStates][numTopics];
			sumNumTopic = new int[numStates];
			numTransition = new int[numStates][numStates];
			sumNumTransition = new int[numStates];
			aNormSquared = new double[numStates];
		}

		public PerAuthorRecipientData(PerAuthorRecipientData other) {
			numTopic = new int[other.numTopic.length][];
			for (int i = 0; i < other.numTopic.length; i++)
				numTopic[i] = other.numTopic[i].clone();
			sumNumTopic = other.sumNumTopic.clone();
			numTransition = new int[other.numTransition.length][];
			for (int i = 0; i < other.numTransition.length; i++)
				numTransition[i] = other.numTransition[i].clone();
			sumNumTransition = other.sumNumTransition.clone();
			aNormSquared = other.aNormSquared.clone();
		}
	}

	public class SamplingState extends LDABasedGibbsSamplingModel<ProcessedMessage>.SamplingState {
		private static final long serialVersionUID = -650913939177097270L;

		// latent variable assignments
		protected int[][] z;
		protected int[] c;

		// counts
		protected transient SparseObject2DArray<PerAuthorRecipientData> perAuthorRecipientData;

		protected SamplingState(boolean isQuery, SamplingType type) {
			super(isQuery, type);
		}

		public SamplingState(int numDocuments, boolean isQuery, SamplingType type) {
			super(isQuery, type);

			z = new int[numDocuments][];
			c = new int[numDocuments];

			perAuthorRecipientData = new SparseObject2DArray<PerAuthorRecipientData>(numPersons, numPersons);
		}

		private PerAuthorRecipientData initPerAuthorRecipientData(int author, int recipient) {
			PerAuthorRecipientData ar = perAuthorRecipientData.get(author, recipient);
			if (ar == null) {
				ar = new PerAuthorRecipientData(numStates, numTopics);
				perAuthorRecipientData.set(author, recipient, ar);
			}
			return ar;
		}

		private void initANormSquared() {
			Sparse2DIterator<PerAuthorRecipientData> it = perAuthorRecipientData.sparseIterator();
			while (it.hasNext()) {
				it.advance();
				PerAuthorRecipientData ar = it.getValue();
				for (int i = 0; i < numStates; i++) {
					ar.aNormSquared[i] = 0.0;
					for (int j = 0; j < numTopics; j++) {
						double v = ar.numTopic[i][j] + alpha[j];
						ar.aNormSquared[i] += v * v;
					}
				}
			}
		}

		public Object2DArray<double[][]> estimateTheta() {
			Object2DArray<double[][]> theta = new SparseObject2DArray<>(numPersons, numPersons);

			Sparse2DIterator<PerAuthorRecipientData> it = perAuthorRecipientData.sparseIterator();
			while (it.hasNext()) {
				it.advance();
				PerAuthorRecipientData ar = it.getValue();
				double[][] thetaAutorRecipient = new double[numStates][numTopics];
				for (int i = 0; i < numStates; i++) {
					for (int j = 0; j < numTopics; j++)
						thetaAutorRecipient[i][j] = (ar.numTopic[i][j] + alpha[j]) / (ar.sumNumTopic[i] + sumAlpha);
				}
				theta.set(it.getRow(), it.getColumn(), thetaAutorRecipient);
			}

			return theta;
		}

		/** estimate Markov state transition probabilities */
		public Object2DArray<double[][]> estimatePi() {
			Object2DArray<double[][]> pi = new SparseObject2DArray<>(numPersons, numPersons);

			Sparse2DIterator<PerAuthorRecipientData> it = perAuthorRecipientData.sparseIterator();
			while (it.hasNext()) {
				it.advance();
				PerAuthorRecipientData ar = it.getValue();
				double[][] piAuthorRecipient = new double[numStates][numStates];
				for (int i = 0; i < numStates; i++) {
					for (int j = 0; j < numStates; j++) {
						piAuthorRecipient[i][j] = (ar.numTransition[i][j] + gamma[i][j]) /
							(ar.sumNumTransition[i] + sumGamma[i]);
					}
				}
				pi.set(it.getRow(), it.getColumn(), piAuthorRecipient);
			}

			return pi;
		}

		@Override
		public int[][] getTopicAssignment() {
			return z;
		}

		public int[] getMarkovStateAssignment() {
			return c;
		}
	}

	protected int numPersons, numStates;

	/** Dirichlet prior of Markov state transition probabilities */
	private double[][] gamma;

	private double[] sumGamma;

	/** model state after training (for querying) */
	protected SamplingState modelState;

	/** topic-word distributions */
	protected double[][] phi;

	/** uniform priors */
	public SocialTopicModel(int numWords, int numPersons, int numTopics, int numStates, double alpha,
			double beta, double gamma) {
		super(numWords, numTopics, alpha, beta);
		this.numPersons = numPersons;
		this.numStates = numStates;
		setUniformGamma(gamma);
	}

	/** uniform alpha, beta + non-uniform gamma */
	public SocialTopicModel(int numWords, int numPersons, int numTopics, int numStates, double alpha,
			double beta, double[][] gamma) {
		super(numWords, numTopics, alpha, beta);
		this.numPersons = numPersons;
		this.numStates = numStates;
		setGamma(gamma);
	}

	/** non-uniform priors */
	public SocialTopicModel(int numWords, int numPersons, int numTopics, int numStates, double[] alpha,
			double[] beta, double[][] gamma) {
		super(numWords, numTopics, alpha, beta);
		this.numPersons = numPersons;
		this.numStates = numStates;
		setGamma(gamma);
	}

	/** non-uniform priors, separate priors for each topic */
	public SocialTopicModel(int numWords, int numPersons, int numTopics, int numStates, double[] alpha,
			Double2DArray beta, double[][] gamma) {
		super(numWords, numTopics, alpha, beta);
		this.numPersons = numPersons;
		this.numStates = numStates;
		setGamma(gamma);
	}

	protected SocialTopicModel(SocialTopicModel other) {
		super(other);
		numPersons = other.numPersons;
		numStates = other.numStates;
		setGamma(other.gamma);
	}

	public double[][] getGamma() {
		return gamma;
	}

	public void setUniformGamma(double gamma) {
		this.gamma = new double[numStates][numStates];
		for (double[] g : this.gamma)
			Arrays.fill(g, gamma);
		sumGamma = new double[numStates];
		Arrays.fill(sumGamma, numStates * gamma);
	}

	public void setGamma(double[][] gamma) {
		this.gamma = gamma;
		sumGamma = new double[numStates];
		for (int i = 0; i < numStates; i++) {
			for (double g : gamma[i])
				sumGamma[i] += g;
		}
	}

	@Override
	public SocialTopicModel createSimilar() {
		return new SocialTopicModel(this);
	}

	@Override
	public SamplingState train(Corpus<ProcessedMessage> corpus) {
		// random assignment of words and states to topics
		modelState = new SamplingState(corpus.size(), false, samplingType);
		initSamplingState(sampler, (MessageSequenceCorpus<?, ?>) corpus, modelState, true);

		// sample until required number of iterations for burn in is reached
		burnIn(sampler, corpus, modelState, burnIn);
		phi = modelState.estimatePhi();	// estimate the parameters of the multinomial topic-word distributions
		modelState.compact(keepNumWordTopic);
		return modelState;
	}

	@Override
	public SamplingState resumeTraining(Corpus<ProcessedMessage> corpus, int numRemainingIter) {
		SamplingState assignmentState = modelState;
		modelState = new SamplingState(corpus.size(), false, samplingType);
		modelState.z = assignmentState.z;
		modelState.c = assignmentState.c;
		initSamplingState(sampler, (MessageSequenceCorpus<?, ?>) corpus, modelState, false);

		burnIn(sampler, corpus, modelState, numRemainingIter);
		phi = modelState.estimatePhi();
		modelState.compact(keepNumWordTopic);
		return modelState;
	}

	protected void burnIn(DiscreteDistributionSampler sampler, Corpus<ProcessedMessage> corpus, SamplingState state,
			int numIter) {
		logger.info("Starting sampling");
		int reportIdx = 0;
		for (int i = 0; i < numIter; i++) {
			sample(sampler, corpus, state);
			if (++reportIdx == burnInReportRate) {
				logger.info("Burn in step " + (i + 1) + " of " + numIter);
				reportIdx = 0;
			}
		}
		logger.info("Burn in finished");
	}

	protected void initSamplingState(DiscreteDistributionSampler sampler, MessageSequenceCorpus<?, ?> corpus,
			SamplingState state, boolean assignRandom) {
		Iterator<MessageSequenceCorpus.MessageSequence> seqIt = corpus.sequenceIterator();
		while (seqIt.hasNext()) {
			MessageSequenceCorpus.MessageSequence seq = seqIt.next();
			while (seq.hasNext()) {
				ProcessedMessage message = seq.next();
				int posMessage = message.getIndex();
				if (message.getRecipientIds().length != 1) {
					throw new RuntimeException("message " + message.getIndex() +
							" has an invalid number of recipients: " + message.getRecipientIds().length);
				}
				PerAuthorRecipientData ar = state.initPerAuthorRecipientData(message.getSenderId(),
						message.getRecipientIds()[0]);
				int markovState;
				if (assignRandom) {
					state.z[posMessage] = new int[message.size()];
					markovState = sampler.sample1DUniform(numStates);
					state.c[posMessage] = markovState;
				} else
					markovState = state.c[posMessage];

				int posPrevMessage = message.getPredecessor();	// index of previous message of author-recipient pair
				int prevMarkovState = (posPrevMessage >= 0) ? state.c[posPrevMessage] : 0;
				ar.numTransition[prevMarkovState][markovState]++;
				ar.sumNumTransition[prevMarkovState]++;
				if (message.getSuccessor() < 0) {
					// transition to final state if message is last in chain
					ar.numTransition[markovState][0]++;
					ar.sumNumTransition[markovState]++;
				}

				int posWord = 0;
				for (int wordIndex : message.getWordIds()) {
					int topic;
					if (assignRandom) {
						topic = sampler.sample1DUniform(numTopics);
						state.z[posMessage][posWord] = topic;
					} else
						topic = state.z[posMessage][posWord];

					ar.numTopic[markovState][topic]++;
					ar.sumNumTopic[markovState]++;

					if (!state.isQuery) {
						state.numWordTopic.adjust(wordIndex, topic, 1);
						state.sumNumTopic[topic]++;
					}

					posWord++;
				}
			}
		}

		// initialize data structures for FastLDA
		if (state.type == SamplingType.FASTLDA) {
			if (!state.isQuery)
				state.initFastLDA();
			state.initANormSquared();
		}
	}

	protected void sample(DiscreteDistributionSampler sampler, Corpus<ProcessedMessage> corpus, SamplingState state) {
		for (ProcessedMessage message : corpus) {
			PerAuthorRecipientData ar = state.perAuthorRecipientData.get(message.getSenderId(),
					message.getRecipientIds()[0]);
			if (state.type == SamplingType.FASTLDA)
				sampleTopicsFastLda(sampler, message, ar, state);
			else
				sampleTopicsRegular(sampler, message, ar, state);
		}
	}

	private void sampleTopicsFastLda(DiscreteDistributionSampler sampler, ProcessedMessage message,
			PerAuthorRecipientData ar, SamplingState state) {
		double[] cumulTopicProbs = new double[numTopics];
		double[] cumulStateProbs = new double[numStates];

		int posPrevMessage = message.getPredecessor();
		int prevMarkovState = (posPrevMessage >= 0) ? state.c[posPrevMessage] : 0;
		int posMessage = message.getIndex();
		int markovState = state.c[posMessage];
		int posNextMessage = message.getSuccessor();
		int nextMarkovState = (posNextMessage >= 0) ? state.c[posNextMessage] : 0;

		double[] pTransBase = new double[numStates];
		for (int i = 0; i < numStates; i++) {
			int incTransition = ((prevMarkovState == i) && (i == nextMarkovState)) ? 1 : 0;
			int incState = (prevMarkovState == i) ? 1 : 0;
			pTransBase[i] = (ar.numTransition[prevMarkovState][i] + gamma[prevMarkovState][i]) *
					(ar.numTransition[i][nextMarkovState] + incTransition + gamma[i][nextMarkovState]) /
					(ar.sumNumTransition[i] + incState + sumGamma[i]);
		}

		int posWord = 0;
		for (int wordIndex : message.getWordIds()) {
			int topic = state.z[posMessage][posWord];

			ar.numTopic[markovState][topic]--;
			ar.sumNumTopic[markovState]--;
			ar.aNormSquared[markovState] -= (2.0 * (ar.numTopic[markovState][topic] + alpha[topic])) + 1.0;
			ar.numTransition[prevMarkovState][markovState]--;
			ar.sumNumTransition[prevMarkovState]--;
			ar.numTransition[markovState][nextMarkovState]--;
			ar.sumNumTransition[markovState]--;

			// sample Markov state
			double prevCumulProb = 0.0;
			for (int i = 0; i < numStates; i++) {
				double p = pTransBase[i] * (ar.numTopic[i][topic] + alpha[topic]) / (ar.sumNumTopic[i] + sumAlpha);
				cumulStateProbs[i] = prevCumulProb + p;
				prevCumulProb = cumulStateProbs[i];
			}
			int oldMarkovState = markovState;
			markovState = sampler.sample1D(cumulStateProbs);

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
			prevCumulProb = 0.0;
			double upperBound = 0.0, prevUpperBound;
			double aNormSquared = ar.aNormSquared[markovState];
			double bNormSquared = topicState.bNormSquared[wordIndex];
	sampleTopic:
			for (int i = 0; i < numTopics; i++) {
				int t = topicState.sortedWordTopic.get(wordIndex, i);

				// compute cumulative probabilities up to current topic
				double messageTerm = ar.numTopic[markovState][t] + alpha[t];
				double topicTerm;
				if (!state.isQuery) {
					topicTerm = state.numWordTopic.get(wordIndex, t) + beta.get(t, wordIndex);
					cumulTopicProbs[i] = prevCumulProb + (messageTerm *
							(topicTerm / (state.sumNumTopic[t] + sumBeta[t])));
				} else {
					topicTerm = phi[t][wordIndex] * (topicState.sumNumTopic[t] + sumBeta[t]);
					cumulTopicProbs[i] = prevCumulProb + (messageTerm * phi[t][wordIndex]);
				}
				prevCumulProb = cumulTopicProbs[i];

				// compute new upper bound
				prevUpperBound = upperBound;
				aNormSquared -= messageTerm * messageTerm;
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
			state.z[posMessage][posWord] = topic;
			state.c[posMessage] = markovState;

			ar.numTopic[markovState][topic]++;
			ar.sumNumTopic[markovState]++;
			ar.aNormSquared[markovState] += (2.0 * (ar.numTopic[markovState][topic] + alpha[topic])) - 1.0;
			ar.numTransition[prevMarkovState][markovState]++;
			ar.sumNumTransition[prevMarkovState]++;
			ar.numTransition[markovState][nextMarkovState]++;
			ar.sumNumTransition[markovState]++;
			if (!state.isQuery) {
				state.numWordTopic.adjust(wordIndex, topic, 1);
				state.sumNumTopic[topic]++;
				state.bNormSquared[wordIndex] +=
						(2.0 * (state.numWordTopic.get(wordIndex, topic) + beta.get(topic, wordIndex))) - 1.0;
				state.updateSortInc(topic, wordIndex);
				state.trackWordTopicChange(wordIndex, topic);
			}

			reassignMessageState(state, ar, state.z[posMessage], posWord, oldMarkovState, markovState);
			posWord++;
		}
	}

	private int sampleStateTopicAssignment(int wordIndex, int prevMarkovState, int nextMarkovState,
			PerAuthorRecipientData ar, Int2DArray numWordTopic, int[] sumNumTopic, double[] cumulTopicProbs,
			boolean isQuery) {
		double prevCumulProb = 0.0;
		int idx = 0;
		for (int i = 0; i < numStates; i++) {
			int incTransition = ((prevMarkovState == i) && (i == nextMarkovState)) ? 1 : 0;
			int incState = (prevMarkovState == i) ? 1 : 0;
			double pTrans = (ar.numTransition[prevMarkovState][i] + gamma[prevMarkovState][i]) *
					(ar.numTransition[i][nextMarkovState] + incTransition + gamma[i][nextMarkovState]) /
					(ar.sumNumTransition[i] + incState + sumGamma[i]);

			for (int j = 0; j < numTopics; j++) {
				double p = pTrans * (ar.numTopic[i][j] + alpha[j]) / (ar.sumNumTopic[i] + sumAlpha);
				if (!isQuery)
					p *= (numWordTopic.get(wordIndex, j) + beta.get(j, wordIndex)) / (sumNumTopic[j] + sumBeta[j]);
				else
					p *= phi[j][wordIndex];

				cumulTopicProbs[idx] = prevCumulProb + p;
				prevCumulProb = cumulTopicProbs[idx];
				idx++;
			}
		}
		return sampler.sample1D(cumulTopicProbs);
	}

	private int sampleStateTopicAssignmentQuery(int wordIndex, int prevMarkovState, int nextMarkovState,
			PerAuthorRecipientData ar, double[] cumulTopicProbs) {
		return sampleStateTopicAssignment(wordIndex, prevMarkovState, nextMarkovState, ar, null, null, cumulTopicProbs,
				true);
	}

	private void sampleTopicsRegular(DiscreteDistributionSampler sampler, ProcessedMessage message,
			PerAuthorRecipientData ar, SamplingState state) {
		int posPrevMessage = message.getPredecessor();
		int prevMarkovState = (posPrevMessage >= 0) ? state.c[posPrevMessage] : 0;
		int posMessage = message.getIndex();
		int markovState = state.c[posMessage];
		int posNextMessage = message.getSuccessor();
		int nextMarkovState = (posNextMessage >= 0) ? state.c[posNextMessage] : 0;

		double[] cumulTopicProbs = new double[numTopics * numStates];
		int posWord = 0;
		for (int wordIndex : message.getWordIds()) {
			// do not count currently assigned topic and Markov state
			int topic = state.z[posMessage][posWord];
			ar.numTopic[markovState][topic]--;
			ar.sumNumTopic[markovState]--;
			ar.numTransition[prevMarkovState][markovState]--;
			ar.sumNumTransition[prevMarkovState]--;
			ar.numTransition[markovState][nextMarkovState]--;
			ar.sumNumTransition[markovState]--;
			if (!state.isQuery) {
				state.numWordTopic.adjust(wordIndex, topic, -1);
				state.sumNumTopic[topic]--;
			}

			// sample new topic and Markov state
			int v = sampleStateTopicAssignment(wordIndex, prevMarkovState, nextMarkovState, ar, state.numWordTopic,
					state.sumNumTopic, cumulTopicProbs, state.isQuery);
			int oldMarkovState = markovState;
			markovState = v / numTopics;
			topic = v % numTopics;
			state.z[posMessage][posWord] = topic;
			state.c[posMessage] = markovState;

			ar.numTopic[markovState][topic]++;
			ar.sumNumTopic[markovState]++;
			ar.numTransition[prevMarkovState][markovState]++;
			ar.sumNumTransition[prevMarkovState]++;
			ar.numTransition[markovState][nextMarkovState]++;
			ar.sumNumTransition[markovState]++;
			if (!state.isQuery) {
				state.numWordTopic.adjust(wordIndex, topic, 1);
				state.sumNumTopic[topic]++;
			}

			reassignMessageState(state, ar, state.z[posMessage], posWord, oldMarkovState, markovState);
			posWord++;
		}
	}

	private void reassignMessageState(SamplingState state, PerAuthorRecipientData ar, int[] z, int posWord,
			int oldMarkovState, int markovState) {
		if (oldMarkovState != markovState) {
			for (int i = 0; i < z.length; i++) {
				if (i == posWord)
					continue;

				int curTopic = z[i];
				ar.numTopic[oldMarkovState][curTopic]--;
				ar.sumNumTopic[oldMarkovState]--;
				ar.numTopic[markovState][curTopic]++;
				ar.sumNumTopic[markovState]++;
				if (state.type == SamplingType.FASTLDA) {
					ar.aNormSquared[oldMarkovState] -=
							(2.0 * (ar.numTopic[oldMarkovState][curTopic] + alpha[curTopic])) + 1.0;
					ar.aNormSquared[markovState] +=
							(2.0 * (ar.numTopic[markovState][curTopic] + alpha[curTopic])) - 1.0;
				}
			}
		}
	}

	@Override
	public double[][] getTopicWordDistr() {
		return phi;
	}

	private void sampleWordLR(int posMessage, int posWord, int wordIndex, int[] c, int[] zMsg, int msgLength,
			PerAuthorRecipientData ar, double[] cumulTopicProbs, boolean isResampling) {
		int topic;
		// c is initialized to 0, so this is correct even if the current/next message has not yet been sampled
		int markovState = c[posMessage];
		int prevMarkovState = (posMessage > 0) ? c[posMessage - 1] : 0;
		int nextMarkovState = (posMessage < (c.length - 1)) ? c[posMessage + 1] : 0;

		if (isResampling) {
			topic = zMsg[posWord];
			ar.numTopic[markovState][topic]--;
			ar.sumNumTopic[markovState]--;
		}
		if (isResampling || (posWord > 0)) {
			ar.numTransition[prevMarkovState][markovState]--;
			ar.sumNumTransition[prevMarkovState]--;
			ar.numTransition[markovState][nextMarkovState]--;
			ar.sumNumTransition[markovState]--;
		} else if (posMessage > 0) {
			ar.numTransition[prevMarkovState][markovState]--;
			ar.sumNumTransition[prevMarkovState]--;
		}

		int v = sampleStateTopicAssignmentQuery(wordIndex, prevMarkovState, nextMarkovState, ar, cumulTopicProbs);
		int oldMarkovState = markovState;
		markovState = v / numTopics;
		topic = v % numTopics;

		zMsg[posWord] = topic;
		c[posMessage] = markovState;

		ar.numTopic[markovState][topic]++;
		ar.sumNumTopic[markovState]++;
		ar.numTransition[prevMarkovState][markovState]++;
		ar.sumNumTransition[prevMarkovState]++;
		ar.numTransition[markovState][nextMarkovState]++;
		ar.sumNumTransition[markovState]++;

		if (oldMarkovState != markovState) {
			// reassign all already sampled words, except posWord
			for (int i = 0; i < msgLength; i++) {
				if (i == posWord)
					continue;

				int curTopic = zMsg[i];
				ar.numTopic[oldMarkovState][curTopic]--;
				ar.sumNumTopic[oldMarkovState]--;
				ar.numTopic[markovState][curTopic]++;
				ar.sumNumTopic[markovState]++;
			}
		}
	}

	private void sampleSingleParticleLR(List<ProcessedMessage> seq, double[] tokenProb, ResampleMode mode) {
		int[][] z = new int[seq.size()][];
		int[] c = new int[seq.size()];
		PerAuthorRecipientData ar = new PerAuthorRecipientData(numStates, numTopics);
		double[] cumulTopicProbs = new double[numTopics * numStates];

		double[] cumulMessageProb = null;
		if (mode == ResampleMode.LOG) {
			cumulMessageProb = new double[seq.size()];
			double prevCumulProb = 0.0;
			int idx = 0;
			for (ProcessedMessage msg : seq) {
				cumulMessageProb[idx] = prevCumulProb + msg.size();
				prevCumulProb = cumulMessageProb[idx];
				idx++;
			}
		}

		int posMessage = 0;
		int pos = 0;
		for (ProcessedMessage msg : seq) {
			z[posMessage] = new int[msg.size()];
			int[] y = msg.getWordIds();

			for (int i = 0; i < y.length; i++) {
				// resample state/topic assignments of previous messages/words
				if (pos > 0) {
					if (mode == ResampleMode.LOG) {
						int n = tokenProb.length;
						int numResample = Math.max(1, (int) (Math.log(n) * pos / n));
						int maxMsgIdx = posMessage + ((i > 0) ? 1 : 0);
						Map<Integer, Integer> msgCounts = new HashMap<Integer, Integer>();
						for (int j = 0; j < numResample; j++) {
							int msgIdx = sampler.sample1D(cumulMessageProb, maxMsgIdx);
							Integer cnt = msgCounts.get(msgIdx);
							if (cnt == null)
								cnt = 0;
							msgCounts.put(msgIdx, cnt + 1);
						}

						for (Map.Entry<Integer, Integer> e : msgCounts.entrySet()) {
							int posMessageRes = e.getKey();
							ProcessedMessage msgRes = seq.get(posMessageRes);
							int[] yRes = msgRes.getWordIds();
							int maxIdx = (posMessageRes == posMessage) ? i : yRes.length;

							for (int j = 0; j < e.getValue(); j++) {
								int idx = sampler.prng.nextInt(maxIdx);
								sampleWordLR(posMessageRes, idx, yRes[idx], c, z[posMessageRes], maxIdx, ar,
										cumulTopicProbs, true);
							}
						}
					} else if (mode == ResampleMode.ALL) {
						int posMessageRes = 0;
resamplingLoop:
						for (ProcessedMessage msgRes : seq) {
							boolean inCurrentMessage = (posMessageRes == posMessage);
							int[] yRes = msgRes.getWordIds();
							for (int j = 0; j < yRes.length; j++) {
								if (inCurrentMessage && (j == i))
									break resamplingLoop;

								int maxIdx = inCurrentMessage ? i : yRes.length;
								sampleWordLR(posMessageRes, j, yRes[j], c, z[posMessageRes], maxIdx, ar,
										cumulTopicProbs, true);
							}
							posMessageRes++;
						}
					}	// otherwise, resampling mode is NONE - skip resampling
				}

				// compute predictive probability of current word, given previous state/topic assignments
				int prevMarkovState = (posMessage > 0) ? c[posMessage - 1] : 0;
				double p = 0.0;
				for (int j = 0; j < numStates; j++) {
					for (int k = 0; k < numTopics; k++) {
						double estPi = (ar.numTransition[prevMarkovState][j] + gamma[prevMarkovState][j]) /
								(ar.sumNumTransition[prevMarkovState] + sumGamma[prevMarkovState]);
						double estTheta = (ar.numTopic[j][k] + alpha[k]) / (ar.sumNumTopic[j] + sumAlpha);
						p += Math.pow(estPi, 1.0 / (i + 1)) * estTheta * phi[k][y[i]];
					}
				}
				tokenProb[pos] += p;

				// sample state and topic assignment of current word
				sampleWordLR(posMessage, i, y[i], c, z[posMessage], i + 1, ar, cumulTopicProbs, false);
				pos++;
			}
			posMessage++;
		}
	}

	@Override
	public double queryLogLikelihoodLR(Corpus<ProcessedMessage> corpus, ResampleMode mode) {
		int numParticles = likelihoodSamples;
		double logLikelihood = 0.0;
		Iterator<MessageSequenceCorpus.MessageSequence> seqIt =
				((MessageSequenceCorpus<?, ?>) corpus).sequenceIterator();
		while (seqIt.hasNext()) {
			MessageSequenceCorpus.MessageSequence seq = seqIt.next();
			List<ProcessedMessage> msgs = new ArrayList<ProcessedMessage>();
			int numTokens = 0;
			while (seq.hasNext()) {
				ProcessedMessage msg = seq.next();
				msgs.add(msg);
				numTokens += msg.size();
			}

			double llSeq = 0.0;
			double[][] tokenProb = new double[numParticles][numTokens];
			for (int i = 0; i < numParticles; i++)
				sampleSingleParticleLR(msgs, tokenProb[i], mode);
			for (int i = 0; i < numTokens; i++) {
				double sum = 0.0;
				for (int j = 0; j < numParticles; j++)
					sum += tokenProb[j][i];
				llSeq += Math.log(sum) - Math.log(numParticles);
			}
			logLikelihood += llSeq;
		}
		return logLikelihood;
	}

	public static double estimateLogLikelihood(MessageSequenceCorpus<?, ?> corpus, double[][] phi,
			Object2DArray<double[][]> theta, Object2DArray<double[][]> pi, int numStates) {
		double[][] seqLogPi = new double[numStates][numStates];
		double[] tmpLL = new double[numStates];
		double[] prevLL = new double[numStates];
		double[] curLL = new double[numStates];

		double logLikelihood = LogSpace.ONE;
		Iterator<MessageSequenceCorpus.MessageSequence> seqIt = corpus.sequenceIterator();
		while (seqIt.hasNext()) {
			MessageSequenceCorpus.MessageSequence seq = seqIt.next();

			// forward algorithm with likelihood of message as "output probability"
			boolean isFirstMessage = true;
			double[][] seqTheta = null, seqPi = null;
			while (seq.hasNext()) {
				ProcessedMessage msg = seq.next();
				if (isFirstMessage) {
					seqTheta = theta.get(msg.getSenderId(), msg.getRecipientIds()[0]);
					seqPi = pi.get(msg.getSenderId(), msg.getRecipientIds()[0]);
					if (seq.hasNext()) {
						for (int i = 0; i < numStates; i++) {
							for (int j = 0; j < numStates; j++)
								seqLogPi[i][j] = LogSpace.INST.log(seqPi[i][j]);
						}
					}
				}

				for (int i = 0; i < numStates; i++) {
					if (isFirstMessage) {
						curLL[i] = LogSpace.INST.log(seqPi[0][i]);
					} else {
						for (int j = 0; j < numStates; j++)
							tmpLL[j] = LogSpace.product(prevLL[j], seqLogPi[j][i]);
						curLL[i] = LogSpace.INST.sum(tmpLL);
					}

					double msgLogLikelihood = LogSpace.ONE;
					for (int wordId : msg.getWordIds()) {
						double tokenLikelihood = 0.0;
						for (int j = 0; j < phi.length; j++)
							tokenLikelihood += seqTheta[i][j] * phi[j][wordId];
						msgLogLikelihood = LogSpace.product(msgLogLikelihood, LogSpace.INST.log(tokenLikelihood));
					}
					curLL[i] = LogSpace.product(curLL[i], msgLogLikelihood);
				}

				double[] tmp = prevLL;
				prevLL = curLL;
				curLL = tmp;
				isFirstMessage = false;
			}

			logLikelihood = LogSpace.product(logLikelihood, LogSpace.INST.sum(prevLL));
		}
		return logLikelihood;
	}

}
