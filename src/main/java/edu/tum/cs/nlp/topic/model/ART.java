package edu.tum.cs.nlp.topic.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import com.carrotsearch.hppc.LongObjectScatterMap;
import com.carrotsearch.hppc.cursors.LongObjectCursor;
import com.carrotsearch.hppc.cursors.ObjectCursor;

import edu.tum.cs.math.dist.DiscreteDistributionSampler;
import edu.tum.cs.nlp.corpus.Corpus;
import edu.tum.cs.nlp.corpus.ProcessedMessage;
import edu.tum.cs.util.arrays.Double2DArray;
import edu.tum.cs.util.arrays.Int2DArray;
import edu.tum.cs.util.arrays.Object2DArray;
import edu.tum.cs.util.arrays.Sparse2DIterator;
import edu.tum.cs.util.arrays.SparseObject2DArray;
import edu.tum.cs.util.io.Serializer;

public class ART extends LDABasedGibbsSamplingModel<ProcessedMessage> {

	private static final long serialVersionUID = 2277055559511116400L;
	private static final Logger logger = Logger.getLogger(ART.class.getName());

	private static final int numTopicBits = 8;	// max. 255 topics
	// private static final int numTopicBits = 10;	// max. 1024 topics

	protected static class PerAuthorRecipientData {
		public final int[] numTopic;
		public int sumNumTopic;
		public double aNormSquared;	// FastLDA

		public PerAuthorRecipientData(int numTopics) {
			numTopic = new int[numTopics];
		}

		public PerAuthorRecipientData(PerAuthorRecipientData other) {
			numTopic = other.numTopic.clone();
			sumNumTopic = other.sumNumTopic;
			aNormSquared = other.aNormSquared;
		}
	}

	public interface QueryResult {
		public Object2DArray<double[]> estimateTheta();
		public double[][] estimateAggregatedTheta(boolean aggregateAuthors, Set<Integer> authors,
				boolean aggregateRecipients, Set<Integer> recipients);
	}

	public class MeanQueryResult implements QueryResult {
		public final List<QueryResult> results = new ArrayList<QueryResult>();

		public void add(QueryResult result) {
			results.add(result);
		}

		@Override
		public Object2DArray<double[]> estimateTheta() {
			Object2DArray<double[]> meanThetas = new SparseObject2DArray<>(numPersons, numPersons);
			double z = 1.0 / results.size();
			for (QueryResult result : results) {
				Object2DArray<double[]> thetas = result.estimateTheta();
				Sparse2DIterator<double[]> it = thetas.sparseIterator();
				while (it.hasNext()) {
					it.advance();
					double[] curTheta = it.getValue();
					int row = it.getRow();
					int col = it.getColumn();

					double[] curMeanTheta = meanThetas.get(row, col);
					if (curMeanTheta == null) {
						curMeanTheta = new double[numTopics];
						meanThetas.set(row, col, curMeanTheta);
					}
					for (int i = 0; i < curMeanTheta.length; i++)
						curMeanTheta[i] += z * curTheta[i];
				}
			}
			return meanThetas;
		}

		@Override
		public double[][] estimateAggregatedTheta(boolean aggregateAuthors, Set<Integer> authors,
				boolean aggregateRecipients, Set<Integer> recipients) {
			double[][] meanThetas = null;
			double z = 1.0 / results.size();
			for (QueryResult result : results) {
				double[][] thetas = result.estimateAggregatedTheta(aggregateAuthors, authors, aggregateRecipients,
						recipients);
				if (meanThetas == null)
					meanThetas = new double[thetas.length][thetas[0].length];
				for (int i = 0; i < thetas.length; i++)
					for (int j = 0; j < thetas[i].length; j++)
						meanThetas[i][j] += z * thetas[i][j];
			}
			return meanThetas;
		}
	}

	public class SamplingState extends LDABasedGibbsSamplingModel<ProcessedMessage>.SamplingState
			implements QueryResult {
		private static final long serialVersionUID = -4124079102199592430L;

		/**
		 * Topic and recipient index are stored in the same array element. This is a micro-optimization: Storing the
		 * same data in two 2D arrays of type short would incur an additional top-level array of pointers. Flattening
		 * all multidimensional arrays to 1D arrays might be more efficient.
		 * Of course, it would be best not to store the assignments at all, and store the documents as (sparse) word
		 * histograms. Originally, it was planned to implement a model that relaxes the bag-of-words assumption,
		 * similar to LDA-HMM, which would have required a representation of documents as token sequences.
		 */
		protected int[][] xz;

		protected transient LongObjectScatterMap<PerAuthorRecipientData> perAuthorRecipientData;

		protected SamplingState(boolean isQuery, SamplingType type) {
			super(isQuery, type);
		}

		public SamplingState(int numDocuments, boolean isQuery, SamplingType type) {
			super(isQuery, type);
			xz = new int[numDocuments][];
			perAuthorRecipientData = new LongObjectScatterMap<PerAuthorRecipientData>(numPersons);
		}

		public QueryResult cloneResults() {
			SamplingState resultState = new SamplingState(true, type);
			resultState.perAuthorRecipientData = new LongObjectScatterMap<PerAuthorRecipientData>(numPersons);
			for (LongObjectCursor<PerAuthorRecipientData> csr : perAuthorRecipientData)
				resultState.perAuthorRecipientData.put(csr.key, new PerAuthorRecipientData(csr.value));
			return resultState;
		}

		private PerAuthorRecipientData initPerAuthorRecipientData(int author, int recipient) {
			long key = ((long) author << 32) | recipient;
			PerAuthorRecipientData ar = perAuthorRecipientData.get(key);
			if (ar == null) {
				ar = new PerAuthorRecipientData(numTopics);
				perAuthorRecipientData.put(key, ar);
			}
			return ar;
		}

		private PerAuthorRecipientData getPerAuthorRecipientData(int author, int recipient) {
			long key = ((long) author << 32) | recipient;
			return perAuthorRecipientData.get(key);
		}

		private Iterable<ObjectCursor<PerAuthorRecipientData>> getPerAuthorRecipientIterable() {
			return perAuthorRecipientData.values();
		}

		private void initANormSquared() {
			for (ObjectCursor<PerAuthorRecipientData> csr : getPerAuthorRecipientIterable()) {
				PerAuthorRecipientData ar = csr.value;
				ar.aNormSquared = 0.0;
				for (int k = 0; k < numTopics; k++) {
					double v = ar.numTopic[k] + alpha[k];
					ar.aNormSquared += v * v;
				}
			}
		}

		@Override
		public Object2DArray<double[]> estimateTheta() {
			Object2DArray<double[]> theta = new SparseObject2DArray<>(numPersons, numPersons);

			for (LongObjectCursor<PerAuthorRecipientData> csr : perAuthorRecipientData) {
				PerAuthorRecipientData ar = csr.value;
				double[] thetaAutorRecipient = new double[numTopics];
				for (int k = 0; k < numTopics; k++)
					thetaAutorRecipient[k] = (ar.numTopic[k] + alpha[k]) / (ar.sumNumTopic + sumAlpha);
				theta.set((int) (csr.key >>> 32), (int) csr.key, thetaAutorRecipient);
			}

			return theta;
		}

		@Override
		public double[][] estimateAggregatedTheta(boolean aggregateAuthors, Set<Integer> authors,
				boolean aggregateRecipients, Set<Integer> recipients) {
			double[][] theta;
			if (aggregateAuthors && aggregateRecipients)
				theta = new double[1][numTopics];
			else if (aggregateAuthors || aggregateRecipients)
				theta = new double[numPersons][numTopics];
			else
				throw new IllegalArgumentException("no aggregation specified");

			int[] sumAgg = new int[theta.length];
			for (LongObjectCursor<PerAuthorRecipientData> csr : perAuthorRecipientData) {
				PerAuthorRecipientData ar = csr.value;
				int author = (int) (csr.key >>> 32);
				int recipient = (int) csr.key;

				int aggKey = 0;
				if (aggregateAuthors) {
					if ((authors != null) && !authors.contains(author))
						continue;
					if (!aggregateRecipients)
						aggKey = recipient;
				}
				if (aggregateRecipients) {
					if ((recipients != null) && !recipients.contains(recipient))
						continue;
					if (!aggregateAuthors)
						aggKey = author;
				}

				for (int k = 0; k < numTopics; k++) {
					theta[aggKey][k] += ar.numTopic[k];
					sumAgg[aggKey] += ar.numTopic[k];
				}
			}

			for (int i = 0; i < theta.length; i++)
				for (int j = 0; j < numTopics; j++)
					theta[i][j] = (theta[i][j] + alpha[j]) / (sumAgg[i] + sumAlpha);

			return theta;
		}

		@Override
		public int[][] getTopicAssignment() {
			int[][] z = new int[xz.length][];
			for (int i = 0; i < xz.length; i++) {
				int n = xz[i].length;
				z[i] = new int[n];
				for (int j = 0; j < n; j++)
					z[i][j] = xz[i][j] & ((1 << numTopicBits) - 1);
			}
			return z;
		}

		public int[][] getRecipientAssignment() {
			int[][] x = new int[xz.length][];
			for (int i = 0; i < xz.length; i++) {
				int n = xz[i].length;
				x[i] = new int[n];
				for (int j = 0; j < n; j++)
					x[i][j] = xz[i][j] >>> numTopicBits;
			}
			return x;
		}
	}

	// values taken from PhD thesis of H. Wallach (p.31); min. value is 10^-16 for non-uniform, 10^-13 for uniform case
	private static final int alphaMaxIter = 100;
	private static final double alphaEpsilon = 1.0E-13;

	// parameters
	protected int numPersons;
	public static enum AlphaType { CONSTANT, DATA_DRIVEN_UNIFORM, DATA_DRIVEN_NON_UNIFORM };
	protected AlphaType alphaSamplingType = AlphaType.CONSTANT;

	// model state after training (for querying)
	protected SamplingState modelState;
	protected double[][] phi;

	public ART(int numWords, int numPersons, int numTopics, double alpha, double beta) {
		super(numWords, numTopics, alpha, beta);
		this.numPersons = numPersons;
	}

	public ART(int numWords, int numPersons, int numTopics, double[] alpha, double[] beta) {
		super(numWords, numTopics, alpha, beta);
		this.numPersons = numPersons;
	}

	public ART(int numWords, int numPersons, int numTopics, double[] alpha, Double2DArray beta) {
		super(numWords, numTopics, alpha, beta);
		this.numPersons = numPersons;
	}

	protected ART(ART other) {
		super(other);
		numPersons = other.numPersons;
		alphaSamplingType = other.alphaSamplingType;
		if (alphaSamplingType != AlphaType.CONSTANT)
			alpha = alpha.clone();
	}

	@Override
	public ART createSimilar() {
		return new ART(this);
	}

	public void setAlphaSamplingType(AlphaType type) {
		alphaSamplingType = type;
	}

	@Override
	public SamplingState train(Corpus<ProcessedMessage> corpus) {
		modelState = new SamplingState(corpus.size(), false, samplingType);
		initSamplingState(sampler, corpus, modelState, true);
		burnIn(sampler, corpus, modelState, burnIn);
		estimateParameters(modelState);
		return modelState;
	}

	/** Call on a model obtained via {@link Serializer#loadObjectFromFile(String)} to continue training. */
	@Override
	public SamplingState resumeTraining(Corpus<ProcessedMessage> corpus, int numRemainingIter) {
		SamplingState assignmentState = modelState;
		modelState = new SamplingState(corpus.size(), false, samplingType);
		modelState.xz = assignmentState.xz;
		initSamplingState(sampler, corpus, modelState, false);

		burnIn(sampler, corpus, modelState, numRemainingIter);
		estimateParameters(modelState);
		return modelState;
	}

	private void estimateParameters(SamplingState state) {
		phi = state.estimatePhi();
		state.compact(keepNumWordTopic);
	}

	protected void burnIn(DiscreteDistributionSampler sampler, Corpus<ProcessedMessage> corpus, SamplingState state,
			int numIter) {
		// sample until required number of iterations for burn in is reached
		logger.info("Starting sampling with type " + samplingType);
		long t = System.currentTimeMillis();
		int reportIdx = 0;
		int saveIdx = 0;
		for (int i = 0; i < numIter; i++) {
			if (!state.isQuery)
				sampleHyperParameters(state, alphaSamplingType);

			sample(sampler, corpus, state);

			if (++reportIdx == burnInReportRate) {
				logger.info("Burn in step " + (i + 1) + " of " + numIter);
				reportIdx = 0;
			}
			if (!state.isQuery && (++saveIdx == burnInSaveRate)) {
				logger.info("Saving ART model");
				Serializer.saveObjectToFile(this, intermediateStateFile);
				logger.info("Done saving ART model");
				saveIdx = 0;
			}
		}

		if (numIter > 0) {
			t = System.currentTimeMillis() - t;
			float ips = numIter / ((float) t / 1000);
			logger.info(ips + " iterations per second");
		}
	}

	protected void initSamplingState(DiscreteDistributionSampler sampler, Corpus<ProcessedMessage> corpus,
			final SamplingState state, boolean assignRandom) {
		for (ProcessedMessage message : corpus) {
			int posMessage = message.getIndex();
			int author = message.getSenderId();
			int[] recipients = message.getRecipientIds();

			PerAuthorRecipientData[] ar = new PerAuthorRecipientData[recipients.length];
			for (int i = 0; i < ar.length; i++)
				ar[i] = state.initPerAuthorRecipientData(author, recipients[i]);

			if (assignRandom)
				state.xz[posMessage] = new int[message.size()];
			else if (state.xz[posMessage].length != message.size())
				throw new RuntimeException("length mismatch for message " + posMessage +
						": message size = " + message.size() + ", saved state = " + state.xz[posMessage].length);

			int posWord = 0;
			for (int wordIndex : message.getWordIds()) {
				int recipientIdx;
				int topic;
				if (assignRandom) {
					// assign recipients and topics randomly
					recipientIdx = sampler.sample1DUniform(recipients.length);
					topic = sampler.sample1DUniform(numTopics);
					state.xz[posMessage][posWord] = (recipients[recipientIdx] << numTopicBits) | topic;
				} else {
					int recipient = state.xz[posMessage][posWord] >>> numTopicBits;
					for (recipientIdx = 0; recipientIdx < recipients.length; recipientIdx++)
						if (recipients[recipientIdx] == recipient)
							break;
					if (recipientIdx == recipients.length) {
						logger.warning("invalid recipient " + recipient);
						recipientIdx = 0;
					}
					topic = state.xz[posMessage][posWord] & ((1 << numTopicBits) - 1);
				}

				ar[recipientIdx].numTopic[topic]++;
				ar[recipientIdx].sumNumTopic++;

				if (!state.isQuery) {
					state.numWordTopic.adjust(wordIndex, topic, 1);
					state.sumNumTopic[topic]++;
				}

				posWord++;
			}
		}

		// initialize data structures for FastLDA
		if (state.type != SamplingType.REGULAR) {
			if (!state.isQuery)
				state.initFastLDA();
			state.initANormSquared();
		}
	}

	private void sampleHyperParameters(SamplingState state, AlphaType alphaType) {
		if (alphaType != AlphaType.CONSTANT) {
			// Wallach suggests averaging the update over multiple MCMC samples, but according to Andrieu et al.,
			// "An Introduction to MCMC for Machine Learning" (p.24), using a single sample is valid as well.
			if (alphaType == AlphaType.DATA_DRIVEN_UNIFORM)
				sampleAlphaUniform(state);
			else if (alphaType == AlphaType.DATA_DRIVEN_NON_UNIFORM)
				sampleAlphaNonUniform(state);

			if (state.type == SamplingType.FASTLDA)
				state.initANormSquared();	// recompute aNormSquared
		}
	}

	private void sampleAlphaUniform(SamplingState state) {
		double curAlpha = alpha[0];	// assume uniform alpha
		double prevAlpha;
		int numIter = 0;
		do {
			prevAlpha = curAlpha;

			double num = 0.0;
			double denom = 0.0;
			for (ObjectCursor<PerAuthorRecipientData> csr : state.getPerAuthorRecipientIterable()) {
				int sumN = 0;
				for (int n : csr.value.numTopic) {
					for (int i = 0; i < (n - 1); i++)
						num += 1.0 / (curAlpha + i);
					sumN += n;
				}

				for (int i = 0; i < (sumN - 1); i++)
					denom += 1.0 / ((curAlpha * numTopics) + i);
			}
			num /= numTopics;
			curAlpha *= num / denom;
		} while ((++numIter < alphaMaxIter) && (Math.abs(curAlpha - prevAlpha) > alphaEpsilon));
		setUniformAlpha(curAlpha);	// updates sumAlpha
	}

	private void sampleAlphaNonUniform(SamplingState state) {
		List<Map<Integer, Integer>> topicBins = new ArrayList<Map<Integer, Integer>>(numTopics);
		for (int i = 0; i < numTopics; i++)
			topicBins.add(new HashMap<Integer, Integer>());
		Map<Integer, Integer> sumBin = new HashMap<Integer, Integer>();
		int[] maxN = new int[numTopics + 1];

		for (ObjectCursor<PerAuthorRecipientData> csr : state.getPerAuthorRecipientIterable()) {
			int sum = 0;
			for (int i = 0; i < numTopics; i++) {
				int n = csr.value.numTopic[i];
				sum += n;
				if (n > 1) {
					Map<Integer, Integer> topicBin = topicBins.get(i);
					Integer v = topicBin.get(n);
					if (v == null)
						topicBin.put(n, 1);
					else
						topicBin.put(n, v + 1);
					if (n > maxN[i])
						maxN[i] = n;
				}
			}

			if (sum > 1) {
				Integer v = sumBin.get(sum);
				if (v == null)
					sumBin.put(sum, 1);
				else
					sumBin.put(sum, v + 1);
				if (sum > maxN[numTopics])
					maxN[numTopics] = sum;
			}
		}

		double maxAlphaDiff;
		int numIter = 0;
		do {
			double denom = 0.0;
			double d = 0.0;
			for (int i = 2; i <= maxN[numTopics]; i++) {
				d += 1.0 / (sumAlpha + (i - 2));
				Integer n = sumBin.get(i);
				if (n != null)
					denom += n * d;
			}

			maxAlphaDiff = 0.0;
			for (int i = 0; i < numTopics; i++) {
				Map<Integer, Integer> topicBin = topicBins.get(i);
				double num = 0.0;
				d = 0.0;
				for (int j = 2; j <= maxN[i]; j++) {
					d += 1.0 / (alpha[i] + (j - 2));
					Integer n = topicBin.get(j);
					if (n != null)
						num += n * d;
				}

				double newAlpha;
				if (num == 0.0)
					newAlpha = alphaEpsilon;	// alpha must not be 0, so set to smallest possible value
				else
					newAlpha = alpha[i] * (num / denom);
				double alphaDiff = Math.abs(alpha[i] - newAlpha);
				if (alphaDiff > maxAlphaDiff)
					maxAlphaDiff = alphaDiff;
				alpha[i] = newAlpha;
			}

			sumAlpha = 0.0;
			for (double a : alpha)
				sumAlpha += a;
		} while ((++numIter < alphaMaxIter) && (maxAlphaDiff > alphaEpsilon));
	}

	protected void sample(DiscreteDistributionSampler sampler, Corpus<ProcessedMessage> corpus, SamplingState state) {
		for (ProcessedMessage message : corpus)
			sampleMessage(sampler, message, state);
	}

	private void sampleMessage(DiscreteDistributionSampler sampler, ProcessedMessage message, SamplingState state) {
		switch (state.type) {
		case FASTLDA:
			sampleMessageFastLda(sampler, message, state);
			break;
		case REGULAR:
			sampleMessageRegular(sampler, message, state);
			break;
		}
	}

	private int sampleRecipientTopicAssignment(int wordIndex, PerAuthorRecipientData[] ar, Int2DArray numWordTopic,
			int[] sumNumTopic, double[] cumulRecipientTopicProbs, boolean isQuery) {
		double prevCumulProb = 0.0;
		int idx = 0;
		for (int i = 0; i < ar.length; i++) {
			for (int j = 0; j < numTopics; j++) {
				if (!isQuery) {
					cumulRecipientTopicProbs[idx] = prevCumulProb +
							(((ar[i].numTopic[j] + alpha[j]) / (ar[i].sumNumTopic + sumAlpha)) *
							((numWordTopic.get(wordIndex, j) + beta.get(j, wordIndex)) /
							 (sumNumTopic[j] + sumBeta[j])));
				} else {
					cumulRecipientTopicProbs[idx] = prevCumulProb +
							(((ar[i].numTopic[j] + alpha[j]) / (ar[i].sumNumTopic + sumAlpha)) * phi[j][wordIndex]);
				}
				prevCumulProb = cumulRecipientTopicProbs[idx];
				idx++;
			}
		}
		return sampler.sample1D(cumulRecipientTopicProbs, idx);
	}

	private int sampleRecipientTopicAssignmentQuery(int wordIndex, PerAuthorRecipientData[] ar,
			double[] cumulRecipientTopicProbs) {
		return sampleRecipientTopicAssignment(wordIndex, ar, null, null, cumulRecipientTopicProbs, true);
	}

	private void sampleMessageRegular(DiscreteDistributionSampler sampler, ProcessedMessage message,
			SamplingState state) {
		int author = message.getSenderId();
		int[] recipients = message.getRecipientIds();
		PerAuthorRecipientData[] ar = new PerAuthorRecipientData[recipients.length];
		for (int i = 0; i < ar.length; i++)
			ar[i] = state.getPerAuthorRecipientData(author, recipients[i]);

		double[] cumulRecipientTopicProbs = new double[recipients.length * numTopics];
		int posMessage = message.getIndex();
		int posWord = 0;
		for (int wordIndex : message.getWordIds()) {
			// do not count currently assigned recipient and topic
			int recipient = state.xz[posMessage][posWord] >>> numTopicBits;
			int topic = state.xz[posMessage][posWord] & ((1 << numTopicBits) - 1);

			// TODO: should store index into per-message recipient list in xz
			PerAuthorRecipientData arCur = state.getPerAuthorRecipientData(author, recipient);
			arCur.numTopic[topic]--;
			arCur.sumNumTopic--;
			if (!state.isQuery) {
				state.numWordTopic.adjust(wordIndex, topic, -1);
				state.sumNumTopic[topic]--;
			}

			// compute joint probability of new recipient and topic and sample
			int v = sampleRecipientTopicAssignment(wordIndex, ar, state.numWordTopic, state.sumNumTopic,
					cumulRecipientTopicProbs, state.isQuery);
			recipient = recipients[v / numTopics];
			topic = v % numTopics;

			// assign new recipient and topic and update counts
			state.xz[posMessage][posWord] = (recipient << numTopicBits) | topic;
			arCur = state.getPerAuthorRecipientData(author, recipient);
			arCur.numTopic[topic]++;
			arCur.sumNumTopic++;
			if (!state.isQuery) {
				state.numWordTopic.adjust(wordIndex, topic, 1);
				state.sumNumTopic[topic]++;
			}

			posWord++;
		}
	}

	private void sampleMessageFastLda(DiscreteDistributionSampler sampler, ProcessedMessage message,
			SamplingState state) {
		int author = message.getSenderId();
		int[] recipients = message.getRecipientIds();

		double[] cumulTopicProbs = new double[numTopics];
		double[] cumulRecipientProbs = new double[recipients.length];
		int posMessage = message.getIndex();
		int posWord = 0;
		for (int wordIndex : message.getWordIds()) {
			int recipient = state.xz[posMessage][posWord] >>> numTopicBits;
			int topic = state.xz[posMessage][posWord] & ((1 << numTopicBits) - 1);
			PerAuthorRecipientData ar = state.getPerAuthorRecipientData(author, recipient);

			/*
			// test norm invariants; really expensive, only for debugging!
			double actualANorm = 0.0;
			double actualBNorm = 0.0;
			for (int k = 0; k < numTopics; k++) {
				int actualK = state.sortedWordTopic.get(wordIndex, k);
				double v = ar.numTopic[actualK] + alpha[actualK];
				actualANorm += v * v;
			}
			for (int k = 0; k < numTopics; k++) {
				int n = 0;
				for (int a = 0; a < numTopics; a++) {
					if (state.sortedWordTopic.get(wordIndex, a) == k)
						n++;
					if (n > 1)
						throw new RuntimeException("sort order");
				}
				int actualK = state.sortedWordTopic.get(wordIndex, k);
				double v = state.numWordTopic.get(wordIndex, actualK) + beta.get(actualK, wordIndex);
				actualBNorm += v * v;
			}
			double cachedANorm = ar.aNormSquared;
			double cachedBNorm = state.bNormSquared[wordIndex];
			if (Math.abs(cachedANorm - actualANorm) > 0.1)
				System.err.println("error in a norm, expected " + actualANorm + ", got " + cachedANorm);
			if (Math.abs(cachedBNorm - actualBNorm) > 0.1)
				System.err.println("error in b norm, expected " + actualBNorm + ", got " + cachedBNorm);
			*/


			// sample recipient "x" by simple 1D binary search
			// -----------------------------------------------

			// do not count currently assigned recipient
			ar.numTopic[topic]--;
			ar.sumNumTopic--;
			ar.aNormSquared -= 2.0 * (ar.numTopic[topic] + alpha[topic] + 0.5);

			if (recipients.length > 1) {
				// compute probability distribution of recipients
				double prevCumulProb = 0.0;
				for (int i = 0; i < recipients.length; i++) {
					ar = state.getPerAuthorRecipientData(author, recipients[i]);
					cumulRecipientProbs[i] = prevCumulProb +
							((ar.numTopic[topic] + alpha[topic]) / (ar.sumNumTopic + sumAlpha));
					prevCumulProb = cumulRecipientProbs[i];
				}

				// sample new recipient - result is cached, xz is updated later
				int recipientIdx = sampler.sample1D(cumulRecipientProbs);
				recipient = recipients[recipientIdx];

				// from here on, recipient is fixed
				ar = state.getPerAuthorRecipientData(author, recipient);
			}


			// sample topic "z" by FastLDA algorithm
			// -------------------------------------

			SamplingState topicState = state;
			if (!state.isQuery) {
				// do not count currently assigned topic
				int n = state.numWordTopic.adjust(wordIndex, topic, -1);
				state.sumNumTopic[topic]--;
				state.bNormSquared[wordIndex] -= 2.0 * (n + beta.get(topic, wordIndex) + 0.5);
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
			double aNormSquared = ar.aNormSquared;
			double bNormSquared = topicState.bNormSquared[wordIndex];
	sampleTopic:
			for (int i = 0; i < numTopics; i++) {
				int t = topicState.sortedWordTopic.get(wordIndex, i);

				// compute cumulative probabilities up to current topic
				double messageTerm = ar.numTopic[t] + alpha[t];
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

			// got new recipient/topic, update counts
			ar.numTopic[topic]++;
			ar.sumNumTopic++;
			ar.aNormSquared += 2.0 * (ar.numTopic[topic] + alpha[topic] - 0.5);
			if (!state.isQuery) {
				int n = state.numWordTopic.adjust(wordIndex, topic, 1);
				state.sumNumTopic[topic]++;
				state.bNormSquared[wordIndex] += 2.0 * (n + beta.get(topic, wordIndex) - 0.5);
				state.updateSortInc(topic, wordIndex);
				state.trackWordTopicChange(wordIndex, topic);
			}

			state.xz[posMessage][posWord] = (recipient << numTopicBits) | topic;

			posWord++;
		}
	}

	@Override
	public double[][] getTopicWordDistr() {
		return phi;
	}

	public QueryResult query(Corpus<ProcessedMessage> corpus) {
		SamplingState state = new SamplingState(corpus.size(), true, samplingType);
		initSamplingState(sampler, corpus, state, true);
		burnIn(sampler, corpus, state, burnIn);

		if (querySamples == 1)
			return state;

		MeanQueryResult result = new MeanQueryResult();
		result.add(state);
		for (int i = 1; i < querySamples; i++) {
			for (int j = 0; j < queryLag + 1; j++)
				sample(sampler, corpus, state);
			result.add(state.cloneResults());
		}
		return result;
	}

	private PerAuthorRecipientData[] collectAuthorRecipientData(ProcessedMessage msg,
			Map<Integer, PerAuthorRecipientData> recipientData) {
		PerAuthorRecipientData[] ar = new PerAuthorRecipientData[msg.getRecipientIds().length];
		int idx = 0;
		for (int recipient : msg.getRecipientIds())
			ar[idx++] = recipientData.get(recipient);
		return ar;
	}

	private void sampleWordLR(int pos, ProcessedMessage msg, int wordIndex, int[] xz,
			Map<Integer, PerAuthorRecipientData> recipientData, PerAuthorRecipientData[] ar, double[] cumulTopicProbs,
			boolean isNewWord) {
		int recipient = xz[pos] >>> numTopicBits;
		PerAuthorRecipientData arCur = recipientData.get(recipient);

		int topic;
		if (!isNewWord) {
			topic = xz[pos] & ((1 << numTopicBits) - 1);
			arCur.numTopic[topic]--;
			arCur.sumNumTopic--;
		}

		int v = sampleRecipientTopicAssignmentQuery(wordIndex, ar, cumulTopicProbs);
		recipient = msg.getRecipientIds()[v / numTopics];
		topic = v % numTopics;
		xz[pos] = (recipient << numTopicBits) | topic;

		arCur = recipientData.get(recipient);
		arCur.numTopic[topic]++;
		arCur.sumNumTopic++;
	}

	private void sampleSingleParticleLR(ArrayList<ProcessedMessage> messages, Set<Integer> recipients,
			double[] tokenProb, ResampleMode mode) {
		int n = tokenProb.length;
		int[] xz = new int[n];
		Map<Integer, PerAuthorRecipientData> recipientData = new HashMap<Integer, PerAuthorRecipientData>();
		for (Integer recipient : recipients)
			recipientData.put(recipient, new PerAuthorRecipientData(numTopics));
		double[] cumulTopicProbs = new double[recipients.size() * numTopics];

		double[] cumulMessageProb = null;
		int[] msgStartPos = null;
		if (mode == ResampleMode.LOG) {
			cumulMessageProb = new double[messages.size()];
			msgStartPos = new int[messages.size()];
			double prevCumulProb = 0.0;
			int idx = 0;
			for (ProcessedMessage msg : messages) {
				if ((idx + 1) < msgStartPos.length)
					msgStartPos[idx + 1] = msgStartPos[idx] + msg.size();
				cumulMessageProb[idx] = prevCumulProb + msg.size();
				prevCumulProb = cumulMessageProb[idx];
				idx++;
			}
		}

		int pos = 0;
		int posMessage = 0;
		for (ProcessedMessage msg : messages) {
			int[] y = msg.getWordIds();
			PerAuthorRecipientData[] ar = collectAuthorRecipientData(msg, recipientData);

			for (int i = 0; i < y.length; i++) {
				if (pos > 0) {
					// if resampling mode is not "NONE", resample topic assignments of previous words...
					if (mode == ResampleMode.LOG) {
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
							ProcessedMessage msgPrev = messages.get(posMessageRes);
							int[] yPrev = msgPrev.getWordIds();
							int maxIdx = (posMessageRes == posMessage) ? i : yPrev.length;
							PerAuthorRecipientData[] arRes = collectAuthorRecipientData(msgPrev, recipientData);

							for (int j = 0; j < e.getValue(); j++) {
								int idx = sampler.prng.nextInt(maxIdx);
								sampleWordLR(msgStartPos[posMessageRes] + idx, msgPrev, yPrev[idx], xz, recipientData,
										arRes, cumulTopicProbs, false);
							}
						}
					} else if (mode == ResampleMode.ALL) {
						int posPrev = 0;
resamplingLoop:
						for (ProcessedMessage msgPrev : messages) {
							int[] yPrev = msgPrev.getWordIds();
							PerAuthorRecipientData[] arRes = collectAuthorRecipientData(msgPrev, recipientData);

							for (int j = 0; j < yPrev.length; j++) {
								sampleWordLR(posPrev, msgPrev, yPrev[j], xz, recipientData, arRes, cumulTopicProbs,
										false);

								if (++posPrev == pos)
									break resamplingLoop;
							}
						}
					}	// ...otherwise, skip resampling
				}

				// compute predictive probability of current word, given the topic assignments of previous words
				double p = 0.0;
				for (int j = 0; j < msg.getRecipientIds().length; j++) {
					for (int k = 0; k < numTopics; k++) {
						PerAuthorRecipientData arCur = recipientData.get(msg.getRecipientIds()[j]);
						p += (arCur.numTopic[k] + alpha[k]) / (arCur.sumNumTopic + sumAlpha) * phi[k][y[i]];
					}
				}
				tokenProb[pos] += p / msg.getRecipientIds().length;

				// sample topic assignment of current word
				sampleWordLR(pos, msg, y[i], xz, recipientData, ar, cumulTopicProbs, true);

				pos++;
			}
			posMessage++;
		}
	}

	private double sampleParticlesSenderLR(List<ProcessedMessage> senderMessages, ResampleMode mode) {
		int numParticles = likelihoodSamples;
		double logLikelihood = 0.0;
		while (!senderMessages.isEmpty()) {
			// TODO: This partitioning should be moved into class MessageCorpus so that ParallelART can also benefit.
			// Find a minimal set of messages that share no recipients with messages outside of the set.
			ArrayList<ProcessedMessage> curMessages = new ArrayList<ProcessedMessage>();
			Set<Integer> curRecipients = new HashSet<Integer>();
			int numTokens = 0;
			boolean setChanged;
			do {
				setChanged = false;
				Iterator<ProcessedMessage> it = senderMessages.iterator();
				while (it.hasNext()) {
					ProcessedMessage msg = it.next();
					boolean hasOverlappingRecipients = false;
					if (curRecipients.isEmpty()) {
						hasOverlappingRecipients = true;
					} else {
						for (int recipient : msg.getRecipientIds()) {
							if (curRecipients.contains(recipient)) {
								hasOverlappingRecipients = true;
								break;
							}
						}
					}

					if (hasOverlappingRecipients) {
						curMessages.add(msg);
						it.remove();
						for (int recipient : msg.getRecipientIds())
							setChanged |= curRecipients.add(recipient);
						numTokens += msg.size();
					}
				}
			} while (setChanged);

			double llSet = 0.0;
			double[][] tokenProb = new double[numParticles][numTokens];
			for (int i = 0; i < numParticles; i++)
				sampleSingleParticleLR(curMessages, curRecipients, tokenProb[i], mode);
			for (int i = 0; i < numTokens; i++) {
				double sum = 0.0;
				for (int j = 0; j < numParticles; j++)
					sum += tokenProb[j][i];
				llSet += Math.log(sum) - Math.log(numParticles);
			}
			logLikelihood += llSet;
		}
		return logLikelihood;
	}

	@Override
	public double queryLogLikelihoodLR(Corpus<ProcessedMessage> corpus, ResampleMode mode) {
		double logLikelihood = 0.0;
		// Assume that the messages in the corpus are sorted by sender. Gather all messages of one sender.
		List<ProcessedMessage> senderMessages = new LinkedList<ProcessedMessage>();
		int curSender = -1;
		for (ProcessedMessage message : corpus) {
			if ((curSender < 0) || (curSender == message.getSenderId())) {
				curSender = message.getSenderId();
				senderMessages.add(message);
			} else {
				logLikelihood += sampleParticlesSenderLR(senderMessages, mode);

				senderMessages.clear();
				senderMessages.add(message);
				curSender = message.getSenderId();
			}
		}
		logLikelihood += sampleParticlesSenderLR(senderMessages, mode);
		return logLikelihood;
	}

	public static double estimateLogLikelihood(Corpus<ProcessedMessage> corpus, double[][] phi,
			Object2DArray<double[]> theta) {
		double logLikelihood = 0.0;

		for (ProcessedMessage message : corpus) {
			int author = message.getSenderId();
			int[] recipients = message.getRecipientIds();

			for (int wordIndex : message.getWordIds()) {
				double sumTopicRecipientProb = 0.0;
				for (int recipient : recipients) {
					double[] authorRecipientTheta = theta.get(author, recipient);
					for (int i = 0; i < phi.length; i++)	// iterate over topics
						sumTopicRecipientProb += authorRecipientTheta[i] * phi[i][wordIndex];
				}
				logLikelihood += Math.log(sumTopicRecipientProb);
			}

			logLikelihood -= message.size() * Math.log(recipients.length);
		}

		return logLikelihood;
	}

}
