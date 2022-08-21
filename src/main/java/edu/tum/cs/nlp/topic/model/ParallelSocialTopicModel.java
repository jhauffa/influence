package edu.tum.cs.nlp.topic.model;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Logger;

import com.carrotsearch.hppc.IntScatterSet;
import com.carrotsearch.hppc.IntSet;
import com.carrotsearch.hppc.cursors.IntCursor;

import edu.tum.cs.math.dist.DiscreteDistributionSampler;
import edu.tum.cs.nlp.corpus.Corpus;
import edu.tum.cs.nlp.corpus.ProcessedDocumentReader;
import edu.tum.cs.nlp.corpus.ProcessedMessage;
import edu.tum.cs.nlp.corpus.SlicedCorpus;
import edu.tum.cs.util.ExperimentConfiguration;
import edu.tum.cs.util.arrays.Double2DArray;
import edu.tum.cs.util.arrays.SparseShadowInt2DArray;

public class ParallelSocialTopicModel extends SocialTopicModel {

	private static final long serialVersionUID = -410008187562179424L;
	private static final Logger logger = Logger.getLogger(ParallelSocialTopicModel.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(ParallelSocialTopicModel.class);

	protected class SharedSamplingState extends SamplingState implements Runnable {
		private static final long serialVersionUID = -4494604289627780842L;
		private static final int numTopicBits = 8;	// max. 255 topics

		private final Corpus<ProcessedMessage> slice;
		private final DiscreteDistributionSampler sampler;

		public SharedSamplingState(SamplingState other, Corpus<ProcessedMessage> slice,
				DiscreteDistributionSampler sampler) {
			super(other.isQuery, other.type);
			this.slice = slice;
			this.sampler = sampler;

			// share arrays if elements are updated disjointly by the worker threads
			z = other.z;
			c = other.c;
			perAuthorRecipientData = other.perAuthorRecipientData;

			// deep clone all other arrays
			if (!other.isQuery) {
				numWordTopic = new SparseShadowInt2DArray(other.numWordTopic);
				sumNumTopic = other.sumNumTopic.clone();

				if (type == SamplingType.FASTLDA) {
					sortedWordTopic = new SparseShadowInt2DArray(other.sortedWordTopic);
					invSortedWordTopic = new SparseShadowInt2DArray(other.invSortedWordTopic);
					sortedTopic = other.sortedTopic.clone();
					invSortedTopic = other.invSortedTopic.clone();
					bNormSquared = other.bNormSquared.clone();
				}
			}
		}

		@Override
		public void trackWordTopicChange(int wordIndex, int topic) {
			synchronized (changesWordTopic) {
				changesWordTopic.add((wordIndex << numTopicBits) | topic);
			}
		}

		@Override
		public void run() {
			ParallelSocialTopicModel.super.sample(sampler, slice, this);
		}
	}

	protected int numThreads;
	private boolean useExternalMemoryCorpus;
	private boolean keepOriginalCorpus;

	private transient ExecutorService pool;
	private transient DiscreteDistributionSampler[] samplers;
	private transient Future<?>[] futures;
	private transient List<SharedSamplingState> localStates;
	private transient IntSet changesWordTopic;

	public ParallelSocialTopicModel(int numWords, int numPersons, int numTopics, int numStates, double alpha,
			double beta, double gamma) {
		super(numWords, numPersons, numTopics, numStates, alpha, beta, gamma);
		init();
	}

	public ParallelSocialTopicModel(int numWords, int numPersons, int numTopics, int numStates, double alpha,
			double beta, double gamma[][]) {
		super(numWords, numPersons, numTopics, numStates, alpha, beta, gamma);
		init();
	}

	public ParallelSocialTopicModel(int numWords, int numPersons, int numTopics, int numStates, double alpha[],
			double beta[], double gamma[][]) {
		super(numWords, numPersons, numTopics, numStates, alpha, beta, gamma);
		init();
	}

	public ParallelSocialTopicModel(int numWords, int numPersons, int numTopics, int numStates, double alpha[],
			Double2DArray beta, double gamma[][]) {
		super(numWords, numPersons, numTopics, numStates, alpha, beta, gamma);
		init();
	}

	protected ParallelSocialTopicModel(ParallelSocialTopicModel other) {
		super(other);
		useExternalMemoryCorpus = other.useExternalMemoryCorpus;
		init();
	}

	@Override
	public ParallelSocialTopicModel createSimilar() {
		return new ParallelSocialTopicModel(this);
	}

	private void readObject(ObjectInputStream stream) throws IOException, ClassNotFoundException {
		stream.defaultReadObject();
		init();
	}

	private void init() {
		assert (numTopics < (1 << SharedSamplingState.numTopicBits));
		assert (numWords < (1 << (32 - SharedSamplingState.numTopicBits)));

		changesWordTopic = new IntScatterSet();

		int numThreads = cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_THREADS,
				Runtime.getRuntime().availableProcessors());
		logger.info("using " + numThreads + " threads");
		setNumThreads(numThreads);
	}

	public void setUseExternalMemoryCorpus(boolean useExternalMemoryCorpus) {
		this.useExternalMemoryCorpus = useExternalMemoryCorpus;
	}

	public void setKeepOriginalCorpus(boolean keepOriginalCorpus) {
		this.keepOriginalCorpus = keepOriginalCorpus;
	}

	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
		localStates = new ArrayList<SharedSamplingState>(numThreads);
		futures = new Future<?>[numThreads];

		// DiscreteDistributionSampler is not thread safe due to keeping its own PRNG state, so we need one instance per
		// worker.
		samplers = new DiscreteDistributionSampler[numThreads];
		for (int i = 0; i < numThreads; i++)
			samplers[i] = new DiscreteDistributionSampler(sampler);
	}

	@Override
	protected void burnIn(DiscreteDistributionSampler sampler, Corpus<ProcessedMessage> corpus, SamplingState state,
			int numIter) {
		pool = Executors.newFixedThreadPool(numThreads);
		ProcessedDocumentReader<ProcessedMessage> reader = new ProcessedMessage.Reader();
		SlicedCorpus<ProcessedMessage> slices = new SlicedCorpus<ProcessedMessage>(corpus, reader, numThreads,
				useExternalMemoryCorpus, keepOriginalCorpus);
		super.burnIn(sampler, slices, state, numIter);
		pool.shutdown();
	}

	@Override
	protected void sample(DiscreteDistributionSampler sampler, Corpus<ProcessedMessage> corpus, SamplingState state) {
		if (pool.isShutdown()) {	// calling sample after burn-in
			super.sample(sampler, corpus, state);
			return;
		}

		List<Corpus<ProcessedMessage>> slices = corpus.slice(numThreads, keepOriginalCorpus);
		for (int i = 0; i < numThreads; i++) {
			SharedSamplingState localState = new SharedSamplingState(state, slices.get(i), samplers[i]);
			localStates.add(localState);
			futures[i] = pool.submit(localState);
		}

		for (Future<?> future : futures) {
			try {
				future.get();
			} catch (Exception ex) {
				throw new RuntimeException(ex);
			}
		}

		// update numTopicWord and sumNumTopic, sort order and bNormSquared
		if (!state.isQuery) {
			for (IntCursor csr : changesWordTopic) {
				int wordIndex = csr.value >>> SharedSamplingState.numTopicBits;
				int topic = csr.value & ((1 << SharedSamplingState.numTopicBits) - 1);

				int prevNum = state.numWordTopic.get(wordIndex, topic);
				int sumDiff = 0;
				for (SharedSamplingState localState : localStates)
					sumDiff += localState.numWordTopic.get(wordIndex, topic) - prevNum;
				if (sumDiff == 0)
					continue;

				state.numWordTopic.adjust(wordIndex, topic, sumDiff);
				state.sumNumTopic[topic] += sumDiff;

				if (state.type == SamplingType.FASTLDA) {
					if (sumDiff > 0)
						state.updateSortInc(topic, wordIndex);
					else
						state.updateSortDec(topic, wordIndex);

					double old = prevNum + beta.get(topic, wordIndex);
					state.bNormSquared[wordIndex] += ((old + sumDiff) * (old + sumDiff)) - (old * old);
				}
			}

			changesWordTopic.clear();
		}

		localStates.clear();
	}

	// TODO: In theory, queryLogLikelihoodLR can be trivially parallelized, like it is done in ParallelLDA and
	//	ParallelART. In practice, slicing a MessageSequenceCorpus invalidates ProcessedMessage.prev/nextIndex, which the
	//	sequence iterator uses for navigation.

}
