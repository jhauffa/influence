package edu.tum.cs.nlp.topic.model;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
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

public class ParallelART extends ART {

	private static final long serialVersionUID = -759416728172529228L;
	private static final Logger logger = Logger.getLogger(ParallelART.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(ParallelART.class);

	protected class SharedSamplingState extends SamplingState implements Runnable {
		private static final long serialVersionUID = -214695780592998006L;
		private static final int numTopicBits = 8;	// max. 255 topics
		// private static final int numTopicBits = 10;	// max. 1024 topics

		private final Corpus<ProcessedMessage> slice;
		private final DiscreteDistributionSampler sampler;

		public SharedSamplingState(SamplingState other, Corpus<ProcessedMessage> slice,
				DiscreteDistributionSampler sampler, int expectedWordTopicShadowSize, int expectedSortShadowSize) {
			super(other.isQuery, other.type);
			this.slice = slice;
			this.sampler = sampler;

			// share arrays if elements are updated disjointly by the worker threads
			xz = other.xz;
			perAuthorRecipientData = other.perAuthorRecipientData;

			// deep clone all other arrays
			if (!other.isQuery) {
				numWordTopic = new SparseShadowInt2DArray(other.numWordTopic, expectedWordTopicShadowSize);
				sumNumTopic = other.sumNumTopic.clone();

				if (other.type != SamplingType.REGULAR) {
					sortedWordTopic = new SparseShadowInt2DArray(other.sortedWordTopic, expectedSortShadowSize);
					invSortedWordTopic = new SparseShadowInt2DArray(other.invSortedWordTopic, expectedSortShadowSize);
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
			ParallelART.super.sample(sampler, slice, this);
		}
	}

	/** maximum number of threads to use */
	protected int maxThreads;

	private boolean useExternalMemoryCorpus;
	private boolean keepOriginalCorpus;

	/**
	 * Increments the number of threads based an the free memory in the start up phase. This reduces the memory peak at
	 * the beginning.
	 */
	private boolean useMemorySavingStartup;

	/** number of threads in the last iteration */
	private transient int prevThreads;

	/** memory consumption of a localState */
	private transient long mbPerWorker;

	/** shadow size of the WordTopicShadowArray hash map to reduce the need of rehashes */
	private transient int prevWordTopicShadowSize;

	/** shadow size of the SortShadowArray hash map to reduce the need of rehashes */
	private transient int prevSortShadowSize;

	private transient ExecutorService pool;
	private transient Future<?>[] futures;
	private transient DiscreteDistributionSampler[] samplers;
	private transient List<SharedSamplingState> localStates;
	private transient IntSet changesWordTopic;

	public ParallelART(int numWords, int numPersons, int numTopics, double alpha, double beta) {
		super(numWords, numPersons, numTopics, alpha, beta);
		init();
	}

	public ParallelART(int numWords, int numPersons, int numTopics, double alpha[], double beta[]) {
		super(numWords, numPersons, numTopics, alpha, beta);
		init();
	}

	public ParallelART(int numWords, int numPersons, int numTopics, double alpha[], Double2DArray beta) {
		super(numWords, numPersons, numTopics, alpha, beta);
		init();
	}

	protected ParallelART(ParallelART other) {
		super(other);
		useMemorySavingStartup = other.useMemorySavingStartup;
		useExternalMemoryCorpus = other.useExternalMemoryCorpus;
		init();
	}

	@Override
	public ParallelART createSimilar() {
		return new ParallelART(this);
	}

	private void readObject(ObjectInputStream stream) throws IOException, ClassNotFoundException {
		stream.defaultReadObject();
		init();
	}

	private void init() {
		assert (numTopics < (1 << SharedSamplingState.numTopicBits));
		assert (numWords < (1 << (32 - SharedSamplingState.numTopicBits)));

		changesWordTopic = new IntScatterSet();

		int maxThreads = cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_THREADS,
				Runtime.getRuntime().availableProcessors());
		logger.info("using max. " + maxThreads + " threads");
		setMaxThreads(maxThreads);
		useMemorySavingStartup = cfg.getLocalBooleanProperty("memorySavingStartup", useMemorySavingStartup);
	}

	public void setMaxThreads(int maxThreads) {
		this.maxThreads = maxThreads;
		futures = new Future[maxThreads];
		localStates = new ArrayList<SharedSamplingState>(maxThreads);

		// DiscreteDistributionSampler is not thread safe due to keeping its own PRNG state, so we need one instance per
		// worker.
		samplers = new DiscreteDistributionSampler[maxThreads];
		for (int i = 0; i < maxThreads; i++)
			samplers[i] = new DiscreteDistributionSampler(sampler);
	}

	public void setUseExternalMemoryCorpus(boolean useExternalMemoryCorpus) {
		this.useExternalMemoryCorpus = useExternalMemoryCorpus;
	}

	public void setKeepOriginalCorpus(boolean keepOriginalCorpus) {
		this.keepOriginalCorpus = keepOriginalCorpus;
	}

	public void setMemorySavingStartup(boolean useMemorySavingStartup) {
		this.useMemorySavingStartup = useMemorySavingStartup;
	}

	@Override
	protected void burnIn(DiscreteDistributionSampler sampler, Corpus<ProcessedMessage> corpus, SamplingState state,
			int numIter) {
		prevThreads = 0;
		prevWordTopicShadowSize = 0;
		prevSortShadowSize = 0;

		pool = Executors.newFixedThreadPool(maxThreads);
		ProcessedDocumentReader<ProcessedMessage> reader = new ProcessedMessage.Reader();
		SlicedCorpus<ProcessedMessage> slicedCorpus;
		if (useMemorySavingStartup) {
			slicedCorpus = new SlicedCorpus<ProcessedMessage>(corpus, reader, 1, useExternalMemoryCorpus, true);
		} else {
			slicedCorpus = new SlicedCorpus<ProcessedMessage>(corpus, reader, maxThreads, useExternalMemoryCorpus,
					keepOriginalCorpus);
		}
		super.burnIn(sampler, slicedCorpus, state, numIter);
		pool.shutdown();
	}

	@Override
	protected void sample(DiscreteDistributionSampler sampler, Corpus<ProcessedMessage> corpus, SamplingState state) {
		if (pool.isShutdown()) {	// calling sample after burn-in
			super.sample(sampler, corpus, state);
			return;
		}

		long usedMBs = 0;
		if (useMemorySavingStartup && (prevThreads < maxThreads))
			usedMBs = getUsedMemory();

		int numThreads;
		if ((prevThreads == maxThreads) || !useMemorySavingStartup) {
			numThreads = maxThreads;
		} else if (prevThreads == 0) {
			numThreads = 1;
		} else {
			// No need for calling System.gc(), has already been done be getUsedMemory.
			long freeMB = Runtime.getRuntime().freeMemory() / (1024L * 1024L);
			numThreads = Math.max((int) Math.min(Math.min(maxThreads, (freeMB / mbPerWorker) - 1), prevThreads * 2), 1);
			if (numThreads > prevThreads) {
				logger.info("freeMB: " + freeMB + " mbPerWorker: " + mbPerWorker);
				logger.info("Updating threads from " + prevThreads + " to " + numThreads + " (max " + maxThreads + ")");
			} else
				numThreads = prevThreads;
		}

		// start the work
		List<Corpus<ProcessedMessage>> slices = corpus.slice(numThreads, keepOriginalCorpus ||
				(useMemorySavingStartup && (numThreads < maxThreads)));
		for (int i = 0; i < numThreads; i++) {
			SharedSamplingState localState = new SharedSamplingState(state, slices.get(i), samplers[i],
					prevWordTopicShadowSize, prevSortShadowSize);
			localStates.add(localState);
			futures[i] = pool.submit(localState);
		}

		// wait for completion of all work items
		try {
			for (int i = 0; i < numThreads; i++)
				futures[i].get();
		} catch (Exception ex) {
			throw new RuntimeException(ex);
		}

		if (useMemorySavingStartup && (numThreads < maxThreads)) {
			long usedMbAfterStates = getUsedMemory();	// memory calculation including local states
			mbPerWorker = Math.max(1, (usedMbAfterStates - usedMBs) / numThreads);	// prevent division by zero
		}

		// update numTopicWord, sumNumTopic, sort order and bNormSquared
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

				if (state.type != SamplingType.REGULAR) {
					if (sumDiff > 0)
						state.updateSortInc(topic, wordIndex);
					else
						state.updateSortDec(topic, wordIndex);

					double old = prevNum + beta.get(topic, wordIndex);
					state.bNormSquared[wordIndex] += ((old + sumDiff) * (old + sumDiff)) - (old * old);
				}
			}

			changesWordTopic.clear();

			if (state.type != SamplingType.REGULAR) {
				prevWordTopicShadowSize = 0;
				prevSortShadowSize = 0;
				for (SharedSamplingState localState : localStates) {
					prevWordTopicShadowSize = Math.max(prevWordTopicShadowSize, localState.numWordTopic.internalSize());
					prevSortShadowSize = Math.max(prevSortShadowSize, localState.sortedWordTopic.internalSize());
				}
			}
		}

		prevThreads = numThreads;

		localStates.clear();
	}

	@Override
	public void configureForQuerying() {
		setMemorySavingStartup(false);
		setUseExternalMemoryCorpus(false);
		super.configureForQuerying();
	}

	@Override
	public double queryLogLikelihoodLR(Corpus<ProcessedMessage> corpus, final ResampleMode mode) {
		double logLikelihood = 0.0;
		List<Corpus<ProcessedMessage>> slices = corpus.slice(maxThreads, keepOriginalCorpus);
		List<Future<Double>> futures = new ArrayList<Future<Double>>(maxThreads);
		ExecutorService pool = Executors.newFixedThreadPool(maxThreads);
		try {
			for (final Corpus<ProcessedMessage> slice : slices) {
				futures.add(pool.submit(new Callable<Double>() {
					@Override
					public Double call() throws Exception {
						return ParallelART.super.queryLogLikelihoodLR(slice, mode);
					}
				}));
			}
			for (Future<Double> f : futures) {
				try {
					logLikelihood += f.get();
				} catch (Exception ex) {
					throw new RuntimeException(ex);
				}
			}
		} finally {
			pool.shutdown();
		}
		return logLikelihood;
	}

	/**
	 * Calls System.gc() and returns amount of used memory in MB.
	 */
	private static long getUsedMemory() {
	    Runtime runtime = Runtime.getRuntime();
		runtime.gc();
		return (runtime.maxMemory() - runtime.freeMemory()) / (1024L * 1024L);
	}

}
