package edu.tum.cs.nlp.topic.model;

import edu.tum.cs.math.dist.DiscreteDistributionSampler;
import edu.tum.cs.nlp.corpus.Corpus;
import edu.tum.cs.nlp.corpus.ProcessedDocument;
import edu.tum.cs.nlp.corpus.ProcessedDocumentReader;
import edu.tum.cs.nlp.corpus.SlicedCorpus;
import edu.tum.cs.util.ExperimentConfiguration;
import edu.tum.cs.util.arrays.Double2DArray;
import edu.tum.cs.util.arrays.SparseShadowInt2DArray;

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

public class ParallelLDA extends LDA {

	private static final long serialVersionUID = 6060028930747974973L;
	private static final Logger logger = Logger.getLogger(ParallelLDA.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(ParallelLDA.class);

	protected class SharedSamplingState extends SamplingState implements Runnable {
		private static final long serialVersionUID = 1885559150882199348L;
		private static final int numTopicBits = 8;	// max. 255 topics

		private final Corpus<ProcessedDocument> slice;
		private final DiscreteDistributionSampler sampler;

		public SharedSamplingState(SamplingState other, Corpus<ProcessedDocument> slice,
				DiscreteDistributionSampler sampler) {
			super(other.isQuery, other.type);
			this.slice = slice;
			this.sampler = sampler;

			// share arrays if elements are updated disjointly by the worker threads
			z = other.z;
			numDocumentTopic = other.numDocumentTopic;
			aNormSquared = other.aNormSquared;

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
			ParallelLDA.super.sample(sampler, slice, this);
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

	public ParallelLDA(int numWords, int numTopics, double alpha, double beta) {
		super(numWords, numTopics, alpha, beta);
		init();
	}

	public ParallelLDA(int numWords, int numTopics, double alpha[], double beta[]) {
		super(numWords, numTopics, alpha, beta);
		init();
	}

	public ParallelLDA(int numWords, int numTopics, double alpha[], Double2DArray beta) {
		super(numWords, numTopics, alpha, beta);
		init();
	}

	protected ParallelLDA(ParallelLDA other) {
		super(other);
		useExternalMemoryCorpus = other.useExternalMemoryCorpus;
		init();
	}

	@Override
	public ParallelLDA createSimilar() {
		return new ParallelLDA(this);
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
	protected SamplingState burnIn(DiscreteDistributionSampler sampler, Corpus<ProcessedDocument> corpus,
			boolean isQuery) {
		pool = Executors.newFixedThreadPool(numThreads);
		ProcessedDocumentReader<ProcessedDocument> reader = new ProcessedDocument.Reader();
		SlicedCorpus<ProcessedDocument> slices = new SlicedCorpus<ProcessedDocument>(corpus, reader, numThreads,
				useExternalMemoryCorpus, keepOriginalCorpus);
		SamplingState state = super.burnIn(sampler, slices, isQuery);
		pool.shutdown();
		return state;
	}

	@Override
	protected void sample(DiscreteDistributionSampler sampler, Corpus<ProcessedDocument> corpus, SamplingState state) {
		if (pool.isShutdown()) {	// calling sample after burn-in
			super.sample(sampler, corpus, state);
			return;
		}

		List<Corpus<ProcessedDocument>> slices = corpus.slice(numThreads, keepOriginalCorpus);
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

	@Override
	public double queryLogLikelihoodLR(Corpus<ProcessedDocument> corpus, final ResampleMode mode) {
		double logLikelihood = 0.0;
		List<Corpus<ProcessedDocument>> slices = corpus.slice(numThreads, keepOriginalCorpus);
		List<Future<Double>> futures = new ArrayList<Future<Double>>(numThreads);
		ExecutorService pool = Executors.newFixedThreadPool(numThreads);
		try {
			for (final Corpus<ProcessedDocument> slice : slices) {
				futures.add(pool.submit(new Callable<Double>() {
					@Override
					public Double call() throws Exception {
						return ParallelLDA.super.queryLogLikelihoodLR(slice, mode);
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
}
