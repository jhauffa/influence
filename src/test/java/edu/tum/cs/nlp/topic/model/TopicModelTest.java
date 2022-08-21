package edu.tum.cs.nlp.topic.model;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.distribution.ZipfDistribution;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.simple.RandomSource;
import org.junit.Test;

import edu.tum.cs.math.dist.DirichletDistribution;
import edu.tum.cs.math.dist.DiscreteDistribution;
import edu.tum.cs.math.dist.DiscreteDistributionSampler;
import edu.tum.cs.nlp.corpus.Document;
import edu.tum.cs.nlp.corpus.DocumentCorpus;
import edu.tum.cs.nlp.corpus.Index;
import edu.tum.cs.nlp.corpus.Message;
import edu.tum.cs.nlp.corpus.MessageCorpus;
import edu.tum.cs.nlp.corpus.MessageSequenceCorpus;
import edu.tum.cs.nlp.corpus.ProcessedDocument;
import edu.tum.cs.nlp.corpus.TokenizedMessage;
import edu.tum.cs.nlp.topic.model.LDABasedGibbsSamplingModel.SamplingType;
import edu.tum.cs.util.arrays.Object2DArray;
import edu.tum.cs.util.arrays.Sparse2DIterator;
import edu.tum.cs.util.arrays.SparseObject2DArray;
import edu.tum.cs.util.io.Serializer;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class TopicModelTest {

	private static final Logger logger = Logger.getLogger(TopicModelTest.class.getName());

	public static class RandomDocumentBuilder {
		private final long author;
		private final List<Long> recipients;
		private final Date date;
		private final int docLength;
		private final List<String> words;

		private RandomDocumentBuilder(UniformRandomProvider rand, int meanLength, int author, Integer[] recipients,
				long timeStamp) {
			this.author = author;
			this.recipients = new ArrayList<Long>(recipients.length);
			for (Integer recipient : recipients)
				this.recipients.add((long) recipient);
			this.date = new Date(timeStamp);

			// sample document length from a Poisson distribution
			PoissonDistribution p = new PoissonDistribution(meanLength);
			p.reseedRandomGenerator(rand.nextLong());
			docLength = p.sample();
			words = new ArrayList<String>(docLength);
		}

		/** generate a document suitable for LDA */
		public RandomDocumentBuilder(UniformRandomProvider rand, int meanLength, double[] thetaCumul,
				double[][] phiCumul) {
			this(rand, meanLength, 0, new Integer[] { 0 }, 0);

			// sample topic t from discrete distribution theta, then sample word index from discrete distribution phi_t
			for (int i = 0; i < docLength; i++) {
				int topic = DiscreteDistributionSampler.sample1D(rand.nextDouble(), thetaCumul);
				words.add(Integer.toString(DiscreteDistributionSampler.sample1D(rand.nextDouble(), phiCumul[topic])));
			}
		}

		/** generate a message suitable for ART */
		public RandomDocumentBuilder(UniformRandomProvider rand, int meanLength, int author, Integer[] recipients,
				long timeStamp, Object2DArray<double[]> thetaCumul, double[][] phiCumul) {
			this(rand, meanLength, author, recipients, timeStamp);

			// sample recipient uniformly, then sample topic t from theta, then sample word index from phi
			for (int i = 0; i < docLength; i++) {
				int recipientIdx = rand.nextInt(recipients.length);
				int topic = DiscreteDistributionSampler.sample1D(rand.nextDouble(),
						thetaCumul.get(author, recipients[recipientIdx]));
				words.add(Integer.toString(DiscreteDistributionSampler.sample1D(rand.nextDouble(), phiCumul[topic])));
			}
		}

		public Message<Long> buildMessage() {
			return new TokenizedMessage(words, author, recipients, date);
		}
	}

	public static class RandomCorpusBuilder {
		private static final int meanDocumentLength = 300;
		private static final int vocabularyFactor = 10;	// usually 10 - 100
		private static final double vocabularyExponent = 0.4;	// usually 0.4 - 0.6
		private static final double recipientsExponent = 1.86;	// Engel, "Clusters, recipients and reciprocity", 2011

		private final UniformRandomProvider rand;
		public final Index<String> bagOfWords;
		public final Index<Long> bagOfPersons;
		private final List<RandomDocumentBuilder> documents = new ArrayList<RandomDocumentBuilder>();
		private final double alpha;
		private final double beta;
		public final double[][] phi;

		public double[][] thetaLDA;
		public Object2DArray<double[]> thetaART;
		public Object2DArray<double[][]> thetaSocial;
		public Object2DArray<double[][]> piSocial;

		public RandomCorpusBuilder(UniformRandomProvider rand, int numDocuments, int numPersons, int numTopics,
				double alpha, double beta, boolean fillBag) {
			this.rand = rand;
			this.alpha = alpha;
			this.beta = beta;
			bagOfWords = new Index<String>();
			bagOfPersons = new Index<Long>();

			// compute a plausible number of unique words (vocabulary size) according to Heaps' law
			int numWords = (int) Math.round(vocabularyFactor *
					Math.pow(numDocuments * meanDocumentLength, vocabularyExponent));
			logger.info("generating random corpus with " + numWords + " unique words");

			if (fillBag) {
				for (int i = 0; i < numWords; i++)
					bagOfWords.add(Integer.toString(i));
				for (int i = 0; i < numPersons; i++)
					bagOfPersons.add((long) i);
			}

			// generate topics by repeated sampling from a Dirichlet distribution
			DirichletDistribution dist = new DirichletDistribution(numWords, beta, rand);
			phi = new double[numTopics][];
			for (int i = 0; i < numTopics; i++)
				phi[i] = dist.sample();
		}

		/** Generate new corpus with same topics and hyper-parameters as an existing corpus. */
		public RandomCorpusBuilder(RandomCorpusBuilder parent) {
			bagOfWords = parent.bagOfWords;
			bagOfPersons = parent.bagOfPersons;
			rand = parent.rand;
			alpha = parent.alpha;
			beta = parent.beta;
			phi = parent.phi;
		}

		public void generateDocuments(int numDocuments) {
			double[][] phiCumul = cumulative2D(phi);

			DirichletDistribution dist = new DirichletDistribution(phiCumul.length, alpha, rand);
			thetaLDA = new double[numDocuments][];
			thetaART = null;
			thetaSocial = null;

			for (int i = 0; i < numDocuments; i++) {
				// sample topic proportions from a Dirichlet distribution
				thetaLDA[i] = dist.sample();
				double[] thetaCumul = DiscreteDistribution.cumulative(thetaLDA[i]);

				documents.add(new RandomDocumentBuilder(rand, meanDocumentLength, thetaCumul, phiCumul));
			}
		}

		public void generateMessages(int numDocuments, int numPersons) {
			double[][] phiCumul = cumulative2D(phi);

			DiscreteDistributionSampler samp = new DiscreteDistributionSampler(rand.nextLong());
			DirichletDistribution dist = new DirichletDistribution(phiCumul.length, alpha, rand);
			thetaLDA = null;
			thetaSocial = null;
			thetaART = new SparseObject2DArray<double[]>(numPersons, numPersons);
			Object2DArray<double[]> thetaCumulART = new SparseObject2DArray<double[]>(numPersons, numPersons);

			for (int i = 0; i < numDocuments; i++) {
				int author = rand.nextInt(numPersons);

				// sample number of recipients from Zipf distribution
				ZipfDistribution z = new ZipfDistribution(numPersons, recipientsExponent);
				z.reseedRandomGenerator(rand.nextLong());
				int numRecipients = z.sample();

				// sample recipients uniformly
				int[] r = samp.sample1DUniformWithoutReplacement(numPersons, numRecipients);
				Integer[] recipients = new Integer[numRecipients];
				for (int j = 0; j < numRecipients; j++) {
					recipients[j] = r[j];

					// sample topic proportions for author-recipient pair if not yet present
					double[] theta = thetaART.get(author, recipients[j]);
					if (theta == null) {
						theta = dist.sample();
						thetaART.set(author, recipients[j], theta);
						thetaCumulART.set(author, recipients[j], DiscreteDistribution.cumulative(theta));
					}
				}

				documents.add(new RandomDocumentBuilder(rand, meanDocumentLength, author, recipients, 0, thetaCumulART,
						phiCumul));
			}
		}

		public void generateMessageSequences(int numDocuments, int numPersons, int numPairs, int numStates,
				double gamma) {
			double[][] phiCumul = cumulative2D(phi);

			DirichletDistribution distTheta = new DirichletDistribution(phiCumul.length, alpha, rand);
			DirichletDistribution distPi = new DirichletDistribution(numStates, gamma, rand);

			thetaLDA = null;
			thetaART = null;
			thetaSocial = new SparseObject2DArray<double[][]>(numPersons, numPersons);
			piSocial = new SparseObject2DArray<double[][]>(numPersons, numPersons);

			for (int i = 0; i < numPairs; i++) {
				int author = rand.nextInt(numPersons);
				int recipient = rand.nextInt(numPersons);

				double[][] theta = new double[numStates][];
				double[][] thetaCumul = new double[numStates][];
				for (int j = 0; j < theta.length; j++) {
					theta[j] = distTheta.sample();
					thetaCumul[j] = DiscreteDistribution.cumulative(theta[j]);
				}
				thetaSocial.set(author, recipient, theta);

				double[][] pi = new double[numStates][];
				double[][] piCumul = new double[numStates][];
				for (int j = 0; j < pi.length; j++) {
					pi[j] = distPi.sample();
					piCumul[j] = DiscreteDistribution.cumulative(pi[j]);
				}
				piSocial.set(author, recipient, pi);

				int state = 0;
				long date = 0;
				for (int j = 0; j < (numDocuments / numPairs); j++) {
					state = DiscreteDistributionSampler.sample1D(rand.nextDouble(), piCumul[state]);

					Object2DArray<double[]> thetaCumulLocal = new SparseObject2DArray<double[]>(numPersons, numPersons);
					thetaCumulLocal.set(author, recipient, thetaCumul[state]);
					documents.add(new RandomDocumentBuilder(rand, meanDocumentLength, author,
							new Integer[] { recipient }, date, thetaCumulLocal, phiCumul));
					date += 1000;
				}
			}
		}

		public int getNumWords() {
			return phi[0].length;
		}

		public MessageSequenceCorpus<Long, Message<Long>> buildMessageSequenceCorpus() {
			List<Message<Long>> rawMessages = new ArrayList<Message<Long>>(documents.size());
			for (RandomDocumentBuilder doc : documents)
				rawMessages.add(doc.buildMessage());
			return new MessageSequenceCorpus<Long, Message<Long>>(bagOfWords, bagOfPersons, rawMessages);
		}

		public MessageCorpus<Long, Message<Long>> buildMessageCorpus() {
			List<Message<Long>> rawMessages = new ArrayList<Message<Long>>(documents.size());
			for (RandomDocumentBuilder doc : documents)
				rawMessages.add(doc.buildMessage());
			return new MessageCorpus<Long, Message<Long>>(bagOfWords, bagOfPersons, rawMessages);
		}

		public DocumentCorpus buildDocumentCorpus() {
			List<Document> rawDocuments = new ArrayList<Document>(documents.size());
			for (RandomDocumentBuilder doc : documents)
				rawDocuments.add(doc.buildMessage());
			return new DocumentCorpus(bagOfWords, rawDocuments);
		}
	}

	private static final long seed = 4889128L;

	private static void assertBelowThreshold(double value, double threshold) {
		assertTrue(value + " exceeds threshold of " + threshold, value < threshold);
	}

	private static final double eps = 1E-12;

	public static void checkProbabilityDistr(double[] p) {
		double sum = 0.0;
		for (double v : p) {
			assertTrue(Double.toString(v), (v >= 0.0) && (v <= 1.0));
			sum += v;
		}
		assertEquals(1.0, sum, eps);
	}

	private static final int numDocuments = 400;
	private static final int numPersons = 20;
	private static final int numTopics = 10;
	private static final int numStates = 2;
	private static final double alpha = 50.0 / numTopics;
	private static final double beta = 0.01;
	private static final double gamma = 1.0;
	private static final int numIter = 200;

	private static final double maxAvgPhiJSD = 0.5;
	private static final double maxThetaJSD = 0.3;
	private static final double maxAvgThetaJSD = 0.15;
	private static final double maxAvgPiJSD = 0.2;
	private static final double maxPerplexityDiff = 6.0;
	private static final double minPerplexityRatio = 300.0;

	private static double[][] cumulative2D(double[][] p) {
		double[][] cumul = new double[p.length][];
		for (int i = 0; i < p.length; i++)
			cumul[i] = DiscreteDistribution.cumulative(p[i]);
		return cumul;
	}

	private static int[] matchGreedy(double[][] o1, double[][] o2) {
		int numItems = o1.length;
		int[] invMap = new int[numItems];
		Arrays.fill(invMap, -1);

		for (int i = 0; i < numItems; i++) {
			int bestTopicIdx = -1;
			double minDist = Double.MAX_VALUE;
			for (int j = 0; j < numItems; j++) {
				if (invMap[j] != -1)
					continue;

				double dist = DiscreteDistribution.distJS2(o1[i], o2[j]);
				if (dist < minDist) {
					minDist = dist;
					bestTopicIdx = j;
				}
			}

			invMap[bestTopicIdx] = i;
		}

		int[] map = new int[numItems];
		for (int i = 0; i < numItems; i++)
			map[invMap[i]] = i;

		// test for consistency: expect bijective mapping
		boolean[] isMapped = new boolean[numItems];
		for (int i = 0; i < numItems; i++) {
			int m = map[i];
			assertTrue((m >= 0) && (m < numItems));
			assertFalse(isMapped[m]);
			isMapped[m] = true;
		}

		return map;
	}

	private static double[] mapVector(double[] v, int[] map) {
		double[] mapped = new double[v.length];
		for (int i = 0; i < v.length; i++)
			mapped[i] = v[map[i]];
		return mapped;
	}

	private static void testLDAInstance(LDA lda, RandomCorpusBuilder refCorpus, UniformRandomProvider rand)
			throws IOException {
		DocumentCorpus corpus = refCorpus.buildDocumentCorpus();
		int numDocuments = corpus.size();

		// test training
		lda.setSampler(new DiscreteDistributionSampler(rand.nextLong()));
		lda.setBurnIn(numIter);
		lda.setBurnInReportRate(numIter / 10);
		lda.setBurnInSaveRate(-1);
		LDA.SamplingState state = lda.train(corpus);

		// save and reload; resuming training is not implemented
		File tempFile = File.createTempFile("lda-test-", ".ser");
		tempFile.deleteOnExit();
		Serializer.saveObjectToFile(lda, tempFile);
		lda = Serializer.loadObjectFromFile(tempFile);

		double[][] estPhi = lda.getTopicWordDistr();
		for (double[] topicDistr : estPhi)
			checkProbabilityDistr(topicDistr);
		int[] topicMap = matchGreedy(refCorpus.phi, estPhi);
		RandomCorpusBuilder corpusRand = new RandomCorpusBuilder(rand, numDocuments, 0, numTopics, alpha, beta, true);
		corpusRand.generateDocuments(numDocuments);
		double sumPhiUnmappedJSD = 0.0;
		double sumPhiMappedJSD = 0.0;
		double sumPhiRandJSD = 0.0;
		for (int i = 0; i < numTopics; i++) {
			// not testing individual phi: unstable due to low number of iterations
			sumPhiUnmappedJSD += DiscreteDistribution.distJS2(refCorpus.phi[i], estPhi[i]);
			sumPhiMappedJSD += DiscreteDistribution.distJS2(refCorpus.phi[i], estPhi[topicMap[i]]);
			sumPhiRandJSD += DiscreteDistribution.distJS2(refCorpus.phi[i], corpusRand.phi[i]);
		}
		assertTrue(sumPhiMappedJSD < sumPhiUnmappedJSD);	// greedy matching should improve avg. JSD
		assertTrue(sumPhiMappedJSD < sumPhiRandJSD);
		assertBelowThreshold(sumPhiMappedJSD / numTopics, maxAvgPhiJSD);

		double[][] estTheta = state.estimateTheta();
		for (double[] docTheta : estTheta)
			checkProbabilityDistr(docTheta);
		double sumThetaJSD = 0.0;
		double sumThetaRandJSD = 0.0;
		for (int i = 0; i < numDocuments; i++) {
			double d = DiscreteDistribution.distJS2(refCorpus.thetaLDA[i], mapVector(estTheta[i], topicMap));
			assertBelowThreshold(d, maxThetaJSD);
			sumThetaJSD += d;
			sumThetaRandJSD += DiscreteDistribution.distJS2(refCorpus.thetaLDA[i], corpusRand.thetaLDA[i]);
		}
		assertTrue(sumThetaJSD < sumThetaRandJSD);
		assertBelowThreshold(sumThetaJSD / numDocuments, maxAvgThetaJSD);

		int numTokens = corpus.countTokens();
		double trueLL = LDA.estimateLogLikelihood(corpus, refCorpus.phi, refCorpus.thetaLDA);
		double estLL = LDA.estimateLogLikelihood(corpus, estPhi, estTheta);
		assertTrue((trueLL < 0.0) && (estLL < 0.0));
		double truePerplexity = DiscreteDistribution.perplexity(trueLL, numTokens);
		double estPerplexity = DiscreteDistribution.perplexity(estLL, numTokens);
		assertBelowThreshold(Math.abs(truePerplexity - estPerplexity), maxPerplexityDiff);

		// test querying
		RandomCorpusBuilder refCorpusQuery = new RandomCorpusBuilder(refCorpus);
		refCorpusQuery.generateDocuments(numDocuments / 10);
		DocumentCorpus corpusQuery = refCorpusQuery.buildDocumentCorpus();
		numDocuments = corpusQuery.size();
		numTokens = corpusQuery.countTokens();

		lda.setBurnIn(numIter / 2);
		lda.setQuerySamples(10);
		LDA.QueryResult result = lda.query(corpusQuery);
		estTheta = result.estimateTheta();
		for (double[] docTheta : estTheta)
			checkProbabilityDistr(docTheta);
		sumThetaJSD = 0.0;
		for (int i = 0; i < numDocuments; i++) {
			double d = DiscreteDistribution.distJS2(refCorpusQuery.thetaLDA[i], mapVector(estTheta[i], topicMap));
			assertBelowThreshold(d, maxThetaJSD);
			sumThetaJSD += d;
		}
		assertBelowThreshold(sumThetaJSD / numDocuments, maxAvgThetaJSD);

		RandomCorpusBuilder refCorpusQueryRand = new RandomCorpusBuilder(rand, numDocuments, 0, numTopics, alpha, beta,
				true);
		refCorpusQueryRand.generateDocuments(numDocuments / 10);
		DocumentCorpus corpusQueryRand = refCorpusQueryRand.buildDocumentCorpus();
		lda.setLikelihoodSamples(10);
		estLL = lda.queryLogLikelihood(corpusQuery);
		double estLLRandPhi = lda.queryLogLikelihood(corpusQueryRand);
		assertTrue((estLL < 0.0) && (estLLRandPhi < 0.0));
		estPerplexity = DiscreteDistribution.perplexity(estLL, corpusQuery.countTokens());
		double estPerplexityRandPhi = DiscreteDistribution.perplexity(estLLRandPhi, corpusQueryRand.countTokens());
		double r = estPerplexityRandPhi / estPerplexity;
		assertTrue(Double.toString(r), r > minPerplexityRatio);
	}

	@Test
	public void testLDA() throws IOException {
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
		RandomCorpusBuilder corpus = new RandomCorpusBuilder(rand, numDocuments, 0, numTopics, alpha, beta, true);
		corpus.generateDocuments(numDocuments);

		LDA lda = new LDA(corpus.getNumWords(), numTopics, alpha, beta);
		logger.info("testing LDA with regular sampling");
		lda.setSamplingType(SamplingType.REGULAR);
		testLDAInstance(lda, corpus, rand);
		logger.info("testing LDA with FastLDA sampling");
		lda.setSamplingType(SamplingType.FASTLDA);
		testLDAInstance(lda, corpus, rand);

		ParallelLDA plda = new ParallelLDA(corpus.getNumWords(), numTopics, alpha, beta);
		plda.setKeepOriginalCorpus(true);
		logger.info("testing parallel LDA");
		plda.setNumThreads(2);	// reduce overhead and variability caused by having many small slices
		testLDAInstance(plda, corpus, rand);
	}

	private static void testARTInstance(ART art, RandomCorpusBuilder refCorpus, int numIter, UniformRandomProvider rand,
			boolean dataDrivenTheta) throws IOException {
		MessageCorpus<Long, Message<Long>> corpus = refCorpus.buildMessageCorpus();

		// test training
		art.setSampler(new DiscreteDistributionSampler(rand.nextLong()));
		art.setBurnIn(numIter);
		art.setBurnInReportRate(numIter / 10);
		art.setBurnInSaveRate(-1);
		ART.SamplingState state = art.train(corpus);

		// save, reload, and resume training for a few iterations; should have no adverse effect
		File tempFile = File.createTempFile("art-test-", ".ser");
		tempFile.deleteOnExit();
		Serializer.saveObjectToFile(art, tempFile);
		art = Serializer.loadObjectFromFile(tempFile);
		art.resumeTraining(corpus, 10);

		double[][] estPhi = art.getTopicWordDistr();
		for (double[] topicDistr : estPhi)
			checkProbabilityDistr(topicDistr);
		int[] topicMap = matchGreedy(refCorpus.phi, estPhi);
		RandomCorpusBuilder corpusRand = new RandomCorpusBuilder(rand, numDocuments, numPersons, numTopics, alpha, beta,
				true);
		corpusRand.generateMessages(numDocuments, numPersons);
		double sumPhiJSD = 0.0;
		double sumPhiRandJSD = 0.0;
		for (int i = 0; i < numTopics; i++) {
			// not testing individual phi: unstable due to low number of iterations
			sumPhiJSD += DiscreteDistribution.distJS2(refCorpus.phi[i], estPhi[topicMap[i]]);
			sumPhiRandJSD += DiscreteDistribution.distJS2(refCorpus.phi[i], corpusRand.phi[i]);
		}
		assertTrue(sumPhiJSD < sumPhiRandJSD);
		assertBelowThreshold(sumPhiJSD / numTopics, maxAvgPhiJSD);

		Object2DArray<double[]> estTheta = state.estimateTheta();
		Sparse2DIterator<double[]> it = estTheta.sparseIterator();
		while (it.hasNext()) {
			it.advance();
			checkProbabilityDistr(it.getValue());
		}
		double sumThetaJSD = 0.0;
		Sparse2DIterator<double[]> itCorpus = refCorpus.thetaART.sparseIterator();
		while (itCorpus.hasNext()) {
			itCorpus.advance();
			double[] curEstTheta = estTheta.get(itCorpus.getRow(), itCorpus.getColumn());
			double d = DiscreteDistribution.distJS2(itCorpus.getValue(), mapVector(curEstTheta, topicMap));
			assertBelowThreshold(d, maxThetaJSD);
			sumThetaJSD += d;
		}
		// not comparing to random as random corpus has different author-recipient pairs
		assertBelowThreshold(sumThetaJSD / estTheta.size(), maxAvgThetaJSD);

		// test the aggregate topic distributions
		double[][] aggTheta = state.estimateAggregatedTheta(true, null, true, null);	// global
		checkProbabilityDistr(aggTheta[0]);
		aggTheta = state.estimateAggregatedTheta(false, null, true, null);	// per author
		for (double[] theta : aggTheta)
			checkProbabilityDistr(theta);
		double[][] aggThetaSingle = state.estimateAggregatedTheta(true, new HashSet<Integer>(Arrays.asList(0)),
				true, null);	// global, selecting first author only
		checkProbabilityDistr(aggThetaSingle[0]);
		for (int i = 0; i < aggThetaSingle[0].length; i++)
			assertEquals(aggTheta[0][i], aggThetaSingle[0][i], eps);
		aggTheta = state.estimateAggregatedTheta(true, null, false, null);	// per recipient
		for (double[] theta : aggTheta)
			checkProbabilityDistr(theta);

		double trueLL = ART.estimateLogLikelihood(corpus, refCorpus.phi, refCorpus.thetaART);
		double estLL = ART.estimateLogLikelihood(corpus, estPhi, estTheta);
		assertTrue((trueLL < 0.0) && (estLL < 0.0));
		int numTokens = corpus.countTokens();
		double truePerplexity = DiscreteDistribution.perplexity(trueLL, numTokens);
		double estPerplexity = DiscreteDistribution.perplexity(estLL, numTokens);
		assertBelowThreshold(Math.abs(truePerplexity - estPerplexity), maxPerplexityDiff);

		// test querying
		RandomCorpusBuilder refCorpusQuery = new RandomCorpusBuilder(refCorpus);
		refCorpusQuery.generateMessages(numDocuments / 10, numTopics);
		MessageCorpus<Long, Message<Long>> corpusQuery = refCorpusQuery.buildMessageCorpus();
		art.setBurnIn(numIter / 2);
		art.setQuerySamples(10);
		ART.QueryResult result = art.query(corpusQuery);
		estTheta = result.estimateTheta();
		it = estTheta.sparseIterator();
		while (it.hasNext()) {
			it.advance();
			checkProbabilityDistr(it.getValue());
		}
		sumThetaJSD = 0.0;
		itCorpus = refCorpusQuery.thetaART.sparseIterator();
		while (itCorpus.hasNext()) {
			itCorpus.advance();
			double[] curEstTheta = estTheta.get(itCorpus.getRow(), itCorpus.getColumn());
			double d = DiscreteDistribution.distJS2(itCorpus.getValue(), mapVector(curEstTheta, topicMap));
			assertBelowThreshold(d, maxThetaJSD);
			sumThetaJSD += d;
		}
		assertBelowThreshold(sumThetaJSD / estTheta.size(), maxAvgThetaJSD);

		RandomCorpusBuilder refCorpusQueryRand = new RandomCorpusBuilder(rand, numDocuments, numPersons, numTopics,
				alpha, beta, true);
		refCorpusQueryRand.generateMessages(numDocuments / 10, numTopics);
		MessageCorpus<Long, Message<Long>> corpusQueryRand = refCorpusQueryRand.buildMessageCorpus();
		art.setLikelihoodSamples(10);
		estLL = art.queryLogLikelihood(corpusQuery);
		double estLLRandPhi = art.queryLogLikelihood(corpusQueryRand);
		assertTrue((estLL < 0.0) && (estLLRandPhi < 0.0));
		estPerplexity = DiscreteDistribution.perplexity(estLL, corpusQuery.countTokens());
		double estPerplexityRandPhi = DiscreteDistribution.perplexity(estLLRandPhi, corpusQueryRand.countTokens());
		double r = estPerplexityRandPhi / estPerplexity;
		assertTrue(Double.toString(r), r > minPerplexityRatio);
	}

	@Test
	public void testART() throws IOException {
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
		RandomCorpusBuilder corpus = new RandomCorpusBuilder(rand, numDocuments, numPersons, numTopics, alpha, beta,
				true);
		corpus.generateMessages(numDocuments, numPersons);

		ART art = new ART(corpus.getNumWords(), numPersons, numTopics, alpha, beta);
		logger.info("testing ART with regular sampling");
		art.setSamplingType(SamplingType.REGULAR);
		testARTInstance(art, corpus, numIter, rand, false);
		logger.info("testing ART with FastLDA sampling");
		art.setSamplingType(SamplingType.FASTLDA);
		testARTInstance(art, corpus, numIter, rand, false);
		logger.info("testing ART with FastLDA sampling and data driven uniform alpha");
		art.setSamplingType(SamplingType.FASTLDA);
		art.setAlphaSamplingType(ART.AlphaType.DATA_DRIVEN_UNIFORM);
		testARTInstance(art, corpus, numIter, rand, true);
		logger.info("testing ART with FastLDA sampling and data driven non-uniform alpha");
		art.setSamplingType(SamplingType.FASTLDA);
		art.setAlphaSamplingType(ART.AlphaType.DATA_DRIVEN_NON_UNIFORM);
		testARTInstance(art, corpus, numIter, rand, true);

		ParallelART part = new ParallelART(corpus.getNumWords(), numPersons, numTopics, alpha, beta);
		part.setKeepOriginalCorpus(true);
		logger.info("testing parallel ART");
		part.setMaxThreads(2);
		testARTInstance(part, corpus, numIter, rand, false);
		logger.info("testing parallel ART with memory saving startup");
		part.setMaxThreads(10);
		part.setMemorySavingStartup(true);
		testARTInstance(part, corpus, numIter, rand, false);
	}

	private static void testSocialInstance(SocialTopicModel soc, RandomCorpusBuilder refCorpus, int numIter,
			UniformRandomProvider rand) throws IOException {
		MessageSequenceCorpus<Long, Message<Long>> corpus = refCorpus.buildMessageSequenceCorpus();

		// test training
		soc.setSampler(new DiscreteDistributionSampler(rand.nextLong()));
		soc.setBurnIn(numIter);
		soc.setBurnInReportRate(numIter / 10);
		soc.setBurnInSaveRate(-1);
		SocialTopicModel.SamplingState state = soc.train(corpus);

		// save, reload, and resume training for a few iterations; should have no adverse effect
		File tempFile = File.createTempFile("soc-test-", ".ser");
		tempFile.deleteOnExit();
		Serializer.saveObjectToFile(soc, tempFile);
		soc = Serializer.loadObjectFromFile(tempFile);
		soc.resumeTraining(corpus, 10);

		double[][] estPhi = soc.getTopicWordDistr();
		for (double[] topicDistr : estPhi)
			checkProbabilityDistr(topicDistr);
		int[] topicMap = matchGreedy(refCorpus.phi, estPhi);
		RandomCorpusBuilder corpusRand = new RandomCorpusBuilder(rand, numDocuments, numPersons, numTopics, alpha, beta,
				true);
		corpusRand.generateMessageSequences(numDocuments, numPersons, numPersons / 2, numStates, gamma);
		double sumPhiJSD = 0.0;
		double sumPhiRandJSD = 0.0;
		for (int i = 0; i < numTopics; i++) {
			// not testing individual phi: unstable due to low number of iterations
			sumPhiJSD += DiscreteDistribution.distJS2(refCorpus.phi[i], estPhi[topicMap[i]]);
			sumPhiRandJSD += DiscreteDistribution.distJS2(refCorpus.phi[i], corpusRand.phi[i]);
		}
		assertTrue(sumPhiJSD < sumPhiRandJSD);
		assertBelowThreshold(sumPhiJSD / numTopics, maxAvgPhiJSD);

		Object2DArray<double[][]> estTheta = state.estimateTheta();
		Sparse2DIterator<double[][]> it = estTheta.sparseIterator();
		while (it.hasNext()) {
			it.advance();
			for (double[] p : it.getValue())
				checkProbabilityDistr(p);
		}
		Object2DArray<int[]> stateMap = new SparseObject2DArray<int[]>(numPersons, numPersons);
		double sumThetaJSD = 0.0;
		Sparse2DIterator<double[][]> itCorpus = refCorpus.thetaSocial.sparseIterator();
		while (itCorpus.hasNext()) {
			itCorpus.advance();
			double[][] curEstTheta = estTheta.get(itCorpus.getRow(), itCorpus.getColumn());
			int[] curStateMap = matchGreedy(curEstTheta, itCorpus.getValue());
			for (int i = 0; i < numStates; i++) {
				double d = DiscreteDistribution.distJS2(itCorpus.getValue()[i],
						mapVector(curEstTheta[curStateMap[i]], topicMap));
				assertBelowThreshold(d, maxThetaJSD);
				sumThetaJSD += d;
			}
			stateMap.set(itCorpus.getRow(), itCorpus.getColumn(), curStateMap);
		}
		// not comparing to random as random corpus has different author-recipient pairs
		assertBelowThreshold(sumThetaJSD / (estTheta.size() * numStates), maxAvgThetaJSD);

		Object2DArray<double[][]> estPi = state.estimatePi();
		it = estPi.sparseIterator();
		while (it.hasNext()) {
			it.advance();
			for (double[] p : it.getValue())
				checkProbabilityDistr(p);
		}
		double sumPiJSD = 0.0;
		itCorpus = refCorpus.piSocial.sparseIterator();
		while (itCorpus.hasNext()) {
			itCorpus.advance();
			double[][] curEstPi = estPi.get(itCorpus.getRow(), itCorpus.getColumn());
			int[] curStateMap = stateMap.get(itCorpus.getRow(), itCorpus.getColumn());
			for (int i = 0; i < numStates; i++) {
				double d = DiscreteDistribution.distJS2(itCorpus.getValue()[i],
						mapVector(curEstPi[curStateMap[i]], curStateMap));
				sumPiJSD += d;
			}
		}
		assertBelowThreshold(sumPiJSD / (estPi.size() * numStates), maxAvgPiJSD);

		double trueLL = SocialTopicModel.estimateLogLikelihood(corpus, refCorpus.phi, refCorpus.thetaSocial,
				refCorpus.piSocial, numStates);
		double estLL = SocialTopicModel.estimateLogLikelihood(corpus, estPhi, estTheta, estPi, numStates);
		assertTrue((trueLL < 0.0) && (estLL < 0.0));
		int numTokens = corpus.countTokens();
		double truePerplexity = DiscreteDistribution.perplexity(trueLL, numTokens);
		double estPerplexity = DiscreteDistribution.perplexity(estLL, numTokens);
		assertBelowThreshold(Math.abs(truePerplexity - estPerplexity), maxPerplexityDiff);

		// querying not implemented yet; only test estimation of held-out likelihood
		RandomCorpusBuilder refCorpusQuery = new RandomCorpusBuilder(refCorpus);
		refCorpusQuery.generateMessageSequences(numDocuments / 10, numPersons, numPersons / 2, numStates, gamma);
		MessageSequenceCorpus<Long, Message<Long>> corpusQuery = refCorpusQuery.buildMessageSequenceCorpus();
		RandomCorpusBuilder refCorpusQueryRand = new RandomCorpusBuilder(rand, numDocuments, numPersons, numTopics,
				alpha, beta, true);
		refCorpusQueryRand.generateMessageSequences(numDocuments / 10, numPersons, numPersons / 2, numStates, gamma);
		MessageSequenceCorpus<Long, Message<Long>> corpusQueryRand = refCorpusQueryRand.buildMessageSequenceCorpus();

		soc.setLikelihoodSamples(10);
		estLL = soc.queryLogLikelihood(corpusQuery);
		double estLLRandPhi = soc.queryLogLikelihood(corpusQueryRand);
		assertTrue("LL = " + estLL + ", LL (rand.) = " + estLLRandPhi, (estLL < 0.0) && (estLLRandPhi < 0.0));
		estPerplexity = DiscreteDistribution.perplexity(estLL, corpusQuery.countTokens());
		double estPerplexityRandPhi = DiscreteDistribution.perplexity(estLLRandPhi, corpusQueryRand.countTokens());
		double r = estPerplexityRandPhi / estPerplexity;
		assertTrue(Double.toString(r), r > minPerplexityRatio);
	}

	@Test
	public void testSocial() throws IOException {
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);

		RandomCorpusBuilder corpus = new RandomCorpusBuilder(rand, numDocuments, numPersons, numTopics, alpha, beta,
				true);
		corpus.generateMessageSequences(numDocuments, numPersons, numPersons / 2, numStates, gamma);

		SocialTopicModel soc = new SocialTopicModel(corpus.getNumWords(), numPersons, numTopics, numStates,
				alpha, beta, gamma);
		logger.info("testing social topic model with regular sampling");
		soc.setSamplingType(SamplingType.REGULAR);
		testSocialInstance(soc, corpus, numIter, rand);
		logger.info("testing social topic model with FastLDA sampling");
		soc.setSamplingType(SamplingType.FASTLDA);
		testSocialInstance(soc, corpus, numIter, rand);

		ParallelSocialTopicModel psoc = new ParallelSocialTopicModel(corpus.getNumWords(), numPersons, numTopics,
				numStates, alpha, beta, gamma);
		psoc.setKeepOriginalCorpus(true);
		logger.info("testing parallel social topic model");
		psoc.setNumThreads(2);
		testSocialInstance(psoc, corpus, numIter, rand);
	}

	@Test
	public void testMessageSequenceShuffle() {
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);

		RandomCorpusBuilder builder = new RandomCorpusBuilder(rand, numDocuments, numPersons, numTopics, alpha, beta,
				true);
		builder.generateMessageSequences(numDocuments, numPersons, numPersons / 2, numStates, gamma);

		MessageSequenceCorpus<Long, Message<Long>> corpus = builder.buildMessageSequenceCorpus();
		int numRefSeq = 0, numRefMsg = 0;
		Iterator<MessageSequenceCorpus.MessageSequence> seqIt = corpus.sequenceIterator();
		while (seqIt.hasNext()) {
			MessageSequenceCorpus.MessageSequence seq = seqIt.next();
			while (seq.hasNext()) {
				seq.next();
				numRefMsg++;
			}
			numRefSeq++;
		}

		corpus.shuffle();
		int numSeq = 0, numMsg = 0;
		seqIt = corpus.sequenceIterator();
		while (seqIt.hasNext()) {
			MessageSequenceCorpus.MessageSequence seq = seqIt.next();
			while (seq.hasNext()) {
				seq.next();
				numMsg++;
				assertTrue(numMsg <= numRefMsg);
			}
			numSeq++;
			assertTrue(numSeq <= numRefSeq);
		}
		assertEquals(numRefMsg, numMsg);
		assertEquals(numRefSeq, numSeq);
	}

	@Test
	public void testOnlineTopicModel() {
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
		RandomCorpusBuilder refCorpus1 = new RandomCorpusBuilder(rand, numDocuments, numPersons, numTopics, alpha, beta,
				false);
		refCorpus1.generateDocuments(numDocuments);
		DocumentCorpus corpus1 = refCorpus1.buildDocumentCorpus();
		RandomCorpusBuilder refCorpus2 = new RandomCorpusBuilder(refCorpus1);
		refCorpus2.generateDocuments(numDocuments);
		DocumentCorpus corpus2 = refCorpus2.buildDocumentCorpus();

		LDA lda = new LDA(refCorpus1.getNumWords(), numTopics, alpha, beta);
		lda.setBurnIn(numIter * 2);
		lda.setLikelihoodSamples(numIter);
		OnlineTopicModel<ProcessedDocument, LDA> olda = new OnlineTopicModel<ProcessedDocument, LDA>(lda,
				new double[] { 1.0 }, false);
		olda.update(corpus1);
		olda.update(corpus2);
		assertEquals(olda.getCurrentModel(), olda.getModelHistory().get(0));
		assertEquals(2, olda.getModelHistory().size());

		// LL of corpus 1 should be higher for model 1 than for model 2, and vice versa.
		double llModel1 = olda.getModelHistory().get(1).queryLogLikelihoodLR(corpus1, LDA.ResampleMode.LOG);
		double llModel2 = olda.getModelHistory().get(0).queryLogLikelihoodLR(corpus1, LDA.ResampleMode.LOG);
		assertTrue(llModel1 > llModel2);
		llModel1 = olda.getModelHistory().get(1).queryLogLikelihoodLR(corpus2, LDA.ResampleMode.LOG);
		llModel2 = olda.getModelHistory().get(0).queryLogLikelihoodLR(corpus2, LDA.ResampleMode.LOG);
		assertTrue(llModel2 > llModel1);

		// Topics of model 2 should be more similar to those of model 1 than the topics of an individual model fit to
		// corpus 2.
		double[][] phiModel1 = olda.getModelHistory().get(1).getTopicWordDistr();
		double[][] phiModel2 = olda.getModelHistory().get(0).getTopicWordDistr();
		LDA ldaRef = new LDA(refCorpus1.getNumWords(), numTopics, alpha, beta);
		ldaRef.setBurnIn(numIter * 2);
		ldaRef.train(corpus2);
		double[][] phiRef = ldaRef.getTopicWordDistr();

		double avgDist12 = 0.0;
		for (int i = 0; i < numTopics; i++)
			avgDist12 += DiscreteDistribution.distJS2(phiModel1[i], phiModel2[i]);
		avgDist12 /= numTopics;

		int[] topicMap = matchGreedy(phiModel2, phiRef);
		double avgDist2Ref = 0.0;
		for (int i = 0; i < numTopics; i++)
			avgDist2Ref += DiscreteDistribution.distJS2(phiModel2[i], phiRef[topicMap[i]]);
		avgDist2Ref /= numTopics;

		assertTrue(avgDist12 < avgDist2Ref);
	}

	private static final int numTests = 10;

	public static void main(String[] args) throws Exception {
		UniformRandomProvider rand = RandomSource.create(RandomSource.XOR_SHIFT_1024_S, seed);
		RandomCorpusBuilder corpus = new RandomCorpusBuilder(rand, numDocuments, numPersons, numTopics, alpha, beta,
				true);
		corpus.generateMessages(numDocuments, numPersons);

		/* Uncomment the GC calls for monitoring memory usage with VisualVM or similar tools. */
		// System.gc();

		long[] t = new long[numTests];
		for (int i = 0; i < numTests; i++) {
			ParallelART art = new ParallelART(corpus.getNumWords(), numPersons, numTopics, alpha, beta);
			art.setKeepOriginalCorpus(true);
			art.setMaxThreads(4);
			art.setSamplingType(SamplingType.FASTLDA);

			logger.info("test " + (i + 1) + " of " + numTests);
			t[i] = System.currentTimeMillis();
			testARTInstance(art, corpus, numIter, rand, false);
			t[i] = System.currentTimeMillis() - t[i];

			/*
			art = null;
			System.gc();
			*/
		}

		double avg = 0.0;
		for (int i = 0; i < numTests; i++)
			avg += t[i];
		avg /= numTests;
		double stdDev = 0.0;
		for (int i = 0; i < numTests; i++)
			stdDev += Math.pow(t[i] - avg, 2.0);
		stdDev = Math.sqrt((1.0 / (numTests - 1)) * stdDev);
		System.out.println("avg. = " + avg + " ms, std.dev = " + stdDev + " ms");
	}

}
