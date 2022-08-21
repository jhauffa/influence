package edu.tum.cs.nlp.topic;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.SocialMediaDaoFactory;
import edu.tum.cs.db.loader.MessageLoader;
import edu.tum.cs.db.loader.MessageLoaderFactory;
import edu.tum.cs.math.dist.DiscreteDistribution;
import edu.tum.cs.math.dist.Histogram;
import edu.tum.cs.math.dist.HistogramImpl;
import edu.tum.cs.nlp.corpus.Corpus;
import edu.tum.cs.nlp.corpus.Index;
import edu.tum.cs.nlp.corpus.MessageSequenceCorpus;
import edu.tum.cs.nlp.corpus.ProcessedMessage;
import edu.tum.cs.nlp.corpus.TokenizedMessage;
import edu.tum.cs.nlp.topic.model.ART;
import edu.tum.cs.nlp.topic.model.ParallelART;
import edu.tum.cs.nlp.topic.model.ParallelSocialTopicModel;
import edu.tum.cs.nlp.topic.model.SocialTopicModel;
import edu.tum.cs.time.IntervalCalculator;
import edu.tum.cs.util.ExperimentConfiguration;
import edu.tum.cs.util.arrays.Object2DArray;
import edu.tum.cs.util.io.SerializableIntIntScatterMap;
import edu.tum.cs.util.io.Serializer;

public class FitSocialTopicModel {

	public static class SocialTopicModelParameters implements Serializable {
		private static final long serialVersionUID = -1141750214627181099L;

		public final double[][] phi;
		public final Object2DArray<double[][]> pi;
		public final Object2DArray<double[][]> theta;
		public final SerializableIntIntScatterMap seqLengths;

		public SocialTopicModelParameters(double[][] phi, Object2DArray<double[][]> pi,
				Object2DArray<double[][]> theta, SerializableIntIntScatterMap seqLengths) {
			this.phi = phi;
			this.pi = pi;
			this.theta = theta;
			this.seqLengths = seqLengths;
		}
	}

	private static final Logger logger = Logger.getLogger(FitSocialTopicModel.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(FitSocialTopicModel.class);

	private static final int numTopics = cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_TOPICS);
	private static final int numStates = cfg.getLocalIntProperty("numStates");
	private static final String modelPath = cfg.getLocalProperty("modelPath", ".");
	private static final double alpha = 50.0 / numTopics;
	private static final double beta = 0.01;
	private static final double gammaOffDiag = cfg.getLocalDoubleProperty("transitionPrior", 1.0);
	private static final double gammaDiag = cfg.getLocalDoubleProperty("selfTransitionPrior", 1.0);
	private static final int numIter = 2000;
	private static final double testRatio = 0.2;

	private static int numWords, numPersons;

	private static MessageSequenceCorpus<Long, TokenizedMessage> buildCorpus(SocialMediaDao dao, Set<Long> userIds,
			Date startDate, Date endDate, Index<String> bagOfWords, Index<Long> bagOfPersons) throws Exception {
		MessageLoader loader = MessageLoaderFactory.createMessageLoader(dao, userIds);
		List<TokenizedMessage> messages = new ArrayList<TokenizedMessage>();
		int idx = 0;
		for (Long userId : userIds) {
			messages.addAll(loader.loadMessages(userId, startDate, endDate));
			if (++idx % 1000 == 0)
				logger.info("Done with " + idx);
		}
		return new MessageSequenceCorpus<Long, TokenizedMessage>(bagOfWords, bagOfPersons, messages);
	}

	private static List<MessageSequenceCorpus<Long, TokenizedMessage>> buildCorpora(Date startDate, Date endDate)
			throws Exception {
		SocialMediaDao dao = SocialMediaDaoFactory.createDao();
		List<Long> userIds = new ArrayList<Long>(dao.getUserIds(false));
		int numUserIds = userIds.size();
		Collections.shuffle(userIds);
		int splitIdx = (int) (numUserIds * testRatio);
		Set<Long> userIdsTest = new HashSet<Long>(userIds.subList(0, splitIdx));
		Set<Long> userIdsTrain = new HashSet<Long>(userIds.subList(splitIdx, numUserIds));

		logger.info("Loading messages from " + startDate + " to " + endDate);
		Index<String> bagOfWords = new Index<String>();
		Index<Long> bagOfPersons = new Index<Long>();
		List<MessageSequenceCorpus<Long, TokenizedMessage>> corpora =
				new ArrayList<MessageSequenceCorpus<Long, TokenizedMessage>>(2);
		corpora.add(buildCorpus(dao, userIdsTrain, startDate, endDate, bagOfWords, bagOfPersons));
		corpora.add(buildCorpus(dao, userIdsTest, startDate, endDate, bagOfWords, bagOfPersons));
		Serializer.saveObjectToFile(corpora.get(0), new File(modelPath, "corpusTrain.ser.gz"));
		Serializer.saveObjectToFile(corpora.get(1), new File(modelPath, "corpusTest.ser.gz"));

		numWords = bagOfWords.getNumUniqueElements();
		numPersons = bagOfPersons.getNumUniqueElements();
		Serializer.saveObjectToFile(bagOfWords, new File(modelPath, FitART.bowFileName));
		Serializer.saveObjectToFile(bagOfPersons, new File(modelPath, FitART.bopFileName));
		logger.info("Finished loading messages - processed " + numUserIds + " users");
		return corpora;
	}

	private static List<MessageSequenceCorpus<Long, TokenizedMessage>> loadCorpora(String trainFile, String testFile)
			throws Exception {
		Index<String> bagOfWords = Serializer.loadObjectFromFile(new File(modelPath, FitART.bowFileName));
		Index<Long> bagOfPersons = Serializer.loadObjectFromFile(new File(modelPath, FitART.bopFileName));
		numWords = bagOfWords.getNumUniqueElements();
		numPersons = bagOfPersons.getNumUniqueElements();

		logger.info("Loading messages from file");
		List<MessageSequenceCorpus<Long, TokenizedMessage>> corpora =
				new ArrayList<MessageSequenceCorpus<Long, TokenizedMessage>>(2);
		corpora.add(Serializer.<MessageSequenceCorpus<Long, TokenizedMessage>>loadObjectFromFile(new File(trainFile)));
		corpora.add(Serializer.<MessageSequenceCorpus<Long, TokenizedMessage>>loadObjectFromFile(new File(testFile)));
		logger.info("Finished loading messages");
		return corpora;
	}

	private static void evaluateTopicModel(double ll, int numParam, Corpus<ProcessedMessage> corpus, String id) {
		int numObs = corpus.countTokens();
		double bic = (-2.0 * ll) + (numParam * Math.log(numObs));
		double perplexity = DiscreteDistribution.perplexity(ll, numObs);
		System.out.println(id + " corpus: LL = " + ll + ", BIC = " + bic + ", perplexity = " + perplexity);
	}

	private static void fitART(List<MessageSequenceCorpus<Long, TokenizedMessage>> corpora, String id) {
		ParallelART part = new ParallelART(numWords, numPersons, numTopics, alpha, beta);
		part.setBurnIn(numIter);
		part.setBurnInReportRate(50);
		part.setMemorySavingStartup(true);
		part.setUseExternalMemoryCorpus(true);
		part.setKeepOriginalCorpus(true);

		logger.info("Training ART model on " + id + " data with " + numIter + " iterations");
		ART.SamplingState state = part.train(corpora.get(0));
		Serializer.saveObjectToFile(part, new File(modelPath, "art-model-final-" + id + ".ser.gz"));
		part.getIntermediateStateFile().delete();
		logger.info("Finished training ART model on " + id + " data");

		int numUserPairs = state.estimateTheta().size();
		int numParam = (numUserPairs * (numTopics - 1)) + (numTopics * (numWords - 1));
		evaluateTopicModel(part.queryLogLikelihoodLR(corpora.get(1), ART.ResampleMode.LOG), numParam, corpora.get(1),
				"test");
		MessageSequenceCorpus<Long, TokenizedMessage> subCorpus =
				new MessageSequenceCorpus<Long, TokenizedMessage>(corpora.get(1), true);
		evaluateTopicModel(part.queryLogLikelihoodLR(subCorpus, ART.ResampleMode.LOG), numParam, subCorpus,
				"test (length > 1)");
	}

	private static void fitSocialTopicModel(List<MessageSequenceCorpus<Long, TokenizedMessage>> corpora, String id,
			SerializableIntIntScatterMap seqLengths) {
		double[][] gamma = new double[numStates][numStates];
		for (int i = 0; i < numStates; i++) {
			for (int j = 0; j < numStates; j++)
				gamma[i][j] = (i == j) ? gammaDiag : gammaOffDiag;
		}

		ParallelSocialTopicModel psoc = new ParallelSocialTopicModel(numWords, numPersons, numTopics, numStates,
				alpha, beta, gamma);
		psoc.setBurnIn(numIter);
		psoc.setBurnInReportRate(50);
		psoc.setUseExternalMemoryCorpus(true);
		psoc.setKeepOriginalCorpus(true);

		logger.info("Training social topic model on " + id + " data with " + numIter + " iterations");
		logger.info(numStates + " states, gamma (diag.) = " + gammaDiag + ", gamma (off-diag.) = " + gammaOffDiag);
		SocialTopicModel.SamplingState state = psoc.train(corpora.get(0));
		Serializer.saveObjectToFile(psoc, new File(modelPath, "soc-model-final-" + id + "-" + numStates + ".ser.gz"));
		psoc.getIntermediateStateFile().delete();
		logger.info("Finished training social topic model on " + id + " data");

		Serializer.saveObjectToFile(new SocialTopicModelParameters(psoc.getTopicWordDistr(), state.estimatePi(),
				state.estimateTheta(), seqLengths),
				new File(modelPath, "soc-param-" + id + "-" + numStates + ".ser.gz"));

		int numUserPairs = state.estimateTheta().size();
		int numParam = (numStates * numUserPairs * (numTopics + numStates - 2)) + (numTopics * (numWords - 1));
		evaluateTopicModel(psoc.queryLogLikelihoodLR(corpora.get(1), SocialTopicModel.ResampleMode.LOG), numParam,
				corpora.get(1), "test");
		MessageSequenceCorpus<Long, TokenizedMessage> subCorpus =
				new MessageSequenceCorpus<Long, TokenizedMessage>(corpora.get(1), true);
		evaluateTopicModel(psoc.queryLogLikelihoodLR(subCorpus, SocialTopicModel.ResampleMode.LOG), numParam, subCorpus,
				"test (length > 1)");
	}

	private static void printSeqLengthHisto(MessageSequenceCorpus<Long, TokenizedMessage> corpus,
			SerializableIntIntScatterMap seqLengths) {
		Histogram<Integer> seqLengthHisto = new HistogramImpl.IdentityHistogram();
		Iterator<MessageSequenceCorpus.MessageSequence> seqIt = corpus.sequenceIterator();
nextSequence:
		while (seqIt.hasNext()) {
			int pos = 0;
			MessageSequenceCorpus.MessageSequence seq = seqIt.next();
			int sender = -1, recipient = -1;
			while (seq.hasNext()) {
				ProcessedMessage msg = seq.next();
				if (pos == 0) {
					sender = msg.getSenderId();
					recipient = msg.getRecipientIds()[0];
					if (sender == recipient)
						continue nextSequence;
				}
				pos++;
			}

			seqLengthHisto.addValue(pos);
			seqLengths.put(((sender << 16) | recipient), pos);
		}
		System.out.println(seqLengthHisto.toString());
	}

	public static void main(String[] args) throws Exception {
		String trainFile = null, testFile = null;
		if (args.length > 1) {
			trainFile = args[0];
			testFile = args[1];
		}

		Date endDate = cfg.getDateProperty(ExperimentConfiguration.PROP_END_DATE);
		IntervalCalculator ic = new IntervalCalculator(endDate);
		Date startDate = ic.getStartDate();

		List<MessageSequenceCorpus<Long, TokenizedMessage>> corpora;
		if (trainFile == null)
			corpora = buildCorpora(startDate, endDate);
		else
			corpora = loadCorpora(trainFile, testFile);
		SerializableIntIntScatterMap seqLengths = new SerializableIntIntScatterMap();
		printSeqLengthHisto(corpora.get(0), seqLengths);

		fitART(corpora, "original");
		fitSocialTopicModel(corpora, "original", seqLengths);

		for (MessageSequenceCorpus<Long, TokenizedMessage> corpus : corpora)
			corpus.shuffle();
		fitSocialTopicModel(corpora, "shuffled", seqLengths);
	}

}
