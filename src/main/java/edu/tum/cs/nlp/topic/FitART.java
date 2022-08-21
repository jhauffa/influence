package edu.tum.cs.nlp.topic;

import java.io.File;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.SocialMediaDaoFactory;
import edu.tum.cs.db.loader.MessageLoader;
import edu.tum.cs.db.loader.MessageLoaderFactory;
import edu.tum.cs.nlp.corpus.Corpus;
import edu.tum.cs.nlp.corpus.Index;
import edu.tum.cs.nlp.corpus.MessageCorpus;
import edu.tum.cs.nlp.corpus.ProcessedMessage;
import edu.tum.cs.nlp.corpus.TokenizedMessage;
import edu.tum.cs.nlp.topic.model.ART;
import edu.tum.cs.nlp.topic.model.OnlineTopicModel;
import edu.tum.cs.nlp.topic.model.ParallelART;
import edu.tum.cs.time.IntervalCalculator;
import edu.tum.cs.util.ExperimentConfiguration;
import edu.tum.cs.util.io.Serializer;

public class FitART {

	private static final Logger logger = Logger.getLogger(FitART.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(FitART.class);

	public static final String bowFileName = "bagOfWords.ser.gz";
	public static final String bopFileName = "bagOfPersons.ser.gz";

	private static final int numTopics = cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_TOPICS);
	private static final ART.AlphaType alphaType = ART.AlphaType.DATA_DRIVEN_NON_UNIFORM;
	private static final double alpha = 50.0 / numTopics;
	private static final double beta = 0.01;
	private static final int numIter = 2000;

	private static final double[] weights = { 1.0 };
	private static final boolean normalizeWeights = true;

	private static Corpus<ProcessedMessage> loadCorpusSlice(int i, Date[] dates, MessageLoader loader,
			Set<Long> userIdSet, Index<String> bagOfWords, Index<Long> bagOfPersons) throws Exception {
		logger.info("Loading messages from " + dates[i] + " to " + dates[i + 1]);
		List<TokenizedMessage> messages = new ArrayList<TokenizedMessage>();
		int idx = 0;
		for (Long userId : userIdSet) {
			messages.addAll(loader.loadMessages(userId, dates[i], dates[i + 1]));
			if (++idx % 1000 == 0)
				logger.info("Done with " + idx);
		}
		return new MessageCorpus<Long, TokenizedMessage>(bagOfWords, bagOfPersons, messages);
	}

	public static void main(String[] args) throws Exception {
		Date endDate = cfg.getDateProperty(ExperimentConfiguration.PROP_END_DATE);
		IntervalCalculator ic = new IntervalCalculator(endDate);
		Date[] dates;
		boolean useOnlineTopicModel = cfg.getBooleanProperty(ExperimentConfiguration.PROP_USE_ONLINE_TOPIC_MODEL);
		if (useOnlineTopicModel) {
			// compute dates for the four different time slices
			dates = new Date[5];
			int idx = dates.length - 1;
			dates[idx--] = endDate;
			for (Date centerDate : ic.getCenterDates())
				dates[idx--] = centerDate;
			dates[idx] = ic.getStartDate();
		} else {
			dates = new Date[] { ic.getStartDate(), endDate };
		}
		String modelPath = cfg.getProperty(ExperimentConfiguration.PROP_TOPIC_MODEL_PATH);

		int resumeAtIter = 0, resumeAtCorpus = 0;
		if (args.length > 0) {
			resumeAtIter = Integer.parseInt(args[0]);
			if (args.length > 1)
				resumeAtCorpus = Integer.parseInt(args[1]);
		}

		// load messages
		Index<String> bagOfWords;
		Index<Long> bagOfPersons;
		if ((resumeAtIter == 0) && (resumeAtCorpus == 0)) {
			bagOfWords = new Index<String>();
			bagOfPersons = new Index<Long>();
		} else {
			bagOfWords = Serializer.loadObjectFromFile(new File(modelPath, bowFileName));
			bagOfWords.setReadOnly();
			bagOfPersons = Serializer.loadObjectFromFile(new File(modelPath, bopFileName));
			bagOfPersons.setReadOnly();
		}

		SocialMediaDao dao = SocialMediaDaoFactory.createDao();
		// quality of model directly depends on amount of data, so use all relevant messages within the chosen interval
		Set<Long> userIdSet = dao.getUserIds(false);
		MessageLoader loader = MessageLoaderFactory.createMessageLoader(dao, userIdSet);
		// load all messages to obtain bags of words/persons that are valid for all slices; discard all but the first
		// slice
		Corpus<ProcessedMessage> curCorpus = null;
		int numCorpusSlices = dates.length - 1;
		for (int i = numCorpusSlices - 1; i >= 0; i--) {
			Corpus<ProcessedMessage> corpus = loadCorpusSlice(i, dates, loader, userIdSet, bagOfWords, bagOfPersons);
			if (i == resumeAtCorpus)
				curCorpus = corpus;
		}
		logger.info("Finished loading messages - processed " + userIdSet.size() + " users");

		int numWords = bagOfWords.getNumUniqueElements();
		int numPersons = bagOfPersons.getNumUniqueElements();
		logger.info(numWords + " unique words and " + numPersons + " unique persons");
		Serializer.saveObjectToFile(bagOfWords, new File(modelPath, bowFileName));
		Serializer.saveObjectToFile(bagOfPersons, new File(modelPath, bopFileName));
		if (!useOnlineTopicModel)
			bagOfPersons = null;

		Runtime runtime = Runtime.getRuntime();
		runtime.gc();
		logger.info("Currently used memory " + ((runtime.totalMemory() - runtime.freeMemory()) / 1024L / 1024L) + "MB");

		// start ART parameter fitting
		OnlineTopicModel<ProcessedMessage, ART> onlineModels = null;
		ART art = null;
		if (resumeAtCorpus > 0) {
			List<ART> models = new ArrayList<ART>(resumeAtCorpus);
			List<Corpus<ProcessedMessage>> corpora = new ArrayList<Corpus<ProcessedMessage>>(resumeAtCorpus);
			for (int i = 0; i < resumeAtCorpus; i++) {
				ART model = Serializer.loadObjectFromFile(new File(modelPath, "onlineART-model-" + (i+1) + ".ser.gz"));
				models.add(model);
				corpora.add(loadCorpusSlice(i, dates, loader, userIdSet, bagOfWords, bagOfPersons));
			}
			onlineModels = new OnlineTopicModel<ProcessedMessage, ART>(models, corpora, weights, normalizeWeights);
		} else {
			if (resumeAtIter > 0) {
				art = Serializer.loadObjectFromFile(new File(modelPath, "topic-model-state.ser.gz"));
			} else {
				ParallelART part = new ParallelART(numWords, numPersons, numTopics, alpha, beta);
				part.setBurnIn(numIter);
				part.setBurnInReportRate(50);
				part.setAlphaSamplingType(alphaType);
				part.setMemorySavingStartup(true);
				part.setUseExternalMemoryCorpus(true);
				art = part;
			}
			onlineModels = new OnlineTopicModel<ProcessedMessage, ART>(art, weights, normalizeWeights);
		}

		if (useOnlineTopicModel) {
			for (int i = resumeAtCorpus; i < numCorpusSlices; i++) {
				if (i > resumeAtCorpus)
					curCorpus = loadCorpusSlice(i, dates, loader, userIdSet, bagOfWords, bagOfPersons);
				logger.info("Updating model with corpus " + (i + 1) + ", " + curCorpus.size() + " messages");
				if (i > 0)
					onlineModels.getCurrentModel().setBurnIn(numIter / 4);
				onlineModels.update(curCorpus);

				ART curModel = onlineModels.getCurrentModel();
				Serializer.saveObjectToFile(curModel, new File(modelPath, "onlineART-model-" + (i+1) + ".ser.gz"));
				TopicWordDistribution.saveTopicsCsv(curModel.getTopicWordDistr(), bagOfWords,
						new File(modelPath, "topics-model-" + (i+1) + ".csv"));
				TopicWordDistribution.saveTopicsLatex(curModel.getTopicWordDistr(), bagOfWords,
						new File(modelPath, "topics-model-" + (i+1) + ".tex"));
			}

			Serializer.saveObjectToFile(onlineModels, new File(modelPath, "onlineART-model-full.ser.gz"));
			onlineModels.getCurrentModel().getIntermediateStateFile().delete();
		} else {
			logger.info("Starting ART for " + curCorpus.size() + " messages");
			if (resumeAtIter > 0)
				art.resumeTraining(curCorpus, art.getBurnIn() - resumeAtIter);
			else
				art.train(curCorpus);

			TopicWordDistribution.saveTopicsCsv(art.getTopicWordDistr(), bagOfWords, new File(modelPath, "topics.csv"));
			TopicWordDistribution.saveTopicsLatex(art.getTopicWordDistr(), bagOfWords,
					new File(modelPath, "topics.tex"));

			Serializer.saveObjectToFile(art, new File(modelPath, "art-model-final.ser.gz"));
			art.getIntermediateStateFile().delete();
		}

		logger.info("Done.");
	}

}
