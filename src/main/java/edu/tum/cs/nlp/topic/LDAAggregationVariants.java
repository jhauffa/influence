package edu.tum.cs.nlp.topic;

import java.io.File;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.SocialMediaDaoFactory;
import edu.tum.cs.db.loader.MessageLoader;
import edu.tum.cs.db.loader.MessageLoaderFactory;
import edu.tum.cs.math.cluster.DBSCAN1D;
import edu.tum.cs.math.dist.DiscreteDistribution;
import edu.tum.cs.nlp.corpus.Document;
import edu.tum.cs.nlp.corpus.DocumentCorpus;
import edu.tum.cs.nlp.corpus.Index;
import edu.tum.cs.nlp.corpus.ProcessedDocument;
import edu.tum.cs.nlp.corpus.TokenizedMessage;
import edu.tum.cs.nlp.topic.model.LDA;
import edu.tum.cs.util.ExperimentConfiguration;
import edu.tum.cs.util.io.Serializer;

public class LDAAggregationVariants {

	private class AggregateIterator implements Iterator<String> {
		private final Iterator<? extends Document> it;
		private Iterator<String> itWord;

		public AggregateIterator(Iterator<? extends Document> it) {
			this.it = it;
		}

		@Override
		public boolean hasNext() {
			while ((itWord == null) || !itWord.hasNext()) {
				if (!it.hasNext())
					return false;
				itWord = it.next().iterator();
			}
			return itWord.hasNext();
		}

		@Override
		public String next() {
			return itWord.next();
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
	}

	private class AggregateDocument implements Document {
		private final List<? extends Document> docs;
		private final int numTokens;

		public AggregateDocument(List<? extends Document> docs) {
			this.docs = docs;
			int n = 0;
			for (Document d : docs)
				n += d.size();
			this.numTokens = n;
		}

		@Override
		public Iterator<String> iterator() {
			return new AggregateIterator(docs.iterator());
		}

		@Override
		public int size() {
			return numTokens;
		}

		public List<? extends Document> getDocuments() {
			return docs;
		}
	}

	private class MessageComparator implements Comparator<TokenizedMessage> {
		@Override
		public int compare(TokenizedMessage o1, TokenizedMessage o2) {
			return o1.getDate().compareTo(o2.getDate());
		}
	}

	private static final Logger logger = Logger.getLogger(LDAAggregationVariants.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(LDAAggregationVariants.class);

	// LDA parameters
	private static final int numTopics = cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_TOPICS);
	private static final double alpha = 50.0 / numTopics;
	private static final double beta = 0.01;
	private static final int numBurnInIter = 200;
	private static final int numEval = 5;

	private static final int timePeriodLengthDays = 60;

	private static final int[] documentSliceSize = new int[] { 1, 2, 5, 10 };
	private static final int[] minPtsVars = new int[] { 1, 2, 5 };
	private static final double[] epsilonTimeFactorVars = new double[] { 0.25, 0.5, 0.75, 1.5 };

	// 5 minutes, 15 minutes, 30 minutes, 1 hour, 3 hours, 6 hours, 12 hours, 24 hours
	private static final long[] timeSteps = new long[] {
		5 * 60 * 1000L, 15 * 60 * 1000L, 30 * 60 * 1000L, 60 * 60 * 1000L, 3 * 60 * 60 * 1000L, 6 * 60 * 60 * 1000L,
		12 * 60 * 60 * 1000L, 24 * 60 * 60 * 1000L
	};

	private final Map<Long, List<TokenizedMessage>> userId2Msg;
	private final Set<Long> testUserIds;
	private final Date startDate, endDate;
	private final boolean aggregateTestCorpus, splitBySender, splitByRecipient;
	private final Index<String> bagOfWords;
	private final AtomicInteger numLDAModels = new AtomicInteger();

	public LDAAggregationVariants(Date startDate, Date endDate, Map<Long, List<TokenizedMessage>> userId2Msg,
			boolean aggregateTestCorpus, boolean splitBySender, boolean splitByRecipient) {
		this.startDate = startDate;
		this.endDate = endDate;
		this.userId2Msg = userId2Msg;
		this.aggregateTestCorpus = aggregateTestCorpus;
		this.splitBySender = splitBySender;
		this.splitByRecipient = splitByRecipient;

		logger.info("using " + (aggregateTestCorpus ? "" : "non-") + "aggregated test corpus");
		logger.info("splitting by: " + (splitBySender ? "sender " : "") + (splitByRecipient ? "recipient" : ""));
		logger.info("messages between " + startDate + " and " + endDate);
		logger.info("LDA with " + numTopics + " topics");

		// add 10% of the users to the test set
		this.testUserIds = new HashSet<Long>();
		Iterator<Long> it = userId2Msg.keySet().iterator();
		for (int i = 0; i < (userId2Msg.keySet().size() / 10); i++)
			testUserIds.add(it.next());

		// sort messages of each user by date
		for (List<TokenizedMessage> msgs : userId2Msg.values())
			Collections.sort(msgs, new MessageComparator());

		// build word index
		List<Document> allDocs = new ArrayList<Document>(userId2Msg.keySet().size());
		for (Map.Entry<Long, List<TokenizedMessage>> e : userId2Msg.entrySet())
			allDocs.add(new AggregateDocument(e.getValue()));
		bagOfWords = new Index<String>();
		new DocumentCorpus(bagOfWords, allDocs);
	}

	public void runExperiments() throws Exception {
		int cores = cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_THREADS,
				Runtime.getRuntime().availableProcessors());
		ThreadPoolExecutor executor = new ThreadPoolExecutor(cores, cores, 60, TimeUnit.SECONDS,
				new LinkedBlockingQueue<Runnable>());
		List<Future<?>> tasks = new ArrayList<Future<?>>();

		// fixed document slices
		for (int numOfMsgsPerDocument : documentSliceSize)
			tasks.add(executor.submit(new FixedDocumentSlicesVariant(numOfMsgsPerDocument)));

		// clustering with DBSCAN
		for (int minPts : minPtsVars) {
			for (double epsilonTimeFactor : epsilonTimeFactorVars)
				tasks.add(executor.submit(new DBSCANVariant(minPts, epsilonTimeFactor)));
		}

		// time bins
		for (long timeStep : timeSteps)
			tasks.add(executor.submit(new TimeBinsVariant(timeStep)));

		// all messages of a user (pair)
		if (splitBySender || splitByRecipient)
			tasks.add(executor.submit(new OneDocumentPerUserVariant()));

		logger.info(tasks.size() + " variants to be processed");
		executor.shutdown();

		for (Future<?> task : tasks) {
			try {
				task.get();
			} catch (ExecutionException ex) {
				logger.log(Level.SEVERE, "error during execution of task, results will be incomplete", ex);
			}
		}
		logger.info("done");
	}

	private abstract class AggregationVariant implements Runnable {
		protected abstract String getVariantId();

		/**
		 * @param msgs a temporally sorted list of messages
		 */
		protected abstract Collection<AggregateDocument> performAggregation(List<TokenizedMessage> msgs);

		private Collection<List<TokenizedMessage>> splitBySender() {
			if (!splitBySender) {
				// Messages are already separated by sender, because this is usually the most convenient representation.
				// Aggregate them, but keep training and test data separate, and keep the resulting lists in temporal
				// order.
				List<TokenizedMessage> mergedTrainMsgs = new ArrayList<TokenizedMessage>();
				List<TokenizedMessage> mergedTestMsgs = new ArrayList<TokenizedMessage>();
				for (Map.Entry<Long, List<TokenizedMessage>> e : userId2Msg.entrySet()) {
					if (testUserIds.contains(e.getKey()))
						mergedTestMsgs.addAll(e.getValue());
					else
						mergedTrainMsgs.addAll(e.getValue());
				}
				Collections.sort(mergedTrainMsgs, new MessageComparator());
				Collections.sort(mergedTestMsgs, new MessageComparator());

				Collection<List<TokenizedMessage>> msgs = new ArrayList<List<TokenizedMessage>>(2);
				msgs.add(mergedTrainMsgs);
				msgs.add(mergedTestMsgs);
				return msgs;
			}

			return userId2Msg.values();
		}

		private Collection<List<TokenizedMessage>> splitByRecipient(List<TokenizedMessage> allMsgs) {
			if (!splitByRecipient) {
				Collection<List<TokenizedMessage>> msgs = new ArrayList<List<TokenizedMessage>>(1);
				msgs.add(allMsgs);
				return msgs;
			}

			Map<Long, List<TokenizedMessage>> recip2Msg = new HashMap<Long, List<TokenizedMessage>>();
			for (TokenizedMessage msg : allMsgs) {
				for (Long recipId : msg.getRecipients()) {
					List<TokenizedMessage> msgs = recip2Msg.get(recipId);
					if (msgs == null) {
						msgs = new ArrayList<TokenizedMessage>();
						recip2Msg.put(recipId, msgs);
					}
					msgs.add(msg);
				}
			}
			return recip2Msg.values();
		}

		private void addToCorpus(AggregateDocument doc, boolean isTestData, List<Document> trainDocs,
				List<Document> testDocs) {
			if (isTestData) {
				if (aggregateTestCorpus) {
					testDocs.add(doc);
				} else {
					testDocs.addAll(doc.getDocuments());
				}
			} else {
				trainDocs.add(doc);
			}
		}

		private void fitAndEvalTopicModel(List<Document> trainDocs, List<Document> testDocs, String variantId) {
			logger.log(Level.INFO, "now processing " + variantId + ", total processed: " + numLDAModels.get());

			DocumentCorpus trainCorpus = new DocumentCorpus(bagOfWords, trainDocs);
			int numDocuments = trainCorpus.size();
			DescriptiveStatistics tokenStats = new DescriptiveStatistics();
			for (ProcessedDocument doc : trainCorpus)
				tokenStats.addValue(doc.size());
			DocumentCorpus testCorpus = new DocumentCorpus(bagOfWords, testDocs);
			int numTestTokens = testCorpus.countTokens();

			DescriptiveStatistics llStats = new DescriptiveStatistics();
			DescriptiveStatistics perplexityStats = new DescriptiveStatistics();
			for (int i = 0; i < numEval; i++) {
				LDA lda = new LDA(bagOfWords.getNumUniqueElements(), numTopics, alpha, beta);
				lda.setBurnIn(numBurnInIter);
				lda.train(trainCorpus);

				double llHeldOut = lda.queryLogLikelihoodLR(testCorpus, LDA.ResampleMode.LOG);
				llStats.addValue(llHeldOut);
				perplexityStats.addValue(DiscreteDistribution.perplexity(llHeldOut, numTestTokens));
			}

			System.out.println(variantId + ";" + numDocuments + ";" + tokenStats.getMin() + ";" + tokenStats.getMax() +
						";" + tokenStats.getMean() + ";" + tokenStats.getStandardDeviation() + ";" +
						tokenStats.getPercentile(50) + ";" + llStats.getMean() + ";" + llStats.getStandardDeviation() +
						";" + perplexityStats.getMean() + ";" + perplexityStats.getStandardDeviation());
			numLDAModels.incrementAndGet();
		}

		@Override
		public void run() {
			List<Document> trainDocs = new ArrayList<Document>();
			List<Document> testDocs = new ArrayList<Document>();
			for (List<TokenizedMessage> senderMsgs : splitBySender()) {
				boolean isTestData = testUserIds.contains(senderMsgs.get(0).getSender());
				for (List<TokenizedMessage> msgs : splitByRecipient(senderMsgs)) {
					for (AggregateDocument doc : performAggregation(msgs))
						addToCorpus(doc, isTestData, trainDocs, testDocs);
				}
			}
			fitAndEvalTopicModel(trainDocs, testDocs, getVariantId());
		}
	}

	private class OneDocumentPerUserVariant extends AggregationVariant {
		@Override
		protected String getVariantId() {
			return "one document per sender" + (splitByRecipient ? "-recipient pair" : "");
		}

		@Override
		protected Collection<AggregateDocument> performAggregation(List<TokenizedMessage> msgs) {
			Collection<AggregateDocument> docs = new ArrayList<AggregateDocument>(1);
			docs.add(new AggregateDocument(msgs));
			return docs;
		}
	}

	private class TimeBinsVariant extends AggregationVariant {
		private final long timeStep;

		public TimeBinsVariant(long timeStep) {
			this.timeStep = timeStep;
		}

		@Override
		protected String getVariantId() {
			return "fixed duration (" + (int) (timeStep / 1000) + " s)";
		}

		@Override
		protected Collection<AggregateDocument> performAggregation(List<TokenizedMessage> msgs) {
			int numSlices = (int) (Math.ceil((endDate.getTime() - startDate.getTime()) / (double) timeStep));
			int arrayCapacity = Math.max(msgs.size() / numSlices, 1);
			Collection<AggregateDocument> docs = new ArrayList<AggregateDocument>(numSlices);

			Date currentDate = new Date(startDate.getTime() + timeStep);
			List<TokenizedMessage> currentSlice = new ArrayList<TokenizedMessage>(arrayCapacity);
			for (TokenizedMessage msg : msgs) {
				// search for the next fitting time slot
				while (currentDate.before(msg.getDate())) {
					currentDate = new Date(currentDate.getTime() + timeStep);
					if (!currentSlice.isEmpty()) {
						docs.add(new AggregateDocument(currentSlice));
						currentSlice = new ArrayList<TokenizedMessage>(arrayCapacity);
					}
				}
				currentSlice.add(msg);
			}
			if (!currentSlice.isEmpty())
				docs.add(new AggregateDocument(currentSlice));
			return docs;
		}
	}

	private class DBSCANVariant extends AggregationVariant {
		private final int minPts;
		private final double epsilonTimeFactor;

		public DBSCANVariant(int minPts, double epsilonTimeFactor) {
			this.minPts = minPts;
			this.epsilonTimeFactor = epsilonTimeFactor;
		}

		@Override
		protected String getVariantId() {
			return "DBSCAN (" + minPts + ", " + epsilonTimeFactor + ")";
		}

		@Override
		protected Collection<AggregateDocument> performAggregation(List<TokenizedMessage> msgs) {
			Collection<AggregateDocument> docs = new ArrayList<AggregateDocument>();
			if (msgs.size() >= 2 * minPts) {
				long userObservedMillis = msgs.get(msgs.size() - 1).getDate().getTime() -
						msgs.get(0).getDate().getTime();
				double avgMinutesPerMsg = (userObservedMillis / 1000.0 / 60.0) / msgs.size();
				double epsilon = avgMinutesPerMsg * epsilonTimeFactor;

				for (List<TokenizedMessage> cluster : DBSCAN1D.clusterMessages(msgs, epsilon, minPts))
					docs.add(new AggregateDocument(cluster));
			} else if (!msgs.isEmpty()) {
				docs.add(new AggregateDocument(msgs));
			}
			return docs;
		}
	}

	private class FixedDocumentSlicesVariant extends AggregationVariant {
		private final int numOfMsgsPerDocument;

		public FixedDocumentSlicesVariant(int numOfMsgsPerDocument) {
			this.numOfMsgsPerDocument = numOfMsgsPerDocument;
		}

		@Override
		protected String getVariantId() {
			return "fixed size (" + numOfMsgsPerDocument + ")";
		}

		@Override
		protected Collection<AggregateDocument> performAggregation(List<TokenizedMessage> msgs) {
			Collection<AggregateDocument> docs = new ArrayList<AggregateDocument>(
					(msgs.size() + numOfMsgsPerDocument - 1) / numOfMsgsPerDocument);
			List<TokenizedMessage> currentSlice = new ArrayList<TokenizedMessage>();
			for (TokenizedMessage msg : msgs) {
				currentSlice.add(msg);
				if (currentSlice.size() == numOfMsgsPerDocument) {
					docs.add(new AggregateDocument(currentSlice));
					currentSlice = new ArrayList<TokenizedMessage>(numOfMsgsPerDocument);
				}
			}
			if (!currentSlice.isEmpty())
				docs.add(new AggregateDocument(currentSlice));
			return docs;
		}
	}

	public static void main(String[] args) throws Exception {
		Date endDate = cfg.getDateProperty(ExperimentConfiguration.PROP_END_DATE);
		GregorianCalendar startCal = new GregorianCalendar();
		startCal.setTimeInMillis(endDate.getTime());
		startCal.add(Calendar.DAY_OF_MONTH, -timePeriodLengthDays);
		Date startDate = startCal.getTime();

		boolean aggregateTestCorpus = false;
		boolean splitBySender = false;
		boolean splitByRecipient = false;
		Map<Long, List<TokenizedMessage>> userId2Msg = null;
		for (String arg : args) {
			if (arg.equals("test"))
				aggregateTestCorpus = true;
			else if (arg.equals("sender"))
				splitBySender = true;
			else if (arg.equals("recipient"))
				splitByRecipient = true;
			else
				userId2Msg = Serializer.loadObjectFromFile(new File(arg));
		}

		if (userId2Msg == null) {
			SocialMediaDao dao = SocialMediaDaoFactory.createDao();
			Set<Long> userIds = dao.getUserIds(false);
			MessageLoader loader = MessageLoaderFactory.createMessageLoader(dao, userIds);

			userId2Msg = new HashMap<Long, List<TokenizedMessage>>();
			for (long userId : userIds) {
				List<TokenizedMessage> msg = loader.loadMessages(userId, startDate, endDate);
				if (!msg.isEmpty())
					userId2Msg.put(userId, msg);
			}
			Serializer.saveObjectToFile(userId2Msg, new File("messages.ser.gz"));
		}

		LDAAggregationVariants ldaVariants = new LDAAggregationVariants(startDate, endDate, userId2Msg,
				aggregateTestCorpus, splitBySender, splitByRecipient);
		ldaVariants.runExperiments();
	}

}
