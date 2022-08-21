package edu.tum.cs.influence.micro;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import org.apache.commons.collections15.Predicate;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.SocialMediaDaoFactory;
import edu.tum.cs.db.loader.MessageLoader;
import edu.tum.cs.db.loader.MessageLoaderFactory;
import edu.tum.cs.graph.GraphCreator;
import edu.tum.cs.math.dist.DiscreteDistribution;
import edu.tum.cs.nlp.corpus.Index;
import edu.tum.cs.nlp.corpus.MessageCorpus;
import edu.tum.cs.nlp.corpus.TokenizedMessage;
import edu.tum.cs.nlp.topic.FitART;
import edu.tum.cs.nlp.topic.model.ART;
import edu.tum.cs.time.IntervalCalculator;
import edu.tum.cs.util.ExperimentConfiguration;
import edu.tum.cs.util.arrays.Object2DArray;
import edu.tum.cs.util.io.Serializer;
import edu.uci.ics.jung.algorithms.filters.VertexPredicateFilter;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;

public class InfluenceMeasurement {

	private static final Logger logger = Logger.getLogger(InfluenceMeasurement.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(InfluenceMeasurement.class);

	public static final String PROP_MUTUAL_ONLY = "mutualNetworksOnly";

	private static final String modelPath = cfg.getProperty(ExperimentConfiguration.PROP_TOPIC_MODEL_PATH);
	private static final String outPath = cfg.getLocalProperty("outputPath", ".");
	private static final boolean testMutualNetworksOnly = cfg.getLocalBooleanProperty(PROP_MUTUAL_ONLY);
	private static final boolean generateExplanations = cfg.getLocalBooleanProperty("generateExplanations", false);
	private static final double noiseTopicThreshold = cfg.getLocalDoubleProperty("noiseTopicThreshold");

	private static class ExperimentResult {
		private final TimeSeriesVariant dataVariant;
		private final ModelVariant modelVariant;
		private final String graphFileName;
		private final int numOriginalEdges;
		private final double earlyRejectRate;
		private final double networkCoverage;
		private final double fakeNetworkCoverage;

		public ExperimentResult(TimeSeriesVariant dataVariant, ModelVariant modelVariant, String graphFileName,
				int numOriginalEdges, double earlyRejectRate, double networkCoverage, double fakeNetworkCoverage) {
			this.dataVariant = dataVariant;
			this.modelVariant = modelVariant;
			this.graphFileName = graphFileName;
			this.numOriginalEdges = numOriginalEdges;
			this.earlyRejectRate = earlyRejectRate;
			this.networkCoverage = networkCoverage;
			this.fakeNetworkCoverage = fakeNetworkCoverage;
		}

		@Override
		public String toString() {
			return dataVariant + ";" + modelVariant + ";" + graphFileName + ";" + numOriginalEdges + ";" +
					earlyRejectRate + ";" + networkCoverage + ";" + fakeNetworkCoverage;
		}

		public static String getCsvHeader() {
			return TimeSeriesVariant.getCsvHeader() + ";" + ModelVariant.getCsvHeader() +
					";graph;numOriginalEdges;earlyRejectRate;networkCoverage;fakeNetworkCoverage";
		}
	}

	private static class ProgressReporter extends Thread {
		private static final long reportInterval = 5 * 60 * 1000L;
		private final TemporalInfluenceModel<?> model;

		public ProgressReporter(TemporalInfluenceModel<?> model) {
			this.model = model;
		}

		@Override
		public void run() {
			Runtime runtime = Runtime.getRuntime();
			double prevProgress = 0.0;
			while (!interrupted()) {
				try {
					Thread.sleep(reportInterval);
					long sumNum = 0, sumDenom = 0;
					for (Map.Entry<String, long[]> e : model.reportProgress().entrySet()) {
						logger.info(e.getValue()[0] + " of " + e.getValue()[1] + " " + e.getKey() + " (" +
								(((double) e.getValue()[0] / e.getValue()[1]) * 100) + "%)");
						sumNum += e.getValue()[0];
						sumDenom += e.getValue()[1];
					}
					double progress = ((double) sumNum / sumDenom) * 100;
					logger.info("overall: " + sumNum + " of " + sumDenom + " (" + progress + "%)");
					double rate = (progress - prevProgress) / (reportInterval / (60 * 1000L));
					if (rate > 0.0)
						logger.info("estimated remaining time: " + (int) ((100.0 - progress) / rate) + " m");
					else
						logger.info("remaining time unknown");
					prevProgress = progress;
					runtime.gc();
					logger.info("free memory: " + (runtime.freeMemory() / (1024L * 1024L)) + " MB, used memory: " +
						((runtime.totalMemory() - runtime.freeMemory()) / (1024L * 1024L)) + " MB");
				} catch (InterruptedException ex) {
					break;
				}
			}
		}
	}

	private final Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> networks =
			new HashMap<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>>();
	private final Map<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> fakeNetworks =
			new HashMap<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>>();
	private Map<TimeSeriesVariant, Map<Long, TimeSeriesData>> series;
	private double[] topicPriorParam;

	private Map<TimeSeriesVariant, Map<Long, TimeSeriesData>> buildTimeSeries(SocialMediaDao mediaDao,
			Set<Long> userIds, ART art, Date startDate, Date endDate, int[] timePeriods,
			DirectedSparseGraph<Long, Integer> networkGraph) throws Exception {
		// Load word and person indices
		Index<String> bagOfWords = Serializer.loadObjectFromFile(new File(modelPath, FitART.bowFileName));
		bagOfWords.setReadOnly();
		Index<Long> bagOfPersons = Serializer.loadObjectFromFile(new File(modelPath, FitART.bopFileName));
		bagOfPersons.setReadOnly();

		// Load users' messages
		MessageLoader loader = MessageLoaderFactory.createMessageLoader(mediaDao, userIds);
		GraphCreator.applyClassifier(loader, networkGraph);
		loader.cacheMessages(startDate, endDate);

		Map<TimeSeriesVariant, Map<Long, TimeSeriesData>> series =
				new HashMap<TimeSeriesVariant, Map<Long, TimeSeriesData>>();
		for (int timePeriodLength : timePeriods) {
			// Generate time slices
			TimeSeriesVariant variant = new TimeSeriesVariant(timePeriodLength);
			List<Date[]> sliceDates = new ArrayList<Date[]>();
			GregorianCalendar cal = new GregorianCalendar();
			cal.setTime(startDate);
			Date curStartDate = cal.getTime();
			cal.add(Calendar.DAY_OF_MONTH, timePeriodLength);
			Date curEndDate = cal.getTime();
			while (!curEndDate.after(endDate)) {
				sliceDates.add(new Date[] { curStartDate, curEndDate });
				curStartDate = curEndDate;
				cal.add(Calendar.DAY_OF_MONTH, timePeriodLength);
				curEndDate = cal.getTime();
			}

			// Generate time series data for each slice by querying the ART model
			Map<Long, TimeSeriesData> seriesData = new HashMap<Long, TimeSeriesData>(userIds.size());
			series.put(variant, seriesData);
			for (Long userId : userIds)
				seriesData.put(userId, new TimeSeriesData(sliceDates.size()));

			int numSlices = 0;
			logger.info("Generating slices of length " + timePeriodLength);
			for (Date[] timeSlice : sliceDates) {
				logger.info("Slice from " + timeSlice[0] + " to " + timeSlice[1] + " starting");

				// Build corpus for ART query
				List<TokenizedMessage> messages = new ArrayList<TokenizedMessage>();
				Set<Long> activeUsers = new HashSet<Long>();
				for (long userId : userIds) {
					List<TokenizedMessage> docs = loader.loadMessages(userId, timeSlice[0], timeSlice[1]);
					if (!docs.isEmpty())
						activeUsers.add(userId);
					for (TokenizedMessage doc : docs) {
						// assign all documents to sender loop edge
						doc.getRecipients().clear();
						doc.getRecipients().add(doc.getSender());
						messages.add(doc);
					}
				}
				MessageCorpus<Long, TokenizedMessage> corpus = new MessageCorpus<Long, TokenizedMessage>(bagOfWords,
						bagOfPersons, messages);

				// Perform query and store result
				ART.QueryResult state = art.query(corpus);
				Object2DArray<double[]> thetas = state.estimateTheta();
				for (Long userId : userIds) {
					double[] theta = null;
					if (activeUsers.contains(userId)) {
						int id = bagOfPersons.getElementId(userId);
						theta = thetas.get(id, id);
					}
					seriesData.get(userId).addValue(theta);
				}
				logger.info("Done with " + (++numSlices) + " of " + sliceDates.size());
			}
		}
		return series;
	}

	private Map<TimeSeriesVariant, Map<Long, TimeSeriesData>> filterTimeSeries(
			Map<TimeSeriesVariant, Map<Long, TimeSeriesData>> series) {
		int[] topicMap = new int[topicPriorParam.length];
		int numRemainingTopics = 0;
		for (int i = 0; i < topicPriorParam.length; i++) {
			if (topicPriorParam[i] < noiseTopicThreshold)
				topicMap[numRemainingTopics++] = i;
		}

		Map<TimeSeriesVariant, Map<Long, TimeSeriesData>> filteredSeries =
				new HashMap<TimeSeriesVariant, Map<Long, TimeSeriesData>>();
		for (Map.Entry<TimeSeriesVariant, Map<Long, TimeSeriesData>> e : series.entrySet()) {
			Map<Long, TimeSeriesData> variantSeries = new HashMap<Long, TimeSeriesData>();
			filteredSeries.put(e.getKey(), variantSeries);
			for (Map.Entry<Long, TimeSeriesData> ee : e.getValue().entrySet()) {
				int n = ee.getValue().getLength();
				TimeSeriesData filteredData = new TimeSeriesData(n);
				for (int i = 0; i < n; i++) {
					double[] theta = ee.getValue().getValue(i);
					if (theta != null) {
						double[] thetaTf = new double[numRemainingTopics];
						for (int j = 0; j < numRemainingTopics; j++)
							thetaTf[j] = theta[topicMap[j]];
						DiscreteDistribution.normalize(thetaTf);
						theta = thetaTf;
					}
					filteredData.addValue(theta);
				}
				variantSeries.put(ee.getKey(), filteredData);
			}
		}
		return filteredSeries;
	}

	public void analyzeData(File f) throws Exception {
		// write missing data statistics to file
		PrintWriter w = new PrintWriter(new FileWriter(f));
		try {
			w.println(TimeSeriesVariant.getCsvHeader() + ";mean;std.dev.");
			for (Map.Entry<TimeSeriesVariant, Map<Long, TimeSeriesData>> e : series.entrySet()) {
				SummaryStatistics stats = new SummaryStatistics();
				for (TimeSeriesData s : e.getValue().values())
					stats.addValue((double) s.getNumMissing() / s.getLength());
				w.println(e.getKey() + ";" + stats.getMean() + ";" + stats.getStandardDeviation());
			}
		} finally {
			w.close();
		}

		// print number of nodes/edges and statistics of in-degree distribution of selected networks
		for (Map.Entry<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> e : networks.entrySet()) {
			logger.info(e.getKey() + ": " + e.getValue().getVertexCount() + " V " + e.getValue().getEdgeCount() + " E");

			SummaryStatistics stats = new SummaryStatistics();
			for (Long v : e.getValue().getVertices())
				stats.addValue(e.getValue().inDegree(v));
			logger.info("in-degree distribution of " + e.getKey() + ": " + stats);
		}
	}

	public int[] loadNetworksAndActions(String followNetworkFileName, String commNetworkFileName,
			String seriesFileName) throws Exception {
		SocialMediaDao mediaDao = SocialMediaDaoFactory.createDao();
		final Set<Long> userIds = mediaDao.getUserIds(true);
		VertexPredicateFilter<Long, Integer> filter =
				new VertexPredicateFilter<Long, Integer>(new Predicate<Long>() {
			@Override
			public boolean evaluate(Long id) {
				return userIds.contains(id);
			}
		});

		DirectedSparseGraph<Long, Integer> followerNetwork;
		if (followNetworkFileName != null) {
			followerNetwork = Serializer.loadObjectFromFile(new File(followNetworkFileName));
			followerNetwork = (DirectedSparseGraph<Long, Integer>) filter.transform(followerNetwork);
		} else {
			followerNetwork = GraphCreator.createFollowerGraph(mediaDao, userIds);
			Serializer.saveObjectToFile(followerNetwork, new File(outPath, "follow-graph.ser.gz"));
		}

		// load ART model
		ART art = Serializer.loadObjectFromFile(new File(modelPath, "art-model-final.ser.gz"));
		art.configureForQuerying();
		art.setQuerySamples(art.getBurnIn());	// hopefully a good compromise
		topicPriorParam = art.getAlpha();

		// load messages, compute topic distributions for time series and build communication graph
		int[] timePeriods = cfg.getIntListProperty(ExperimentConfiguration.PROP_TIME_PERIODS);
		DirectedSparseGraph<Long, Integer> communicationNetwork;
		if ((commNetworkFileName != null) && (seriesFileName != null)) {
			communicationNetwork = Serializer.loadObjectFromFile(new File(commNetworkFileName));
			communicationNetwork = (DirectedSparseGraph<Long, Integer>) filter.transform(communicationNetwork);
			series = Serializer.loadObjectFromFile(new File(seriesFileName));
			for (Map<Long, TimeSeriesData> variantSeries : series.values()) {
				Iterator<Long> it = variantSeries.keySet().iterator();
				while (it.hasNext()) {
					Long e = it.next();
					if (!userIds.contains(e))
						it.remove();
				}
			}
		} else {
			communicationNetwork = GraphCreator.createNodeGraph(userIds);
			Date endDate = cfg.getDateProperty(ExperimentConfiguration.PROP_END_DATE);
			IntervalCalculator ic = new IntervalCalculator(endDate);
			series = buildTimeSeries(mediaDao, userIds, art, ic.getStartDate(), endDate, timePeriods,
					communicationNetwork);
			Serializer.saveObjectToFile(communicationNetwork, new File(outPath, "comm-graph.ser.gz"));
			Serializer.saveObjectToFile(series, new File(outPath, "series.ser.gz"));
		}
		if (noiseTopicThreshold > 0.0)
			series = filterTimeSeries(series);

		if (!testMutualNetworksOnly) {
			networks.put(ModelVariant.NetworkType.FOLLOWING, followerNetwork);
			networks.put(ModelVariant.NetworkType.REPLY, communicationNetwork);
		}
		networks.put(ModelVariant.NetworkType.MUTUAL_FOLLOWING, GraphCreator.keepBidirectionalEdges(followerNetwork));
		networks.put(ModelVariant.NetworkType.MUTUAL_REPLY, GraphCreator.keepBidirectionalEdges(communicationNetwork));

		// construct fake network graphs
		for (Map.Entry<ModelVariant.NetworkType, DirectedSparseGraph<Long, Integer>> e : networks.entrySet())
			fakeNetworks.put(e.getKey(), GraphCreator.sampleComplementGraph(e.getValue(), e.getValue().getEdgeCount()));

		return timePeriods;
	}

	private void writeExplanations(File f, WeightedDirectedGraph<Long, Explanation> exp)
			throws IOException {
		PrintWriter w = new PrintWriter(f);
		try {
			for (WeightedEdge<Explanation> e : exp.getEdges()) {
				Pair<Long> v = exp.getEndpoints(e);
				w.println(v.getFirst() + " -> " + v.getSecond());
				w.println(e.weight);
				w.println();
			}
		} finally {
			w.close();
		}
	}

	private Collection<ExperimentResult> evaluateResults(TimeSeriesVariant dataVariant,
			Map<ModelVariant, WeightedDirectedGraph<Long, Double>> influenceGraphs,
			Map<ModelVariant, Integer> numEarlyRejections,
			Map<ModelVariant, WeightedDirectedGraph<Long, Explanation>> exp) throws IOException {
		Collection<ExperimentResult> results = new ArrayList<ExperimentResult>();
		for (Map.Entry<ModelVariant, WeightedDirectedGraph<Long, Double>> e : influenceGraphs.entrySet()) {
			ModelVariant modelVariant = e.getKey();
			if (!modelVariant.isFakeNetwork()) {
				String experimentId = dataVariant + "-" + modelVariant.getFilePrefix();
				String graphFileName = experimentId + ".graphml";
				e.getValue().writeGraphMl(new File(outPath, graphFileName));
				if (exp != null) {
					WeightedDirectedGraph<Long, Explanation> expGraph = exp.get(modelVariant);
					if (expGraph != null)
						writeExplanations(new File(outPath, experimentId + "-exp.txt"), expGraph);
				}

				WeightedDirectedGraph<Long, Double> fakeGraph = influenceGraphs.get(modelVariant.fake());

				DirectedSparseGraph<Long, Integer> network = networks.get(modelVariant.getNetwork());
				DirectedSparseGraph<Long, Integer> fakeNetwork = fakeNetworks.get(modelVariant.getNetwork());

				double earlyRejectRate = -1.0;
				Integer n = numEarlyRejections.get(modelVariant);
				if (n != null)
					earlyRejectRate = (double) n / network.getEdgeCount();

				// coverage of underlying network; assumes that the generated influence graph is a subgraph
				double networkCoverage = (double) e.getValue().getEdgeCount() / network.getEdgeCount();

				// coverage of network with permuted edges
				double fakeNetworkCoverage = (double) fakeGraph.getEdgeCount() / fakeNetwork.getEdgeCount();

				results.add(new ExperimentResult(dataVariant, modelVariant, graphFileName, network.getEdgeCount(),
						earlyRejectRate, networkCoverage, fakeNetworkCoverage));
			}
		}
		return results;
	}

	public void runExperiments(int[] timePeriods, File evalFile) throws Exception {
		long totalTimeGranger = 0, totalTimeSimilarity = 0;
		int numVariantsGranger = 0, numVariantsSimilarity = 0;

		PrintWriter w = new PrintWriter(new FileWriter(evalFile));
		try {
			w.println(ExperimentResult.getCsvHeader());

			// enumerate all data and model variants
			for (int timePeriodLength : timePeriods) {
				TimeSeriesVariant dataVariant = new TimeSeriesVariant(timePeriodLength);
				Map<Long, TimeSeriesData> seriesData = series.get(dataVariant);

				for (ModelVariant.ModelType modelType : ModelVariant.ModelType.values()) {
					if (modelType == ModelVariant.ModelType.GRANGER)	// skip untransformed GC experiments
						continue;

					logger.info("starting experiment " + dataVariant + "-" + modelType);
					Map<ModelVariant, WeightedDirectedGraph<Long, Explanation>> exp = null;
					if (generateExplanations)
						exp = new HashMap<ModelVariant, WeightedDirectedGraph<Long, Explanation>>();
					TemporalInfluenceModel<?> model = TemporalInfluenceModel.createModel(modelType, seriesData,
							topicPriorParam);
					ProgressReporter rep = new ProgressReporter(model);
					rep.start();
					long tStart = System.currentTimeMillis();
					int n = 0;
					try {
						Map<ModelVariant, WeightedDirectedGraph<Long, Double>> influenceGraphs =
								model.buildInfluenceGraphs(networks, fakeNetworks, exp);
						for (ExperimentResult r : evaluateResults(dataVariant, influenceGraphs,
								model.numEarlyRejections, exp)) {
							if (!r.modelVariant.isSubModel())
								n++;
							w.println(r);
						}
					} finally {
						rep.interrupt();
						w.flush();
					}

					PrintWriter ws = new PrintWriter(new FileWriter(new File(outPath,
							dataVariant + "-" + modelType + "-stats.csv")));
					try {
						model.writeStatistics(ws);
					} finally {
						ws.close();
					}

					long t = System.currentTimeMillis() - tStart;
					if (modelType == ModelVariant.ModelType.SIMILARITY) {
						totalTimeSimilarity += t;
						numVariantsSimilarity += n;
					} else {
						totalTimeGranger += t;
						numVariantsGranger += n;
					}
				}
			}
		} finally {
			w.close();
		}

		if (numVariantsSimilarity > 0) {
			logger.info("similarity-based experiments: " + ((totalTimeSimilarity / numVariantsSimilarity) / 1000) +
					"s per experiment");
		}
		if (numVariantsGranger > 0) {
			logger.info("GC-based experiments: " + ((totalTimeGranger / numVariantsGranger) / 1000) +
					"s per experiment");
		}
	}

	public static void main(String[] args) throws Exception {
		String commNetworkFileName = null, seriesFileName = null, followNetworkFileName = null;
		if (args.length > 1) {
			commNetworkFileName = args[0];
			seriesFileName = args[1];
			if (args.length > 2)
				followNetworkFileName = args[2];
		}

		InfluenceMeasurement im = new InfluenceMeasurement();
		int[] timePeriods = im.loadNetworksAndActions(followNetworkFileName, commNetworkFileName, seriesFileName);
		im.analyzeData(new File(outPath, "missing.csv"));
		im.runExperiments(timePeriods, new File(outPath, "eval-pre.csv"));
		logger.info("Finished.");
	}

}
