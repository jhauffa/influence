package edu.tum.cs.influence.meso;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import edu.tum.cs.graph.SocialEdge;
import edu.tum.cs.graph.SocialNode;
import edu.tum.cs.influence.meso.ModelVariant.ModelType;
import edu.tum.cs.math.dist.DirichletDistribution;
import edu.tum.cs.math.dist.DiscreteDistribution;
import edu.tum.cs.nlp.corpus.ProcessedMessage;
import edu.tum.cs.nlp.topic.model.ART;
import edu.tum.cs.nlp.topic.model.OnlineTopicModel;
import edu.tum.cs.time.IntervalCalculator;
import edu.tum.cs.time.TimeInterval;
import edu.tum.cs.util.ExperimentConfiguration;
import edu.tum.cs.util.io.Serializer;

public class InfluenceExperimentsAnova {

	private static final Logger logger = Logger.getLogger(InfluenceExperimentsAnova.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(InfluenceExperimentsAnova.class);
	private static final SimpleDateFormat df = new SimpleDateFormat("yyyy.MM.dd");

	public static final String PROP_OUT_PATH = "outputPath";
	public static final String FILE_EVAL_PREDICTION = "eval-prediction.csv";
	public static final String FILE_EVAL_PARAMETERS = "eval-parameters.csv";

	private static final boolean splitAnova = cfg.getLocalBooleanProperty("splitAnova");
	private static final boolean online = cfg.getBooleanProperty(ExperimentConfiguration.PROP_USE_ONLINE_TOPIC_MODEL);
	private static final String modelPath = cfg.getProperty(ExperimentConfiguration.PROP_TOPIC_MODEL_PATH);
	private static final String outPath = cfg.getLocalProperty(PROP_OUT_PATH, ".");

	private final AtomicInteger numFinishedExperiments = new AtomicInteger();
	private final int numTopics;
	private final Date[] dates;
	private final DirichletDistribution[] artTopicPrior;
	private final double[][] priorTheta;
	private final SocialGraphCache socialGraph;
	private final int longestTimePeriod;

	private class ExperimentTask implements Callable<ExperimentResult[]> {
		private static final double alphaL1 = 0.0001;
		private static final double alphaL2 = 0.0001;	// Elastic Net ratio = 0.5

		private final ModelVariant variant;
		private final int artIdx;
		private final int numCoeff;
		private final PredictiveModel constModel;
		private final List<InfluencedUser> users;

		public ExperimentTask(ModelVariant variant, int artIdx, List<InfluencedUser> users) {
			this.variant = variant;
			this.artIdx = artIdx;
			this.users = users;

			// assume that all InfluencedUsers have the same parameters at this point
			numCoeff = users.get(0).getNumVariables();
			constModel = new PredictiveModel(numTopics, numCoeff, alphaL1, alphaL2);
		}

		private PredictionResult evaluatePrediction(ModelType type, InfluencedUser user, PredictiveModel model) {
			double[] predictedTheta;
			switch (type) {
			case BASELINE_UNIFORM:
				predictedTheta = priorTheta[artIdx];
				break;
			case BASELINE_RANDOM:
				predictedTheta = artTopicPrior[artIdx].sample();
				break;
			case BASELINE_PREVIOUS:
				predictedTheta = user.getBaselineTheta();
				break;
			default:
				predictedTheta = model.predictTheta(user);
				break;
			}

			double[] actualTheta = user.getThetaToPredict();
			double distJs = DiscreteDistribution.distJS2(predictedTheta, actualTheta);
			double distCos = 1.0 - DiscreteDistribution.simCosine(predictedTheta, actualTheta);
			double distZeroOne = DiscreteDistribution.distZeroOne(predictedTheta, actualTheta);
			return new PredictionResult(user.getFromUserId(), user.getToUserId(), distJs, distCos, distZeroOne);
		}

		@Override
		public ExperimentResult[] call() throws Exception {
			try {
				return compute();
			} catch (Throwable ex) {
				throw new Exception("computation of variant '" + variant + "' failed", ex);
			}
		}

		public ExperimentResult[] compute() {
			ExperimentResult[] results;
			ModelType type = variant.getType();
			boolean isScim = ((type == ModelType.SCIM) ||
					(type == ModelType.SCIM_REVERSE) || (type == ModelType.SCIM_SHUFFLE));
			if (isScim) {
				// for reasons of efficiency, compute model variants with constant and variable coefficients in one go
				results = new ExperimentResult[2];
				results[0] = new ExperimentResult(variant);
				ModelVariant variantConst = new ModelVariant(variant);
				variantConst.setType(ModelType.CONST_COEFF);
				results[1] = new ExperimentResult(variantConst);
			} else {
				results = new ExperimentResult[1];
				results[0] = new ExperimentResult(variant);
			}

			// update neighborhood thetas for all entities of the current experiment
			if (isScim) {
				for (InfluencedUser user : users)
					user.computeNeighbourhoodThetas(variant);
			} else if (type == ModelType.EMPTY_NEIGHBORHOOD) {
				for (InfluencedUser user : users)
					user.setEmptyNeighborhood();
			}

			// split into training and test set
			int numTrainUsers = users.size() / 2;
			List<InfluencedUser> trainUsers = users.subList(0, numTrainUsers);
			List<InfluencedUser> testUsers = users.subList(numTrainUsers, users.size());

			PredictiveModel model = null;
			if (isScim || (type == ModelType.EMPTY_NEIGHBORHOOD)) {
				// perform parameter fitting on training set
				model = new PredictiveModel(numTopics, numCoeff, alphaL1, alphaL2);
				SummaryStatistics trainErrorStats = model.fitParameters(trainUsers);
				double trainErrorMean = trainErrorStats.getMean();
				double trainErrorStdDev = trainErrorStats.getStandardDeviation();
				double trainErrorUserMean = 0.0;
				double trainErrorUserStdDev = 0.0;
				if (isScim) {
					SummaryStatistics trainErrorUserStats = model.computePerUserError(trainUsers);
					trainErrorUserMean = trainErrorUserStats.getMean();
					trainErrorUserStdDev = trainErrorUserStats.getStandardDeviation();
				}
				results[0].setParameterFittingResult(new ParameterFittingResult(model, trainErrorMean, trainErrorStdDev,
						trainErrorUserMean, trainErrorUserStdDev));
			}

			for (InfluencedUser user : testUsers) {
				results[0].addPredictionResult(evaluatePrediction(type, user, model));
				if (isScim)
					results[1].addPredictionResult(evaluatePrediction(type, user, constModel));
			}

			// experiment done - report progress
			int n = numFinishedExperiments.incrementAndGet();
			if ((n % 50) == 0)
				logger.info(n + " experiments done.");

			return results;
		}
	}

	private class ModelVariantExperimentTask implements Callable<ModelVariantExperiment> {
		private final Date timePeriodBase;
		private final int timePeriodLength;
		private final int artIdx;
		private final EntityType entity;
		private final ExecutorService pool;
		private final List<InfluencedUser[]> users;

		public ModelVariantExperimentTask(Date timePeriodBase, int timePeriodLength, int artIdx, EntityType entity,
				List<SocialEdge> testEdges, ExecutorService pool) {
			this.timePeriodBase = timePeriodBase;
			this.timePeriodLength = timePeriodLength;
			this.artIdx = artIdx;
			this.entity = entity;
			this.pool = pool;

			users = prepareUsers(testEdges);
		}

		private List<InfluencedUser[]> prepareUsers(List<SocialEdge> testEdges) {
			List<InfluencedUser[]> users = new ArrayList<InfluencedUser[]>();

			for (SocialEdge edge : testEdges) {
				SocialNode fromUserNode = socialGraph.getNode(edge.fromUserId);
				SocialNode toUserNode = socialGraph.getNode(edge.toUserId);

				if (entity == EntityType.NODE) {
					InfluencedUserNode userFw = new InfluencedUserNode(fromUserNode, edge.toUserId, socialGraph,
							timePeriodBase, timePeriodLength, longestTimePeriod, false, priorTheta[artIdx],
							priorTheta[artIdx + 1]);
					userFw.setThetas(fromUserNode);
					InfluencedUserNode userRev = new InfluencedUserNode(fromUserNode, edge.toUserId, socialGraph,
							timePeriodBase, timePeriodLength, longestTimePeriod, true, priorTheta[artIdx],
							priorTheta[artIdx + 1]);
					userRev.setThetas(fromUserNode);
					users.add(new InfluencedUser[] { userFw, userRev });
				} else {
					InfluencedUserEdge userFw = new InfluencedUserEdge(fromUserNode, edge.toUserId, socialGraph,
							timePeriodBase, timePeriodLength, longestTimePeriod, false, priorTheta[artIdx],
							priorTheta[artIdx + 1]);
					userFw.setThetas(fromUserNode, toUserNode, edge,
							socialGraph.getOutgoingEdges(toUserNode.id).get(fromUserNode.id));
					InfluencedUserEdge userRev = new InfluencedUserEdge(fromUserNode, edge.toUserId, socialGraph,
							timePeriodBase, timePeriodLength, longestTimePeriod, true, priorTheta[artIdx],
							priorTheta[artIdx + 1]);
					userRev.setThetas(fromUserNode, toUserNode, edge,
							socialGraph.getOutgoingEdges(toUserNode.id).get(fromUserNode.id));
					users.add(new InfluencedUser[] { userFw, userRev });
				}
			}

			return users;
		}

		private List<InfluencedUser> cloneUsers(boolean reverse, boolean shuffle) {
			List<InfluencedUser> clonedUsers = new ArrayList<InfluencedUser>(users.size());
			for (InfluencedUser[] user : users) {
				if (shuffle)
					reverse = (Math.random() < 0.5);
				if (reverse)
					clonedUsers.add(user[1].clone());
				else
					clonedUsers.add(user[0].clone());
			}
			return clonedUsers;
		}

		@Override
		public ModelVariantExperiment call() throws Exception {
			ModelVariantExperiment experiment = new ModelVariantExperiment(timePeriodBase);

			ExecutorCompletionService<ExperimentResult[]> tasks =
					new ExecutorCompletionService<ExperimentResult[]>(pool);
			int numTasks = 0;
			for (ModelType type : ModelType.values()) {
				if ((type == ModelType.SCIM) || (type == ModelType.SCIM_REVERSE) || (type == ModelType.SCIM_SHUFFLE)) {
					boolean reverse = (type == ModelType.SCIM_REVERSE);
					boolean shuffle = (type == ModelType.SCIM_SHUFFLE);
					for (IndicatorFunction indicator : IndicatorFunction.getApplicable(socialGraph.isExplicit(),
							socialGraph.isSymmetric())) {
						for (WeightFunction weight : WeightFunction.getApplicable(socialGraph.isExplicit(),
								socialGraph.isSymmetric())) {
							tasks.submit(new ExperimentTask(new ModelVariant(timePeriodLength, entity, type, indicator,
									weight), artIdx, cloneUsers(reverse, shuffle)));
							numTasks++;
						}
					}
				// baselines (including EMPTY_NEIGHBORHOOD); CONST_COEFF is computed by the SCIM tasks
				} else if (type != ModelType.CONST_COEFF) {
					tasks.submit(new ExperimentTask(new ModelVariant(timePeriodLength, entity, type, null, null),
							artIdx, cloneUsers(false, false)));
					numTasks++;
				}
			}

			while (numTasks-- > 0) {
				try {
					ExperimentResult[] results = tasks.take().get();
					for (ExperimentResult r : results)
						experiment.addResult(r);
				} catch (ExecutionException ex) {
					logger.log(Level.SEVERE, "error in experiment", ex.getCause());
				}
			}
			return experiment;
		}
	}


	public InfluenceExperimentsAnova(SocialGraphCache socialGraph, Date[] dates, int[] timePeriods,
			double[][] artModelAlpha) {
		this.socialGraph = socialGraph;
		this.dates = dates;

		int maxPeriod = 0;
		for (int p : timePeriods) {
			if (p > maxPeriod)
				maxPeriod = p;
		}
		longestTimePeriod = maxPeriod;

		artTopicPrior = new DirichletDistribution[artModelAlpha.length];
		priorTheta = new double[artModelAlpha.length][];
		for (int i = 0; i < artModelAlpha.length; i++) {
			artTopicPrior[i] = new DirichletDistribution(artModelAlpha[i]);
			priorTheta[i] = artModelAlpha[i].clone();
			DiscreteDistribution.normalize(priorTheta[i]);
		}
		numTopics = priorTheta[0].length;
	}

	public void startExperiments(List<SocialEdge>[][][] testEdges, int[] timePeriods) throws Exception {
		ExperimentSeries resultContainer = new ExperimentSeries();

		numFinishedExperiments.set(0);
		int numThreads = cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_THREADS,
				Runtime.getRuntime().availableProcessors());
		ExecutorService pool = Executors.newFixedThreadPool(numThreads);
		ExecutorCompletionService<ModelVariantExperiment> tasks =
				new ExecutorCompletionService<ModelVariantExperiment>(pool);
		ExecutorService subPool = Executors.newFixedThreadPool(numThreads);
		try {
			for (int i = 0; i < dates.length; i++) {
				logger.info("Starting for date: " + df.format(dates[i]));

				// map center date to index of OART model for observation period; center dates are ordered (dates[0] is
				// last center date)
				int artIdx = online ? (dates.length - (i + 1)) : 0;

				// enumerate all possible model variant experiments
				int numTasks = 0;
				for (int j = 0; j < timePeriods.length; j++) {
					for (EntityType entity : EntityType.getApplicable(socialGraph.isExplicit())) {
						if (!testEdges[i][j][entity.value].isEmpty()) {
							tasks.submit(new ModelVariantExperimentTask(dates[i], timePeriods[j], artIdx, entity,
									testEdges[i][j][entity.value], subPool));
							numTasks++;
						} else {
							logger.warning("no data for experiment " + timePeriods[j] + ";" + entity + ", skipping");
						}
					}
				}

				// wait for experiments to finish and save results
				ModelVariantExperiment experimentsForDate = new ModelVariantExperiment(dates[i]);
				while (numTasks-- > 0)
					experimentsForDate.merge(tasks.take().get());
				resultContainer.addExperiment(experimentsForDate);

				logger.info("Done with date: " + df.format(dates[i]));
			}
		} finally {
			subPool.shutdown();
			pool.shutdown();
		}

		resultContainer.writeParameterFittingResultsCsv(new File(outPath, FILE_EVAL_PARAMETERS));
		resultContainer.writePredictionResultsCsv(new File(outPath, FILE_EVAL_PREDICTION));
		WeightFunction.writeAvgStatsCsv(new File(outPath, "neighborhood-avg.csv"));
		WeightFunction.writeEmptyNeighborhoodStatsCsv(new File(outPath, "neighborhood-empty.csv"));
	}

	private Map<Integer, List<SocialEdge>> classifyEntities(int[] timePeriods, Set<Long> excludedUserIds) {
		// collect nodes/edges that match the criteria
		Map<Integer, List<SocialEdge>> entitiesBase = new HashMap<Integer, List<SocialEdge>>();
		int[] numNodesPerDate = new int[dates.length];
		int[] numEdgesPerDate = new int[dates.length];
		for (SocialNode node : socialGraph.getNodes()) {
			Map<Long, SocialEdge> outgoingEdges = socialGraph.getOutgoingEdges(node.id);

			// select nodes
			int entityMask = 0;
			if (node.isCore && !excludedUserIds.contains(node.id)) {
				for (int i = 0; i < dates.length; i++) {
					if ((node.getTheta(new TimeInterval(dates[i], -longestTimePeriod)) != null) &&
						(node.getTheta(new TimeInterval(dates[i], longestTimePeriod)) != null)) {
						entityMask |= (1 << (EntityType.NODE.value * dates.length + i));
						numNodesPerDate[i]++;
					}
				}
				if (entityMask != 0) {
					List<SocialEdge> edges = entitiesBase.get(entityMask);
					if (edges == null) {
						edges = new ArrayList<SocialEdge>();
						entitiesBase.put(entityMask, edges);
					}
					edges.add(new SocialEdge(node.id, node.id));	// SocialEdge is used as container for node ID only!
				}
			}

			// select edges
			for (SocialEdge edge : outgoingEdges.values()) {
				if (!node.isCore || !socialGraph.getNode(edge.toUserId).isCore ||
					(excludedUserIds.contains(edge.fromUserId) && excludedUserIds.contains(edge.toUserId)))
					continue;

				entityMask = 0;
				for (int i = 0; i < dates.length; i++) {
					// test if addressive communication happened in both date - timePeriod and date + timePeriod
					if ((edge.getTheta(new TimeInterval(dates[i], -longestTimePeriod)) != null) &&
						(edge.getTheta(new TimeInterval(dates[i], longestTimePeriod)) != null)) {
						entityMask |= (1 << (EntityType.EDGE_COMMUNICATION.value * dates.length + i));
						numEdgesPerDate[i]++;

						if (edge.isExplicit)
							entityMask |= (1 << (EntityType.EDGE_COMMUNICATION_EXPLICIT.value * dates.length + i));
					}
				}

				if (entityMask != 0) {
					List<SocialEdge> edges = entitiesBase.get(entityMask);
					if (edges == null) {
						edges = new ArrayList<SocialEdge>();
						entitiesBase.put(entityMask, edges);
					}
					edges.add(edge);
				}
			}
		}
		for (int i = 0; i < dates.length; i++) {
			logger.info("date " + df.format(dates[i]) + ": " + numNodesPerDate[i] + " nodes, " +
					numEdgesPerDate[i] + " edges");
		}
		return entitiesBase;
	}

	@SuppressWarnings("unchecked")
	private List<SocialEdge>[][][] collectTestEdges(int[] timePeriods, Set<Long> excludedUserIds) {
		Map<Integer, List<SocialEdge>> entitiesBase = classifyEntities(timePeriods, excludedUserIds);
		List<SocialEdge>[][][] testSet = new List[dates.length][timePeriods.length][EntityType.values().length];
		for (int i = 0; i < testSet.length; i++) {
			for (int j = 0; j < testSet[0].length; j++)
				for (int k = 0; k < testSet[0][0].length; k++)
					testSet[i][j][k] = new ArrayList<SocialEdge>();
		}

		for (Map.Entry<Integer, List<SocialEdge>> e : entitiesBase.entrySet()) {
			for (int i = 0; i < EntityType.values().length; i++) {
				for (int j = 0; j < dates.length; j++) {
					if ((e.getKey() & (1 << (i * dates.length + j))) != 0) {
						for (int k = 0; k < timePeriods.length; k++)
							testSet[j][k][i].addAll(e.getValue());
					}
				}
			}
		}

		for (int i = 0; i < testSet.length; i++) {
			for (int j = 0; j < testSet[0].length; j++)
				for (int k = 0; k < testSet[0][0].length; k++)
					Collections.shuffle(testSet[i][j][k]);
		}
		return testSet;
	}

	@SuppressWarnings("unchecked")
	private List<SocialEdge>[][][] collectTestEdgesAnova(int[] timePeriods, Set<Long> excludedUserIds) {
		Map<Integer, List<SocialEdge>> entitiesBase = classifyEntities(timePeriods, excludedUserIds);

		// Iterate over all possible masks starting from the least general; distribute entities equally but randomly
		// over all entities and date slots permitted by the mask.
		List<SocialEdge>[][] entitiesSplit = new List[EntityType.values().length][dates.length];
		for (EntityType entity : EntityType.values())
			for (int i = 0; i < dates.length; i++)
				entitiesSplit[entity.value][i] = new ArrayList<SocialEdge>();
		for (int i = 1; i < (1 << (EntityType.values().length * dates.length)); i++) {
			List<SocialEdge> edges = entitiesBase.get(i);
			if (edges == null)
				continue;
			Collections.shuffle(edges);

			int partSize = edges.size() / Integer.bitCount(i);
			int partOffset = 0;
			for (int j = 0; j < EntityType.values().length; j++) {
				for (int k = 0; k < dates.length; k++) {
					if ((i & (1 << (j * dates.length + k))) == 0)
						continue;
					entitiesSplit[j][k].addAll(edges.subList(partOffset, partOffset + partSize));
					partOffset += partSize;
				}
			}
		}

		// partition entity set of each time period / type into equal-sized subsets for the different period lengths
		int sampleSize = Integer.MAX_VALUE;	// TODO: should be computed separately for edges and nodes
		for (EntityType entity : EntityType.getApplicable(socialGraph.isExplicit())) {
			for (int i = 0; i < dates.length; i++) {
				int n = entitiesSplit[entity.value][i].size();
				logger.info("date " + df.format(dates[i]) + ", " + entity + ": " + n);
				if (n < sampleSize)
					sampleSize = n;
			}
		}
		sampleSize /= timePeriods.length;
		if (sampleSize == 0)
			throw new RuntimeException("need at least one sample per combination of factor values");
		logger.info("sample size = " + sampleSize);

		List<SocialEdge>[][][] testSet = new List[dates.length][timePeriods.length][EntityType.values().length];
		for (int i = 0; i < testSet.length; i++) {
			for (int j = 0; j < testSet[0][0].length; j++) {
				Iterator<SocialEdge> it = entitiesSplit[j][i].iterator();
				for (int k = 0; k < testSet[0].length; k++) {
					testSet[i][k][j] = new ArrayList<SocialEdge>();
					if (!it.hasNext())
						continue;

					for (int l = 0; l < sampleSize; l++)
						testSet[i][k][j].add(it.next());

					// shuffle edges to get a random split into training and test data for parameter fitting
					Collections.shuffle(testSet[i][k][j]);
				}
			}
		}
		return testSet;
	}

	public static void main(String[] args) throws Exception {
		Set<Long> excludedUserIds = Collections.emptySet();
		String excludedIdsFile = cfg.getLocalProperty("exclusionList", "");
		if (!excludedIdsFile.isEmpty()) {
			int numExcludedIds = cfg.getLocalIntProperty("numExcludedUsers", Integer.MAX_VALUE);
			excludedUserIds = ExperimentConfiguration.loadUserIds(excludedIdsFile, numExcludedIds);
		}

		IntervalCalculator ic = new IntervalCalculator(cfg.getDateProperty(ExperimentConfiguration.PROP_END_DATE));
		Date[] dates = ic.getCenterDates();
		int[] timePeriods = cfg.getIntListProperty(ExperimentConfiguration.PROP_TIME_PERIODS);

		// load ART model(s) and retrieve final alpha value
		double[][] alpha;
		if (online) {
			OnlineTopicModel<ProcessedMessage, ART> oart =
					Serializer.loadObjectFromFile(new File(modelPath, "onlineART-model-full.ser.gz"));
			alpha = new double[oart.getModelHistory().size()][];
			// store alpha values chronologically: alpha[0] belongs to model of first observation period
			int idx = alpha.length - 1;
			for (ART art : oart.getModelHistory())
				alpha[idx--] = art.getAlpha();
		} else {
			ART art = Serializer.loadObjectFromFile(new File(modelPath, "art-model-final.ser.gz"));
			alpha = new double[][] { art.getAlpha(), art.getAlpha() };
		}

		logger.info("Loading the social graph, its cliques and communities.");
		SocialGraphCache socialGraph = new SocialGraphCache(new File(outPath, "vertexstat.ser"),
				new File(outPath, "communities-perc.ser"), new File(outPath, "communities-edge.ser"));
		logger.info("An explicit social graph is" + (socialGraph.isExplicit() ? " " : " not ") + "present, graph is" +
				(socialGraph.isSymmetric() ? " " : " not ") + "symmetric.");

		InfluenceExperimentsAnova tiv = new InfluenceExperimentsAnova(socialGraph, dates, timePeriods, alpha);
		logger.info("Collecting test edges.");
		List<SocialEdge>[][][] testEdges;
		if (splitAnova)
			testEdges = tiv.collectTestEdgesAnova(timePeriods, excludedUserIds);
		else
			testEdges = tiv.collectTestEdges(timePeriods, excludedUserIds);
		Serializer.saveObjectToFile(testEdges, new File(outPath, "testEdges.ser"));

		logger.info("Starting experiments.");
		tiv.startExperiments(testEdges, timePeriods);
		logger.info("Finished!");
	}

}
