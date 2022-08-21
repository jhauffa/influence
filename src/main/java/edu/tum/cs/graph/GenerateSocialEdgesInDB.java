package edu.tum.cs.graph;

import java.io.File;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import edu.tum.cs.db.SocialGraphDao;
import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.SocialMediaDaoFactory;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.db.entities.SocialMediaUser;
import edu.tum.cs.db.loader.MessageLoader;
import edu.tum.cs.db.loader.MessageLoaderFactory;
import edu.tum.cs.math.dist.DiscreteDistribution;
import edu.tum.cs.nlp.corpus.Index;
import edu.tum.cs.nlp.corpus.MessageCorpus;
import edu.tum.cs.nlp.corpus.ProcessedMessage;
import edu.tum.cs.nlp.corpus.TokenizedMessage;
import edu.tum.cs.nlp.topic.FitART;
import edu.tum.cs.nlp.topic.model.ART;
import edu.tum.cs.nlp.topic.model.OnlineTopicModel;
import edu.tum.cs.time.IntervalCalculator;
import edu.tum.cs.time.TimeInterval;
import edu.tum.cs.util.ExperimentConfiguration;
import edu.tum.cs.util.arrays.Object2DArray;
import edu.tum.cs.util.io.Serializer;
import edu.uci.ics.jung.algorithms.scoring.PageRank;
import edu.uci.ics.jung.graph.DirectedGraph;

public class GenerateSocialEdgesInDB {

	private static final Logger logger = Logger.getLogger(GenerateSocialEdgesInDB.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(GenerateSocialEdgesInDB.class);

	private static final String modelPath = cfg.getProperty(ExperimentConfiguration.PROP_TOPIC_MODEL_PATH);
	private static final boolean online = cfg.getBooleanProperty(ExperimentConfiguration.PROP_USE_ONLINE_TOPIC_MODEL);
	private static final boolean commGraphMetrics = cfg.getLocalBooleanProperty("alwaysUseCommGraphMetrics", false);

	public static class ClassifyMessages implements MessageLoader.MessageTypeCallback {
		private final Set<Long> userIdSet;
		private final Map<Long, SocialNode> nodesMap;
		private final Map<Long, Map<Long, SocialEdge>> outgoingEdgesMap;
		private final Map<Long, Map<Long, SocialEdge>> incomingEdgesMap;
		private TimeInterval timeSlice;

		public ClassifyMessages(Set<Long> userIdSet, Map<Long, SocialNode> nodesMap,
				Map<Long, Map<Long, SocialEdge>> outgoingEdgesMap, Map<Long, Map<Long, SocialEdge>> incomingEdgesMap) {
			this.userIdSet = userIdSet;
			this.nodesMap = nodesMap;
			this.outgoingEdgesMap = outgoingEdgesMap;
			this.incomingEdgesMap = incomingEdgesMap;
		}

		public void setDate(TimeInterval timeSlice) {
			this.timeSlice = timeSlice;
		}

		@Override
		public boolean foundNonAddressiveMessage(SocialMediaMessage message, long senderId,
				Collection<Long> exposedIds) {
			nodesMap.get(senderId).incNumNonAddr(timeSlice);

			if (!exposedIds.isEmpty()) {
				for (long userId : exposedIds) {
					if (userIdSet.contains(userId))
						nodesMap.get(userId).addExposer(timeSlice, senderId);
				}
			}
			return true;
		}

		@Override
		public boolean foundAddressiveMessage(SocialMediaMessage message, long senderId,
				Collection<Long> recipientIds) {
			for (long recipientId : recipientIds) {
				if (userIdSet.contains(recipientId)) {
					SocialEdge edge = getSocialEdge(senderId, recipientId, outgoingEdgesMap, incomingEdgesMap);
					edge.incNumAddr(timeSlice);
				}
			}
			return true;
		}

		@Override
		public boolean foundSharedMessage(SocialMediaMessage message, long senderId, long originalSenderId) {
			if (userIdSet.contains(originalSenderId)) {
				SocialEdge edge = getSocialEdge(senderId, originalSenderId, outgoingEdgesMap, incomingEdgesMap);
				edge.incNumShared(timeSlice);
			}
			return false;
		}
	}

	private static int generateNodeAndExplicitEdges(SocialMediaDao dao, long userId, Set<Long> userIds,
			Map<Long, SocialNode> userId2SocialNode, Map<Long, Map<Long, SocialEdge>> outgoingEdgesMap,
			Map<Long, Map<Long, SocialEdge>> incomingEdgesMap) throws Exception {
		SocialMediaUser user = dao.getUser(userId);

		SocialNode node = new SocialNode(userId, true);
		node.reportedNumIncoming = user.getNumIncomingEdges();
		node.reportedNumOutgoing = user.getNumOutgoingEdges();
		userId2SocialNode.put(userId, node);

		/* Create explicit social edges for current node. Note that the graph represented by outgoingEdgesMap and
		 * incomingEdgesMap, and later stored in the SOCIAL_EDGE table, essentially contains two different types of
		 * edges: Explicit, directed edges represent a declared relationship and are in opposition to the actual flow
		 * of information (e.g. in case of Twitter, a follower subscribes to messages from the followee).
		 * Communication edges (added later by ClassifyMessages), point from the sender towards the recipient, and thus
		 * always agree with the flow of communication. An edge can be both an explicit and a communication edge. Proper
		 * separation of the two edge types is important. The explicit graph can be obtained by discarding all edges
		 * where isExplicit == false, while the communication graph is obtained by discarding edges where getTheta
		 * returns null for the specified time period.
		 */
		int numExplicitEdges = 0;
		user.getOutgoingEdges().retainAll(userIds);
		for (Long otherUserId : user.getOutgoingEdges()) {
			SocialEdge edge = getSocialEdge(userId, otherUserId, outgoingEdgesMap, incomingEdgesMap);
			edge.isExplicit = true;
			numExplicitEdges++;
		}
		user.getIncomingEdges().retainAll(userIds);
		for (Long otherUserId : user.getIncomingEdges()) {
			SocialEdge edge = getSocialEdge(otherUserId, userId, outgoingEdgesMap, incomingEdgesMap);
			edge.isExplicit = true;
			numExplicitEdges++;
		}
		return numExplicitEdges;
	}

	private static SocialEdge getSocialEdge(long fromUserId, Long toUserId,
			Map<Long, Map<Long, SocialEdge>> outgoingEdgesMap, Map<Long, Map<Long, SocialEdge>> incomingEdgesMap) {
		SocialEdge edge;
		Map<Long, SocialEdge> outgoingEdgesForUser = outgoingEdgesMap.get(fromUserId);
		edge = outgoingEdgesForUser.get(toUserId);
		if (edge == null) {
			edge = new SocialEdge(fromUserId, toUserId);
			outgoingEdgesForUser.put(toUserId, edge);
			incomingEdgesMap.get(toUserId).put(fromUserId, edge);
		}
		return edge;
	}

	private static Collection<TimeInterval> getArtTimeSlices(Date[] dates, int[] timePeriods) {
		Set<TimeInterval> timeSlices = new HashSet<TimeInterval>();
		for (Date centerDate : dates) {
			for (int length : timePeriods) {
				timeSlices.add(new TimeInterval(centerDate, -length));
				timeSlices.add(new TimeInterval(centerDate, length));
			}
		}
		return timeSlices;
	}

	private static void computeGraphMetrics(Set<Long> userIds, Map<Long, SocialNode> socialNodes,
			Map<Long, Map<Long, SocialEdge>> outgoingEdgesMap, boolean isExplicit) {
		long endTime, startTime;
		DirectedGraph<Long, Integer> graph = GraphCreator.createGraph(userIds, outgoingEdgesMap, !isExplicit);

		// compute PageRank
		startTime = System.currentTimeMillis();
		PageRank<Long, Integer> pagy = new PageRank<Long, Integer>(graph, 0.2);
		pagy.acceptDisconnectedGraph(true);
		pagy.evaluate();
		endTime = System.currentTimeMillis();
		logger.info("Calculated PageRank for " + graph.getEdgeCount() + " edges among " +
				graph.getVertexCount() + " nodes in " + (endTime - startTime) / 1000.0 + " secs");
		for (Map.Entry<Long, SocialNode> entry : socialNodes.entrySet()) {
			Long userId = entry.getKey();
			if (!userIds.contains(userId))
				continue;
			SocialNode socialNode = entry.getValue();
			socialNode.pageRank = pagy.getVertexScore(userId);
		}

		// compute Betweenness Centrality
		startTime = System.currentTimeMillis();
		BetweennessCentrality<Long, Integer> betweennessCentrality = new BetweennessCentrality<Long, Integer>(graph);
		betweennessCentrality.compute();
		endTime = System.currentTimeMillis();
		logger.info("Calculated betweeness centrality for " + graph.getEdgeCount() + " edges among " +
				graph.getVertexCount() + " nodes in " + (endTime - startTime) / 1000.0 + " secs");
		for (Map.Entry<Long, SocialNode> entry : socialNodes.entrySet()) {
			Long userId = entry.getKey();
			if (!userIds.contains(userId))
				continue;
			SocialNode socialNode = entry.getValue();
			socialNode.betweennessCentrality = betweennessCentrality.getCentrality(userId);
		}
	}

	private static double[] priorTheta;

	private static boolean isPrior(double[] theta) {
		return (DiscreteDistribution.distJS2(theta, priorTheta) < 1E-14);
	}

	public static void main(String[] args) throws Exception {
		SocialGraphDao graphDao = new SocialGraphDao();
		try {
			graphDao.createTables();
		} catch (SQLException ex) {
			logger.severe("tables already exist, aborting to avoid data corruption");
			return;
		}

		Date globalEndDate = cfg.getDateProperty(ExperimentConfiguration.PROP_END_DATE);
		IntervalCalculator ic = new IntervalCalculator(globalEndDate);
		Date globalStartDate = ic.getStartDate();
		logger.info("globalStartDate: " + globalStartDate + " - globalEndDate: " + globalEndDate);

		// Generate time slices
		Date[] dates = ic.getCenterDates();
		Collection<TimeInterval> artTimeSlices = getArtTimeSlices(dates,
				cfg.getIntListProperty(ExperimentConfiguration.PROP_TIME_PERIODS));
		for (TimeInterval iv : artTimeSlices)
			logger.info(iv.toString());

		// Generate nodes and explicit incoming/outgoing edges
		SocialMediaDao mediaDao = SocialMediaDaoFactory.createDao();
		Set<Long> userIdSet = mediaDao.getUserIds(true);
		Map<Long, SocialNode> userId2SocialNode = new HashMap<Long, SocialNode>(userIdSet.size());
		Map<Long, Map<Long, SocialEdge>> userId2OutgoingEdgesMap =
				new HashMap<Long, Map<Long, SocialEdge>>(userIdSet.size());
		Map<Long, Map<Long, SocialEdge>> userId2IncomingEdgesMap =
				new HashMap<Long, Map<Long, SocialEdge>>(userIdSet.size());
		for (long userId : userIdSet) {
			userId2OutgoingEdgesMap.put(userId, new HashMap<Long, SocialEdge>());
			userId2IncomingEdgesMap.put(userId, new HashMap<Long, SocialEdge>());
		}

		logger.info("Generating nodes and explicit edges");
		int numUsers = 0;
		int numExplicitEdges = 0;
		for (long userId : userIdSet) {
			// This counts both incoming and outgoing edges, and will thus overestimate the total number of edges in the
			// graph, but the result is only used to determine whether the graph contains any explicit edges at all.
			numExplicitEdges += generateNodeAndExplicitEdges(mediaDao, userId, userIdSet, userId2SocialNode,
					userId2OutgoingEdgesMap, userId2IncomingEdgesMap);

			if (++numUsers % 1000 == 0)
				logger.info("Processed " + numUsers + " users");
		}
		logger.info("Processed " + numUsers + " users, done");

		// Load existing ART model
		Index<String> bagOfWords = Serializer.loadObjectFromFile(new File(modelPath, FitART.bowFileName));
		bagOfWords.setReadOnly();
		Index<Long> bagOfPersons = Serializer.loadObjectFromFile(new File(modelPath, FitART.bopFileName));
		bagOfPersons.setReadOnly();

		ART art = null;
		Map<TimeInterval, ART> artModels = null;
		if (online) {
			artModels = new HashMap<TimeInterval, ART>();
			OnlineTopicModel<ProcessedMessage, ART> oart =
					Serializer.loadObjectFromFile(new File(modelPath, "onlineART-model-full.ser.gz"));
			// first model in model history is most current
			int dateIdx = dates.length - 1;
			Date modelEndDate = globalEndDate;
			for (ART curModel : oart.getModelHistory()) {
				curModel.configureForQuerying();
				for (TimeInterval timeSlice : artTimeSlices) {
					if (timeSlice.getStartDate().before(modelEndDate))
						artModels.put(timeSlice, curModel);
				}
				if (dateIdx >= 0)
					modelEndDate = dates[dateIdx--];
			}
		} else {
			art = Serializer.loadObjectFromFile(new File(modelPath, "art-model-final.ser.gz"));
			art.configureForQuerying();
			priorTheta = art.getAlpha().clone();
			DiscreteDistribution.normalize(priorTheta);
		}

		SocialNode globalThetaContainer = new SocialNode(graphDao.getGlobalNodeId(), false);
		userId2SocialNode.put(globalThetaContainer.id, globalThetaContainer);

		// Load users' messages
		MessageLoader loader = MessageLoaderFactory.createMessageLoader(mediaDao, userIdSet);
		ClassifyMessages classifier = new ClassifyMessages(userIdSet, userId2SocialNode, userId2OutgoingEdgesMap,
				userId2IncomingEdgesMap);
		loader.setCallback(classifier);
		loader.cacheMessages(globalStartDate, globalEndDate);

		int numSlices = 0;
		for (TimeInterval timeSlice : artTimeSlices) {
			logger.info("Slice " + timeSlice + " starting");
			classifier.setDate(timeSlice);

			// Generate missing communication edges, and build corpora for ART query
			List<TokenizedMessage> addrMessages = new ArrayList<TokenizedMessage>();
			List<TokenizedMessage> nonAddrMessages = new ArrayList<TokenizedMessage>();
			for (long userId : userIdSet) {
				List<TokenizedMessage> docs = loader.loadMessages(userId, timeSlice.getStartDate(),
						timeSlice.getEndDate());
				for (TokenizedMessage doc : docs) {
					if ((doc.getRecipients().size() == 1) &&
						(doc.getRecipients().iterator().next().equals(doc.getSender())))
						nonAddrMessages.add(doc);
					else
						addrMessages.add(doc);
				}
			}
			MessageCorpus<Long, TokenizedMessage> addrCorpus = new MessageCorpus<Long, TokenizedMessage>(bagOfWords,
					bagOfPersons, addrMessages);
			MessageCorpus<Long, TokenizedMessage> nonAddrCorpus = new MessageCorpus<Long, TokenizedMessage>(bagOfWords,
					bagOfPersons, nonAddrMessages);

			if (online) {
				art = artModels.get(timeSlice);
				priorTheta = art.getAlpha().clone();
				DiscreteDistribution.normalize(priorTheta);
			}

			// Query corpus of addressive messages
			ART.QueryResult addrState = art.query(addrCorpus);
			double[] addrGlobalTheta = addrState.estimateAggregatedTheta(true, null, true, null)[0];
			double[][] addrAuthorThetas = addrState.estimateAggregatedTheta(false, null, true, null);
			double[][] addrRecipientThetas = addrState.estimateAggregatedTheta(true, null, false, null);
			Object2DArray<double[]> addrAuthorRecipientThetas = addrState.estimateTheta();
			addrState = null;

			// Query corpus of non-addressive messages
			ART.QueryResult nonAddrState = art.query(nonAddrCorpus);
			double[] nonAddrGlobalTheta = nonAddrState.estimateAggregatedTheta(true, null, true, null)[0];
			double[][] nonAddrAuthorThetas = nonAddrState.estimateAggregatedTheta(false, null, true, null);
			double[][] nonAddrRecipientThetas = new double[bagOfPersons.getNumUniqueElements()][];
			for (Long userId : userIdSet) {
				Set<Integer> outgoingIds = new HashSet<Integer>();
				SocialNode node = userId2SocialNode.get(userId);
				Set<Long> exposingIds = node.getExposers(timeSlice);
				if (exposingIds != null) {
					for (Long exposingUserId : exposingIds) {
						int id = bagOfPersons.getElementId(exposingUserId);
						if (id >= 0)
							outgoingIds.add(id);
					}
				} else {
					for (Map.Entry<Long, SocialEdge> e : userId2OutgoingEdgesMap.get(userId).entrySet()) {
						if (e.getValue().isExplicit) {
							int id = bagOfPersons.getElementId(e.getKey());
							if (id >= 0)
								outgoingIds.add(id);
						}
					}
				}
				int id = bagOfPersons.getElementId(userId);
				if (id >= 0)
					nonAddrRecipientThetas[id] = nonAddrState.estimateAggregatedTheta(true, outgoingIds, true, null)[0];
			}
			nonAddrState = null;

			// Set thetas of nodes and edges
			globalThetaContainer.setTheta(timeSlice, nonAddrGlobalTheta);
			globalThetaContainer.setSenderTheta(timeSlice, addrGlobalTheta);

			for (Long userId : userIdSet) {
				SocialNode node = userId2SocialNode.get(userId);
				int senderId = bagOfPersons.getElementId(userId);
				if (senderId < 0)
					continue;

				// Could enforce an entropy threshold here, but we're already suffering from data sparsity...
				if (!isPrior(nonAddrAuthorThetas[senderId]))
					node.setTheta(timeSlice, nonAddrAuthorThetas[senderId]);
				if (!isPrior(nonAddrRecipientThetas[senderId]))
					node.setExposureTheta(timeSlice, nonAddrRecipientThetas[senderId]);
				if (!isPrior(addrAuthorThetas[senderId]))
					node.setSenderTheta(timeSlice, addrAuthorThetas[senderId]);
				if (!isPrior(addrRecipientThetas[senderId]))
					node.setRecipientTheta(timeSlice, addrRecipientThetas[senderId]);

				Map<Long, SocialEdge> outgoingEdges = userId2OutgoingEdgesMap.get(userId);
				for (SocialEdge edge : outgoingEdges.values()) {
					int recipientId = bagOfPersons.getElementId(edge.toUserId);
					if (recipientId >= 0) {
						double[] edgeTheta = addrAuthorRecipientThetas.get(senderId, recipientId);
						if ((edgeTheta != null) && !isPrior(edgeTheta))
							edge.setTheta(timeSlice, edgeTheta);
					}
				}
			}

			logger.info("Done with " + (++numSlices) + " of " + artTimeSlices.size());
		}
		logger.info("Done with edge and node theta calculations");

		// Save edges to DB
		logger.info("Saving edges to DB");
		int numEdges = 0;
		numUsers = 0;
		for (long userId : userIdSet) {
			Collection<SocialEdge> socialEdges = userId2OutgoingEdgesMap.get(userId).values();
			graphDao.saveSocialEdges(socialEdges);
			numEdges += socialEdges.size();

			if (++numUsers % 1000 == 0)
				logger.info("Saved " + numEdges + " edges of " + numUsers + " users");
		}
		logger.info("Saved " + numEdges + " edges of " + numUsers + " users, done");

		// Calculate "got shared" values
		numUsers = 0;
		for (long userId : userIdSet) {
			Map<Long, SocialEdge> incomingEdges = userId2IncomingEdgesMap.get(userId);
			SocialNode node = userId2SocialNode.get(userId);
			for (TimeInterval timeSlice : artTimeSlices) {
				int gotShared = 0;
				for (SocialEdge incomingEdge : incomingEdges.values())
					gotShared += incomingEdge.getNumShared(timeSlice);
				node.setNumGotShared(timeSlice, gotShared);
			}

			if (++numUsers % 1000 == 0)
				logger.info("Calculated the 'got shared' values for " + numUsers + " users");
		}

		// Compute graph metrics
		logger.info("Computing graph metrics");
		computeGraphMetrics(userIdSet, userId2SocialNode, userId2OutgoingEdgesMap,
				(numExplicitEdges > 0) && !commGraphMetrics);

		// Save nodes to DB
		logger.info("Saving nodes to DB");
		graphDao.saveSocialNodes(userId2SocialNode.values());

		logger.info("Finished");
	}

}
