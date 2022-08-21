package edu.tum.cs.influence.micro;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.SocialMediaDaoFactory;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.db.loader.MessageLoader;
import edu.tum.cs.db.loader.MessageLoaderFactory;
import edu.tum.cs.graph.GraphCreator;
import edu.tum.cs.influence.micro.EvaluateNetworks.EvaluationResult;
import edu.tum.cs.math.dist.DiscreteDistributionSampler;
import edu.tum.cs.time.IntervalCalculator;
import edu.tum.cs.util.ExperimentConfiguration;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;

public class ExtractTweetSample {

	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(InfluenceMeasurement.class);
	private static final int defaultSampleSize = 100;

	private static class EdgeData {
		public final boolean[] isMember;
		public final double[] magnitude;
		public final Set<Integer> timePeriods = new HashSet<Integer>();

		public EdgeData(int numCategories) {
			isMember = new boolean[numCategories];
			magnitude = new double[numCategories];
		}
	}

	private static void addEdgeToSample(Pair<Long> v, Map<Pair<Long>, EdgeData> sampleEdges, Set<Long> vertices, int c,
			int numCategories, double magnitude, int timePeriod) {
		EdgeData data = sampleEdges.get(v);
		if (data == null) {
			data = new EdgeData(numCategories);
			sampleEdges.put(v, data);
		}
		data.isMember[c] = true;
		data.magnitude[c] = magnitude;
		data.timePeriods.add(timePeriod);

		if (vertices != null) {
			vertices.add(v.getFirst());
			vertices.add(v.getSecond());
		}
	}

	private static <E> void sampleGraph(DirectedSparseGraph<Long, E> graph, int n, DiscreteDistributionSampler sampler,
			Map<Pair<Long>, EdgeData> sampleEdges, Set<Long> vertices, int categoryIdx, int numCategories,
			int timePeriod) {
		int[] edgeIndices = sampler.sample1DUniformWithoutReplacement(graph.getEdgeCount(), n);
		int curEdgeIdx = 0, i = 0;
		for (E e : graph.getEdges()) {
			if (edgeIndices[i] == curEdgeIdx) {
				double magnitude = Double.NaN;
				if (e instanceof WeightedEdge)
					magnitude = (Double) ((WeightedEdge<?>) e).weight;
				addEdgeToSample(graph.getEndpoints(e), sampleEdges, vertices, categoryIdx, numCategories, magnitude,
						timePeriod);
				if (++i >= n)
					break;
			}
			curEdgeIdx++;
		}
	}

	private static DirectedSparseGraph<Long, Integer> buildCommunicationNetwork(SocialMediaDao mediaDao,
			Set<Long> userIds, MessageLoader loader, Date startDate, Date endDate) throws Exception {
		DirectedSparseGraph<Long, Integer> graph = GraphCreator.createNodeGraph(userIds);
		GraphCreator.applyClassifier(loader, graph);
		for (long id : userIds)
			loader.loadMessages(id, startDate, endDate);
		return graph;
	}

	private static String[] sampleNetworkGraphs(Map<Pair<Long>, EdgeData> edges, int n,
			DiscreteDistributionSampler sampler, SocialMediaDao mediaDao, Set<Long> userIds, MessageLoader loader,
			Date startDate, Date endDate) throws Exception {
		boolean testMutualNetworksOnly = cfg.getLocalBooleanProperty(InfluenceMeasurement.PROP_MUTUAL_ONLY);
		int numCategories = testMutualNetworksOnly ? 2 : 4;
		String[] categoryNames = new String[numCategories];
		int idx = 0;

		// follower graph
		DirectedSparseGraph<Long, Integer> network = GraphCreator.createFollowerGraph(mediaDao, userIds);
		if (!testMutualNetworksOnly) {
			sampleGraph(network, n, sampler, edges, null, idx, numCategories, 0);
			categoryNames[idx++] = ModelVariant.NetworkType.FOLLOWING.toString();
		}
		sampleGraph(GraphCreator.keepBidirectionalEdges(network), n, sampler, edges, null, idx, numCategories, 0);
		categoryNames[idx++] = ModelVariant.NetworkType.MUTUAL_FOLLOWING.toString();

		// communication graph
		network = buildCommunicationNetwork(mediaDao, userIds, loader, startDate, endDate);
		if (!testMutualNetworksOnly) {
			sampleGraph(network, n, sampler, edges, null, idx, numCategories, 0);
			categoryNames[idx++] = ModelVariant.NetworkType.REPLY.toString();
		}
		sampleGraph(GraphCreator.keepBidirectionalEdges(network), n, sampler, edges, null, idx, numCategories, 0);
		categoryNames[idx++] = ModelVariant.NetworkType.MUTUAL_REPLY.toString();

		return categoryNames;
	}

	private static void findTopEdges(WeightedDirectedGraph<Long, Double> graph, int n,
			Map<Pair<Long>, EdgeData> edges, Set<Long> vertices, int categoryIdx, int numCategories, int timePeriod) {
		PriorityQueue<WeightedEdge<Double>> curTopEdges = new PriorityQueue<WeightedEdge<Double>>(n,
				new Comparator<WeightedEdge<Double>>() {
					@Override
					public int compare(WeightedEdge<Double> o1, WeightedEdge<Double> o2) {
						return Double.compare(o1.weight, o2.weight);
					}
				});
		for (WeightedEdge<Double> e : graph.getEdges()) {
			if (curTopEdges.isEmpty() || (e.weight >= curTopEdges.peek().weight)) {
				if (curTopEdges.size() >= n)
					curTopEdges.poll();
				curTopEdges.add(e);
			}
		}
		for (WeightedEdge<Double> e : curTopEdges)
			addEdgeToSample(graph.getEndpoints(e), edges, vertices, categoryIdx, numCategories, e.weight, timePeriod);
	}

	private static String[] sampleExperimentResultGraphs(Map<Pair<Long>, EdgeData> sampleEdges,
			Map<Pair<Long>, EdgeData> topEdges, Set<Long> vertices, int n, DiscreteDistributionSampler sampler,
			String resultsFileName) throws Exception {
		List<EvaluationResult> results = EvaluateNetworks.EvaluationResult.readResults(resultsFileName);
		File resultsDir = (new File(resultsFileName)).getAbsoluteFile().getParentFile();
		int numCategories = results.size();
		String[] categoryNames = new String[numCategories];
		for (int i = 0; i < numCategories; i++) {
			WeightedDirectedGraph<Long, Double> graph = WeightedDirectedGraph.readInfluenceGraph(new File(resultsDir,
					results.get(i).graphFileName));
			int timePeriod = results.get(i).dataVariant.getTimePeriodLength();
			sampleGraph(graph, n, sampler, sampleEdges, vertices, i, numCategories, timePeriod);
			findTopEdges(graph, n, topEdges, vertices, i, numCategories, timePeriod);
			categoryNames[i] = Integer.toString(i);
		}
		return categoryNames;
	}

	private static void writeMessage(PrintWriter w, Map<Long, String> screenNames, Pair<Long> e,
			SocialMediaMessage msg) {
		String text = msg.getText().replace(';', '_').replace('\n', '_').replace('\r', '_');
		String senderName = screenNames.get(msg.getSenderId());
		if (senderName == null)
			senderName = Long.toString(msg.getSenderId());
		String senderSrc = (msg.getSenderId() == e.getFirst()) ? senderName : "";
		String senderDst = (msg.getSenderId() == e.getSecond()) ? senderName : "";
		w.println(senderSrc + ";" + senderDst + ";" + msg.getDate() + ";" + text);
	}

	private static List<SocialMediaMessage> sortByDate(List<SocialMediaMessage> messages) {
		Collections.sort(messages, new Comparator<SocialMediaMessage>() {
			@Override
			public int compare(SocialMediaMessage o1, SocialMediaMessage o2) {
				return o1.getDate().compareTo(o2.getDate());
			}
		});
		return messages;
	}

	private static void writeFilteredConversation(PrintWriter w, MessageLoader loader, Map<Long, String> screenNames,
			Pair<Long> e, Date startDate, Date endDate, int timePeriod) throws Exception {
		Date curStartDate = startDate;
		GregorianCalendar cal = new GregorianCalendar();
		cal.setTime(startDate);
		cal.add(Calendar.DAY_OF_YEAR, timePeriod);
		int idx = 0;
		while (!cal.getTime().after(endDate)) {
			Date curEndDate = cal.getTime();

			w.println("time period " + idx + ", influencee;;;");
			List<SocialMediaMessage> messages = new ArrayList<SocialMediaMessage>();
			// We want to print the same messages that would be visible to the influence model, e.g. discard retweets,
			// so we go through the medium-specific message loader instead of calling getRawMessage. The result will be
			// ordered by recipient, so it has to be re-sorted by date.
			loader.loadMessages(e.getSecond(), curStartDate, curEndDate, messages);
			for (SocialMediaMessage msg : sortByDate(messages))
				writeMessage(w, screenNames, e, msg);

			w.println("time period " + idx + ", influencer;;;");
			messages = new ArrayList<SocialMediaMessage>();
			loader.loadMessages(e.getFirst(), curStartDate, curEndDate, messages);
			for (SocialMediaMessage msg : sortByDate(messages))
				writeMessage(w, screenNames, e, msg);

			curStartDate = curEndDate;
			cal.add(Calendar.DAY_OF_YEAR, timePeriod);
			idx++;
		}
	}

	private static void writeFullConversation(PrintWriter w, MessageLoader loader, Map<Long, String> screenNames,
			Pair<Long> e, Date startDate, Date endDate) throws Exception {
		List<SocialMediaMessage> messages = new ArrayList<SocialMediaMessage>();
		messages.addAll(loader.getRawMessages(e.getFirst(), startDate, endDate));
		messages.addAll(loader.getRawMessages(e.getSecond(), startDate, endDate));
		for (SocialMediaMessage msg : sortByDate(messages))
			writeMessage(w, screenNames, e, msg);
	}

	private static void writeConversations(Map<Pair<Long>, EdgeData> edges, MessageLoader loader,
			Map<Long, String> screenNames, Date startDate, Date endDate, boolean annotateResults, String fileName)
					throws Exception {
		PrintWriter w = new PrintWriter(fileName, "UTF-8");
		try {
			w.println("sender ID (influencer);sender ID (influencee);date;content");
			for (Map.Entry<Pair<Long>, EdgeData> e : edges.entrySet()) {
				if (annotateResults) {
					for (int timePeriod : e.getValue().timePeriods)
						writeFilteredConversation(w, loader, screenNames, e.getKey(), startDate, endDate, timePeriod);
				} else
					writeFullConversation(w, loader, screenNames, e.getKey(), startDate, endDate);
				w.println(";;;");
			}
		} finally {
			w.close();
		}
	}

	private static void writeAnnotationSheet(Map<Pair<Long>, EdgeData> edges, String[] categoryNames,
			Map<Long, String> screenNames, boolean printMagnitude, String fileName) throws IOException {
		PrintWriter w = new PrintWriter(fileName);
		try {
			w.print("fromId;toId;fromName;toName;");
			for (String name : categoryNames) {
				w.print(name + ";");
				if (printMagnitude)
					w.print("m" + name + ";");
			}
			w.println("annotation");
			for (Map.Entry<Pair<Long>, EdgeData> e : edges.entrySet()) {
				Long v1 = e.getKey().getFirst();
				String v1Name = screenNames.get(v1);
				if (v1Name == null)
					v1Name = v1.toString();
				Long v2 = e.getKey().getSecond();
				String v2Name = screenNames.get(v2);
				if (v2Name == null)
					v2Name = v2.toString();
				w.print(v1 + ";" + v2 + ";" + v1Name + ";" + v2Name + ";");
				for (int i = 0; i < e.getValue().isMember.length; i++) {
					w.print((e.getValue().isMember[i] ? "1" : "0") + ";");
					if (printMagnitude)
						w.print(e.getValue().magnitude[i] + ";");
				}
				w.println();
			}
		} finally {
			w.close();
		}
	}

	public static void main(String[] args) throws Exception {
		int sampleSize = defaultSampleSize;
		String resultsFileName = null;
		if (args.length > 0) {
			sampleSize = Integer.parseInt(args[0]);
			if (args.length > 1)
				resultsFileName = args[1];
		}

		Date endDate = cfg.getDateProperty(ExperimentConfiguration.PROP_END_DATE);
		IntervalCalculator ic = new IntervalCalculator(endDate);
		Date startDate = ic.getStartDate();

		SocialMediaDao mediaDao = SocialMediaDaoFactory.createDao();
		Set<Long> userIds = mediaDao.getUserIds(true);
		Map<Long, String> screenNames = mediaDao.getScreenNames(userIds);

		MessageLoader loader;
		Map<Pair<Long>, EdgeData> sampleEdges = new HashMap<Pair<Long>, EdgeData>();
		Map<Pair<Long>, EdgeData> topEdges = new HashMap<Pair<Long>, EdgeData>();
		String[] categoryNames;
		boolean annotateResults = false;
		DiscreteDistributionSampler sampler = new DiscreteDistributionSampler();
		if (resultsFileName == null) {
			loader = MessageLoaderFactory.createMessageLoader(mediaDao, userIds);
			loader.cacheMessages(startDate, endDate);
			categoryNames = sampleNetworkGraphs(sampleEdges, sampleSize, sampler, mediaDao, userIds, loader, startDate,
					endDate);
		} else {
			Set<Long> trueUserIds = new HashSet<Long>();
			categoryNames = sampleExperimentResultGraphs(sampleEdges, topEdges, trueUserIds, sampleSize, sampler,
					resultsFileName);
			loader = MessageLoaderFactory.createMessageLoader(mediaDao, trueUserIds);
			loader.cacheMessages(startDate, endDate);
			annotateResults = true;
		}

		// Write "conversations" (temporally ordered tweets of two users corresponding to a sampled edge). Generate an
		// annotation sheet containing the user IDs and their network membership.
		writeConversations(sampleEdges, loader, screenNames, startDate, endDate, annotateResults,
				"conversations-sample.csv");
		writeAnnotationSheet(sampleEdges, categoryNames, screenNames, annotateResults, "annotation-sample.csv");
		if (!topEdges.isEmpty()) {
			writeConversations(topEdges, loader, screenNames, startDate, endDate, annotateResults,
					"conversations-top.csv");
			writeAnnotationSheet(topEdges, categoryNames, screenNames, annotateResults, "annotation-top.csv");
		}
		System.err.println("done");
	}

}
