package edu.tum.cs.influence.meso.eval;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;
import java.util.Set;

import edu.tum.cs.db.SocialGraphDao;
import edu.tum.cs.graph.SocialEdge;
import edu.tum.cs.graph.SocialNode;
import edu.tum.cs.influence.meso.InfluenceExperimentsAnova;
import edu.tum.cs.math.dist.DiscreteDistribution;
import edu.tum.cs.time.IntervalCalculator;
import edu.tum.cs.time.TimeInterval;
import edu.tum.cs.util.ExperimentConfiguration;

public class EntropyCalculation {

	private static final Logger logger = Logger.getLogger(EntropyCalculation.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(InfluenceExperimentsAnova.class);

	public static void main(String[] args) throws Exception {
		boolean printThetas = false;
		if (args.length > 0)
			printThetas = Boolean.parseBoolean(args[0]);
		String outPath = cfg.getLocalProperty(InfluenceExperimentsAnova.PROP_OUT_PATH, ".");
		IntervalCalculator ic = new IntervalCalculator(cfg.getDateProperty(ExperimentConfiguration.PROP_END_DATE));
		Date[] dates = ic.getCenterDates();

		SocialGraphDao dao = new SocialGraphDao();
		Set<Long> userIds = dao.getUserIds();
		Map<TimeInterval, Map<SocialEdge, Double>> entropyPast = new HashMap<>();
		Map<TimeInterval, Map<SocialEdge, Double>> entropyFuture = new HashMap<>();
		entropyForEdges(dao, dates, userIds, entropyPast, entropyFuture);
		writeDataToCSV(new File(outPath, "entropyEdgePast.csv"), entropyPast, printThetas);
		writeDataToCSV(new File(outPath, "entropyEdgeFuture.csv"), entropyFuture, printThetas);
		entropyForNodes(dao, dates, userIds, entropyPast, entropyFuture);
		writeDataToCSV(new File(outPath, "entropyNodePast.csv"), entropyPast, printThetas);
		writeDataToCSV(new File(outPath, "entropyNodeFuture.csv"), entropyFuture, printThetas);
	}

	private static void incHistogram(Map<TimeInterval, Integer> histo, TimeInterval period) {
		Integer c = histo.get(period);
		if (c == null)
			c = 0;
		histo.put(period, c + 1);
	}

	private static int getCount(Map<TimeInterval, Integer> histo, TimeInterval period) {
		Integer c = histo.get(period);
		return (c != null) ? c : 0;
	}

	private static void entropyForEdges(SocialGraphDao dao, Date[] dates, Set<Long> userIds,
			Map<TimeInterval, Map<SocialEdge, Double>> entropyPast,
			Map<TimeInterval, Map<SocialEdge, Double>> entropyFuture) throws Exception {
		Map<TimeInterval, Integer> numActivePast = new HashMap<>();
		Map<TimeInterval, Integer> numActiveFuture = new HashMap<>();
		Map<TimeInterval, Integer> numEdges = new HashMap<>();

		for (Long userId : userIds) {
			Map<Long, SocialEdge> outgoingEdges = dao.getOutgoingSocialEdges(userId);

			for (int i = 0; i < dates.length; i++) {
				for (int timePeriod : cfg.getIntListProperty(ExperimentConfiguration.PROP_TIME_PERIODS)) {
					TimeInterval curPeriodPast = new TimeInterval(dates[i], -timePeriod);
					TimeInterval curPeriodFuture = new TimeInterval(dates[i], timePeriod);
					Map<SocialEdge, Double> entropyEdgeUserMapPast = entropyPast.get(curPeriodPast);
					if (entropyEdgeUserMapPast == null) {
						entropyEdgeUserMapPast = new HashMap<>();
						entropyPast.put(curPeriodPast, entropyEdgeUserMapPast);
					}
					Map<SocialEdge, Double> entropyEdgeUserMapFuture = entropyFuture.get(curPeriodFuture);
					if (entropyEdgeUserMapFuture == null) {
						entropyEdgeUserMapFuture = new HashMap<>();
						entropyFuture.put(curPeriodFuture, entropyEdgeUserMapFuture);
					}

					for (SocialEdge edge : outgoingEdges.values()) {
						double[] thetaPast = edge.getTheta(curPeriodPast);
						double[] thetaFuture = edge.getTheta(curPeriodFuture);
						if (thetaPast != null) {
							entropyEdgeUserMapPast.put(edge, DiscreteDistribution.entropy(thetaPast));
							incHistogram(numActivePast, curPeriodFuture);
						}
						if (thetaFuture != null) {
							entropyEdgeUserMapFuture.put(edge, DiscreteDistribution.entropy(thetaFuture));
							incHistogram(numActiveFuture, curPeriodFuture);
						}
						incHistogram(numEdges, curPeriodFuture);
					}
				}
			}
		}

		for (int i = 0; i < dates.length; i++) {
			for (int timePeriod : cfg.getIntListProperty(ExperimentConfiguration.PROP_TIME_PERIODS)) {
				TimeInterval curPeriod = new TimeInterval(dates[i], timePeriod);
				int n = getCount(numEdges, curPeriod);
				double activePast = ((double) getCount(numActivePast, curPeriod) / n) * 100;
				double activeFuture = ((double) getCount(numActiveFuture, curPeriod) / n) * 100;
				logger.info("day " + dates[i] + ", length " + timePeriod + ": " +
						activePast + "% of past edges active, " +
						activeFuture + "% of future edge active (of " + n + ")");
			}
		}
	}

	private static void entropyForNodes(SocialGraphDao dao, Date[] dates, Set<Long> userIds,
			Map<TimeInterval, Map<SocialEdge, Double>> entropyPast,
			Map<TimeInterval, Map<SocialEdge, Double>> entropyFuture) throws Exception {
		Map<Long, SocialNode> socialNodes = dao.getSocialNodes(userIds);
		for (int i = 0; i < dates.length; i++) {
			for (int timePeriod : cfg.getIntListProperty(ExperimentConfiguration.PROP_TIME_PERIODS)) {
				TimeInterval curPeriodPast = new TimeInterval(dates[i], -timePeriod);
				TimeInterval curPeriodFuture = new TimeInterval(dates[i], timePeriod);
				Map<SocialEdge, Double> entropyNodeUserMapPast = new HashMap<>();
				Map<SocialEdge, Double> entropyNodeUserMapFuture = new HashMap<>();
				int numActivePast = 0;
				int numActiveFuture = 0;
				for (Entry<Long, SocialNode> entry : socialNodes.entrySet()) {
					SocialNode node = entry.getValue();
					Long userId = entry.getKey();
					SocialEdge fakeEdge = new SocialEdge(userId, userId);
					double[] thetaPast = node.getTheta(curPeriodPast);
					double[] thetaFuture = node.getTheta(curPeriodFuture);
					fakeEdge.setTheta(curPeriodPast, thetaPast);
					fakeEdge.setTheta(curPeriodFuture, thetaFuture);
					if (thetaPast != null) {
						entropyNodeUserMapPast.put(fakeEdge, DiscreteDistribution.entropy(thetaPast));
						numActivePast++;
					}
					if (thetaFuture != null) {
						entropyNodeUserMapFuture.put(fakeEdge, DiscreteDistribution.entropy(thetaFuture));
						numActiveFuture++;
					}
				}
				entropyPast.put(curPeriodPast, entropyNodeUserMapPast);
				entropyFuture.put(curPeriodFuture, entropyNodeUserMapFuture);

				double activePast = ((double) numActivePast / socialNodes.size()) * 100;
				double activeFuture = ((double) numActiveFuture / socialNodes.size()) * 100;
				logger.info("day " + dates[i] + ", length " + timePeriod + ": " +
						activePast + "% of past nodes active, " +
						activeFuture + "% of future nodes active (of " + socialNodes.size() + ")");
			}
		}
	}

	private static void writeDataToCSV(File f, Map<TimeInterval, Map<SocialEdge, Double>> entropyEdgeTimeSlice,
			boolean printThetas) throws IOException {
		int numTopics = cfg.getIntProperty(ExperimentConfiguration.PROP_NUM_TOPICS);
		BufferedWriter writer = new BufferedWriter(new FileWriter(f));
		try {
			writer.write("timePeriodBase;timePeriodLength;fromUserId;toUserId;entropy");
			if (printThetas) {
				for (int i = 0; i < numTopics; i++)
					writer.write(";t" + (i + 1));
			}
			writer.newLine();

			for (Entry<TimeInterval, Map<SocialEdge, Double>> entry : entropyEdgeTimeSlice.entrySet()) {
				TimeInterval key = entry.getKey();
				Map<SocialEdge, Double> entropyEdgeUserMap = entry.getValue();
				for (Entry<SocialEdge, Double> entryUser : entropyEdgeUserMap.entrySet()) {
					SocialEdge edge = entryUser.getKey();
					double entropy = entryUser.getValue();
					writer.write(key + ";" + edge.fromUserId + ";" + edge.toUserId + ";" + entropy);

					if (printThetas) {
						double[] theta = edge.getTheta(key);
						for (int i = 0; i < theta.length; i++)
							writer.write(";" + theta[i]);
					}
					writer.newLine();
				}
			}
		} finally {
			writer.close();
		}
	}

}
