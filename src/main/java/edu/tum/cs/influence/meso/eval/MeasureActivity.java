package edu.tum.cs.influence.meso.eval;

import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import edu.tum.cs.db.SocialGraphDao;
import edu.tum.cs.graph.SocialEdge;
import edu.tum.cs.graph.SocialNode;
import edu.tum.cs.influence.meso.InfluenceExperimentsAnova;
import edu.tum.cs.time.IntervalCalculator;
import edu.tum.cs.time.TimeInterval;
import edu.tum.cs.util.ExperimentConfiguration;

public class MeasureActivity {

	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(InfluenceExperimentsAnova.class);

	public static void main(String[] args) throws Exception {
		IntervalCalculator ic = new IntervalCalculator(cfg.getDateProperty(ExperimentConfiguration.PROP_END_DATE));
		Date[] dates = ic.getCenterDates();
		int[] periods = cfg.getIntListProperty(ExperimentConfiguration.PROP_TIME_PERIODS);

		SocialGraphDao dao = new SocialGraphDao();
		Set<Long> userIds = dao.getUserIds();
		activity(dao, dates, periods, userIds, false);
		activity(dao, dates, periods, userIds, true);
	}

	private static boolean updateStatistics(SocialNode node, SocialEdge edge, SummaryStatistics stats,
			SummaryStatistics statsActive, Date[] dates, int curPeriod) {
		int numPeriods = 0, numActivePeriods = 0;
		boolean active = false;
		for (Date curDate : dates) {
			double[] pastTheta, futureTheta;
			if (node != null) {
				pastTheta = node.getTheta(new TimeInterval(curDate, -curPeriod));
				futureTheta = node.getTheta(new TimeInterval(curDate, curPeriod));
			} else {
				pastTheta = edge.getTheta(new TimeInterval(curDate, -curPeriod));
				futureTheta = edge.getTheta(new TimeInterval(curDate, curPeriod));
			}

			if ((pastTheta != null) && (futureTheta != null))
				numActivePeriods++;
			active |= (pastTheta != null) || (futureTheta != null);
			numPeriods++;
		}

		double r = (double) numActivePeriods / numPeriods;
		stats.addValue(r);
		if (numActivePeriods > 0)
			statsActive.addValue(r);
		return active;
	}

	private static void printStats(Map<Integer, SummaryStatistics> stats, Map<Integer, Long> numActive,
			boolean analyzeEdges) {
		String type = (analyzeEdges ? " edges" : " nodes");
		for (Map.Entry<Integer, SummaryStatistics> e : stats.entrySet()) {
			long s = e.getValue().getN();
			System.out.println("length " + e.getKey() + ": " + (e.getValue().getMean() * 100) +	"% (sd = " +
					(e.getValue().getStandardDeviation() * 100) + "%) of " + s + type + " active (paired)");

			if (numActive != null) {
				long n = numActive.get(e.getKey());
				System.out.println(n + " of " + s + type + " active at least once (" + (((double) n / s) * 100) + "%)");
			}
		}
	}

	private static void activity(SocialGraphDao dao, Date[] dates, int[] periods, Set<Long> userIds,
			boolean analyzeEdges) throws Exception {
		Map<Integer, SummaryStatistics> perUserStats = new HashMap<Integer, SummaryStatistics>();
		Map<Integer, SummaryStatistics> perUserStatsActive = new HashMap<Integer, SummaryStatistics>();
		Map<Integer, Long> numActive = new HashMap<Integer, Long>();

		Map<Long, SocialNode> socialNodes = null;
		if (!analyzeEdges)
			socialNodes = dao.getSocialNodes(userIds);

		for (int curPeriod : periods) {
			SummaryStatistics stats = new SummaryStatistics();
			perUserStats.put(curPeriod, stats);
			SummaryStatistics statsActive = new SummaryStatistics();
			perUserStatsActive.put(curPeriod, statsActive);

			long n = 0;
			for (Long userId : userIds) {
				if (analyzeEdges) {
					Map<Long, SocialEdge> outgoingEdges = dao.getOutgoingSocialEdges(userId);
					for (SocialEdge edge : outgoingEdges.values()) {
						n += updateStatistics(null, edge, stats, statsActive, dates, curPeriod) ? 1 : 0;
					}
				} else {
					SocialNode node = socialNodes.get(userId);
					n += updateStatistics(node, null, stats, statsActive, dates, curPeriod) ? 1 : 0;
				}
			}
			numActive.put(curPeriod, n);
		}

		System.out.println("all: ");
		printStats(perUserStats, numActive, analyzeEdges);
		System.out.println("active only: ");
		printStats(perUserStatsActive, null, analyzeEdges);
		System.out.println();
	}

}
