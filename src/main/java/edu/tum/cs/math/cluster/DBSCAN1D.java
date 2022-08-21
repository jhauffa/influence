package edu.tum.cs.math.cluster;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.ml.clustering.Cluster;
import org.apache.commons.math3.ml.clustering.Clusterable;
import org.apache.commons.math3.ml.clustering.DBSCANClusterer;
import org.apache.commons.math3.ml.distance.ManhattanDistance;

import edu.tum.cs.nlp.corpus.TokenizedMessage;

public class DBSCAN1D {

	private static final long refDate;
	static {
		try {
			SimpleDateFormat sdf = new SimpleDateFormat("yyyy.MM.dd HH:mm:ss ZZZ");
			refDate = sdf.parse("2010.01.01 00:00:00 UTC").getTime();
		} catch (ParseException ex) {
			throw new RuntimeException(ex);
		}
	}

	public static class MessageClusterPoint implements Clusterable, Comparable<MessageClusterPoint> {
		private final TokenizedMessage msg;
		private final int date;

		public MessageClusterPoint(TokenizedMessage msg) {
			this.msg = msg;
			date = (int) ((msg.getDate().getTime() - refDate) / 1000 / 60);
		}

		@Override
		public double[] getPoint() {
			return new double[] { date };
		}

		public TokenizedMessage getMessage() {
			return msg;
		}

		@Override
		public int compareTo(MessageClusterPoint q) {
			int timeCompare = Integer.compare(date, q.date);
			if (timeCompare == 0)
				return Integer.compare(msg.hashCode(), q.msg.hashCode());
			return timeCompare;
		}
	}

	/**
	 * Cluster messages with DBSCAN and distribute the noise to the nearest clusters.
	 */
	public static List<List<TokenizedMessage>> clusterMessages(List<TokenizedMessage> msgs, double epsilon,
			int minPts) {
		// perform clustering via Apache Math implementation of DBSCAN
		DBSCANClusterer<MessageClusterPoint> clusterer = new DBSCANClusterer<MessageClusterPoint>(epsilon, minPts - 1,
				new ManhattanDistance());
		Set<MessageClusterPoint> points = new HashSet<MessageClusterPoint>(msgs.size());
		for (TokenizedMessage msg : msgs)
			points.add(new MessageClusterPoint(msg));
		List<Cluster<MessageClusterPoint>> clusters = clusterer.cluster(points);

		// all points not assigned to a cluster are noise - turn each point not assigned to a cluster into a cluster of
		// its own
		for (Cluster<MessageClusterPoint> cluster : clusters)
			points.removeAll(cluster.getPoints());
		List<List<MessageClusterPoint>> mergedClusters =
				new ArrayList<List<MessageClusterPoint>>(clusters.size() + points.size());
		for (Cluster<MessageClusterPoint> cluster : clusters)
			mergedClusters.add(new ArrayList<MessageClusterPoint>(cluster.getPoints()));
		for (MessageClusterPoint noise : points) {
			List<MessageClusterPoint> noiseCluster = new ArrayList<MessageClusterPoint>(1);
			noiseCluster.add(noise);
			mergedClusters.add(noiseCluster);
		}

		// ensure that clusters are sorted by date
		for (List<MessageClusterPoint> cluster : mergedClusters)
			Collections.sort(cluster);

		// convert clusters to lists of messages
		List<List<TokenizedMessage>> clusteredMsgs = new ArrayList<List<TokenizedMessage>>(mergedClusters.size());
		for (List<MessageClusterPoint> cluster : mergedClusters) {
			List<TokenizedMessage> msgList = new ArrayList<TokenizedMessage>(cluster.size());
			for (MessageClusterPoint p : cluster)
				msgList.add(p.getMessage());
			clusteredMsgs.add(msgList);
		}
		return clusteredMsgs;
	}

}
