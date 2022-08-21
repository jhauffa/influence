package edu.tum.cs.influence.meso;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;
import java.util.Collection;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import edu.tum.cs.graph.SocialEdge;
import edu.tum.cs.graph.SocialNode;
import edu.tum.cs.math.dist.DiscreteDistribution;
import edu.tum.cs.time.TimeInterval;

public enum WeightFunction {

	NUM_INCOMING {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return candidate.reportedNumIncoming;
		}
	},

	NUM_INCOMING_1P {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return Math.log1p(candidate.reportedNumIncoming);
		}
	},

	NUM_OUTGOING {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return candidate.reportedNumOutgoing;
		}
	},

	NUM_OUTGOING_1P {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return Math.log1p(candidate.reportedNumOutgoing);
		}
	},

	INCOMING_OUTGOING_LOG_RATIO {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			double logNumOutgoing = Math.log1p(candidate.reportedNumOutgoing);
			if (logNumOutgoing == 0.0)
				logNumOutgoing = Math.log1p(1.0);
			return Math.log1p(candidate.reportedNumIncoming) / logNumOutgoing;
		}
	},

	SHARED_RATIO {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			int numNonAddr = Math.max(1, candidate.getNumNonAddr(iv));
			return candidate.getNumGotShared(iv) / numNonAddr;
		}
	},

	NUM_CANDIDATE_GOT_SHARED {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return candidate.getNumGotShared(iv);
		}
	},

	NUM_CLIQUE {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return socialGraph.getNumCliques(candidate.id);
		}
	},

	NUM_CLIQUE_1P {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return Math.log1p(socialGraph.getNumCliques(candidate.id));
		}
	},

	NUM_COMMUNITIES_PERC {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return socialGraph.getNumCommunitiesPerc(candidate.id);
		}
	},

	INVERSE_NUM_COMMUNITIES_PERC {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return 1.0 / (1 + socialGraph.getNumCommunitiesPerc(candidate.id));
		}
	},

	SHARED_COMMUNITIES_PERC {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return countSharedCommunities(socialGraph.getCommunitiesPerc(origin.id), origin.id, candidate.id);
		}
	},

	NUM_COMMUNITIES_EDGE {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return socialGraph.getNumCommunitiesEdge(candidate.id);
		}
	},

	INVERSE_NUM_COMMUNITIES_EDGE {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return 1.0 / (1 + socialGraph.getNumCommunitiesEdge(candidate.id));
		}
	},

	SHARED_COMMUNITIES_EDGE {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return countSharedCommunities(socialGraph.getCommunitiesEdge(origin.id), origin.id, candidate.id);
		}
	},

	SIMILARITY_NODE {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			double[] thetaOrigin = origin.getTheta(iv);
			double[] thetaCandidate = candidate.getTheta(iv);
			if ((thetaOrigin == null) || (thetaCandidate == null))
				return 0.0;
			return 1.0 - DiscreteDistribution.distJS2(thetaOrigin, thetaCandidate);
		}
	},

	SIMILARITY_EDGE {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			double[] thetaOrigin = origin.getSenderTheta(iv);
			double[] thetaCandidate = candidate.getSenderTheta(iv);
			if ((thetaOrigin == null) || (thetaCandidate == null))
				return 0.0;
			return 1.0 - DiscreteDistribution.distJS2(thetaOrigin, thetaCandidate);
		}
	},

	SIMILARITY_TO_MAINSTREAM {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			double[] thetaGlobal = socialGraph.getGlobalNode().getTheta(iv);
			double[] thetaCandidate = candidate.getTheta(iv);
			if (thetaCandidate == null)
				return 0.0;
			return 1.0 - DiscreteDistribution.distJS2(thetaGlobal, thetaCandidate);
		}
	},

	DISTANCE_TO_MAINSTREAM {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			double[] thetaGlobal = socialGraph.getGlobalNode().getTheta(iv);
			double[] thetaCandidate = candidate.getTheta(iv);
			if (thetaCandidate == null)
				return 0.0;	// should be 1.0 if we wanted a true complement of SIMILARITY_TO_MAINSTREAM
			return DiscreteDistribution.distJS2(thetaGlobal, thetaCandidate);
		}
	},

	CANDIDATE_NON_ADDR {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return candidate.getNumNonAddr(iv);
		}
	},

	INVERSE_CANDIDATE_NON_ADDR {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return 1.0 / (1 + candidate.getNumNonAddr(iv));
		}
	},

	CANDIDATE_ADDR {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return ((incomingEdge != null) ? incomingEdge.getNumAddr(iv) : 0) +
					((outgoingEdge != null) ? outgoingEdge.getNumAddr(iv) : 0);
		}
	},

	INVERSE_CANDIDATE_ADDR {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			int n = ((incomingEdge != null) ? incomingEdge.getNumAddr(iv) : 0) +
					((outgoingEdge != null) ? outgoingEdge.getNumAddr(iv) : 0);
			return 1.0 / (1 + n);
		}
	},

	SHARED_NEIGHBORS {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			Map<Long, SocialEdge> candidateNeighbors = socialGraph.getOutgoingEdges(candidate.id);
			int numSharedNeighbors = 0;
			for (SocialEdge edge : socialGraph.getOutgoingEdges(origin.id).values()) {
				if (edge.isExplicit || !socialGraph.isExplicit()) {
					SocialEdge candidateEdge = candidateNeighbors.get(edge.toUserId);
					if ((candidateEdge != null) && (candidateEdge.isExplicit || !socialGraph.isExplicit()))
						numSharedNeighbors++;
				}
			}
			return numSharedNeighbors;
		}
	},

	// focus level is computed via taking the sum of all the probability values in the 90th percentile
	FOCUS_LEVEL {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			double[] thetaNode = candidate.getTheta(iv);
			if (thetaNode == null)
				return 0.0;

			thetaNode = thetaNode.clone();
			Arrays.sort(thetaNode);
			double sumGreater90thPercentile = 0.0;
			for (int i = thetaNode.length - 1; i > thetaNode.length - 1 - (thetaNode.length / 10); i--)
				sumGreater90thPercentile += thetaNode[i];
			// Sum can take on values between 0.1 (uniform distribution) and 1.0 (all probability mass above the 90th
			// percentile); scale to the interval [0,1].
			return (sumGreater90thPercentile - 0.1) / 0.9;
		}
	},

	PAGE_RANK {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return candidate.pageRank;
		}
	},

	BETWEENNESS_CENTRALITY {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return candidate.betweennessCentrality;
		}
	},

	// strong interaction with indicator function: mean value is be around 2 for COMMUNICATION_DISTANCE_ONE, but close
	// to 1 for all other indicators; disable for any experiment that involves separate analysis of indicator and weight
	BIDIRECTIONAL {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			// edges can be null if not present
			if ((incomingEdge != null) && (outgoingEdge != null)) {
				if (incomingEdge.isExplicit && outgoingEdge.isExplicit) {
					if ((incomingEdge.getTheta(iv) != null) && (outgoingEdge.getTheta(iv) != null))
						return 4.0; // bidirectional explicit edge + bidirectional addressive communication
					if ((incomingEdge.getTheta(iv) != null) || (outgoingEdge.getTheta(iv) != null))
						return 2.0; // bidirectional explicit edge + unidirectional addressive communication
					return 1.0; // bidirectional explicit edge
				}
				if ((incomingEdge.getTheta(iv) != null) && (outgoingEdge.getTheta(iv) != null))
					return 2.0;	// unidirectional explicit edge / no explicit edge, but bidirectional addressive comm.
			}
			return 0.5;
		}
	},

	CONSTANT {
		@Override
		public double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv) {
			return 1.0;
		}
	};

	private static long countSharedCommunities(Collection<Collection<Long>> userCommunities, long originId,
			long candidateId) {
		long num = 0;
		if (userCommunities != null)
			for (Collection<Long> community : userCommunities)
				if (community.contains(candidateId))
					num++;
		return num;
	}

	public abstract double eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge,
			SocialEdge outgoingEdge, SocialGraphCache socialGraph, TimeInterval iv);


	private static class WeightStatistics {
		public SummaryStatistics stats = new SummaryStatistics();
		public int numUsers, numUsersEmpty;
	}

	private static Map<ModelVariant, WeightStatistics> variantStats = new HashMap<ModelVariant, WeightStatistics>();

	public static synchronized void addStats(ModelVariant variant, List<Double> weights) {
		WeightStatistics s = variantStats.get(variant);
		if (s == null) {
			s = new WeightStatistics();
			variantStats.put(variant, s);
		}
		s.numUsers++;
		double weightSum = 0.0;
		for (Double weight : weights) {
			weightSum += weight;
			s.stats.addValue(weight);
		}
		if (weightSum == 0.0)
			s.numUsersEmpty++;
	}

	public static void writeAvgStatsCsv(File f) throws IOException {
		Writer w = new FileWriter(f);
		try {
			w.write(ModelVariant.getCsvHeader() + ";mean;variance;standard deviation;minimum value;maximum value;" +
					"sample size;number of users;number of users with empty neighborhood\n");
			for (Map.Entry<ModelVariant, WeightStatistics> e : variantStats.entrySet()) {
				WeightStatistics s = e.getValue();
				w.write(e.getKey() + ";" + s.stats.getMean() + ";" + s.stats.getVariance() + ";" +
						s.stats.getStandardDeviation() + ";" + s.stats.getMin() + ";" + s.stats.getMax() + ";" +
						s.stats.getN() + ";" + s.numUsers + ";" + s.numUsersEmpty + "\n");
			}
		} finally {
			w.close();
		}
	}

	public static void writeEmptyNeighborhoodStatsCsv(File f) throws IOException {
		Map<IndicatorFunction, SummaryStatistics> indicatorStats = new HashMap<IndicatorFunction, SummaryStatistics>();
		Map<WeightFunction, SummaryStatistics> weightStats = new HashMap<WeightFunction, SummaryStatistics>();

		for (Map.Entry<ModelVariant, WeightStatistics> e : variantStats.entrySet()) {
			SummaryStatistics curIndicatorStats = indicatorStats.get(e.getKey().getIndicator());
			if (curIndicatorStats == null) {
				curIndicatorStats = new SummaryStatistics();
				indicatorStats.put(e.getKey().getIndicator(), curIndicatorStats);
			}
			SummaryStatistics curWeightStats = weightStats.get(e.getKey().getWeight());
			if (curWeightStats == null) {
				curWeightStats = new SummaryStatistics();
				weightStats.put(e.getKey().getWeight(), curWeightStats);
			}

			double missing = (double) e.getValue().numUsersEmpty / e.getValue().numUsers;
			curIndicatorStats.addValue(missing);
			curWeightStats.addValue(missing);
		}

		Writer w = new FileWriter(f);
		try {
			w.write("indicator;weight;missing\n");
			for (Map.Entry<IndicatorFunction, SummaryStatistics> e : indicatorStats.entrySet())
				w.write(e.getKey().name() + ";;" + e.getValue().getMean() + "\n");
			for (Map.Entry<WeightFunction, SummaryStatistics> e : weightStats.entrySet())
				w.write(";" + e.getKey().name() + ";" + e.getValue().getMean() + "\n");
		} finally {
			w.close();
		}
	}

	private static final WeightFunction[] requireExplicit = { BIDIRECTIONAL };
	private static final WeightFunction[] requireAsymmetric = {
		NUM_OUTGOING, NUM_OUTGOING_1P, INCOMING_OUTGOING_LOG_RATIO
	};

	public static EnumSet<WeightFunction> getApplicable(boolean isExplicit, boolean isSymmetric) {
		EnumSet<WeightFunction> weights = EnumSet.allOf(WeightFunction.class);
		if (!isExplicit)
			weights.removeAll(Arrays.asList(requireExplicit));
		if (isSymmetric)
			weights.removeAll(Arrays.asList(requireAsymmetric));
		return weights;
	}

}
