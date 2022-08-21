package edu.tum.cs.influence.meso;

import java.util.Arrays;
import java.util.Collection;
import java.util.EnumSet;
import java.util.Set;

import edu.tum.cs.graph.SocialEdge;
import edu.tum.cs.graph.SocialNode;
import edu.tum.cs.time.TimeInterval;

public enum IndicatorFunction {

	EXPLICIT_DISTANCE_ONE {
		@Override
		public int eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv, Set<Long> explicitOutgoingDistanceTwo,
				Set<Long> commDistanceTwo) {
			return ((outgoingEdge != null) && (outgoingEdge.isExplicit)) ? 1 : 0;
		}
	},

	EXPLICIT_DISTANCE_TWO {
		@Override
		public int eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv, Set<Long> explicitOutgoingDistanceTwo,
				Set<Long> commDistanceTwo) {
			return (explicitOutgoingDistanceTwo.contains(candidate.id)) ? 1 : 0;
		}
	},

	COMMUNICATION_DISTANCE_ONE {
		@Override
		public int eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv, Set<Long> explicitOutgoingDistanceTwo,
				Set<Long> commDistanceTwo) {
			return (((outgoingEdge != null) && (outgoingEdge.getTheta(iv) != null)) ||
					((incomingEdge != null) && (incomingEdge.getTheta(iv) != null))) ? 1 : 0;
		}
	},

	COMMUNICATION_DISTANCE_TWO {
		@Override
		public int eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv, Set<Long> explicitOutgoingDistanceTwo,
				Set<Long> commDistanceTwo) {
			return (commDistanceTwo.contains(candidate.id)) ? 1 : 0;
		}
	},

	CLIQUE {
		@Override
		public int eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv, Set<Long> explicitOutgoingDistanceTwo,
				Set<Long> commDistanceTwo) {
			Set<Long> clique = socialGraph.getLargestMaximalClique(origin.id);
			if (clique == null)
				return 0;
			return (clique.contains(candidate.id)) ? 1 : 0;
		}
	},

	COMMUNITIES_PERC {
		@Override
		public int eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv, Set<Long> explicitOutgoingDistanceTwo,
				Set<Long> commDistanceTwo) {
			return (isCommunityMember(socialGraph.getCommunitiesPerc(origin.id), origin.id, candidate.id)) ? 1 : 0;
		}
	},

	COMMUNITIES_EDGE {
		@Override
		public int eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
				SocialGraphCache socialGraph, TimeInterval iv, Set<Long> explicitOutgoingDistanceTwo,
				Set<Long> commDistanceTwo) {
			return (isCommunityMember(socialGraph.getCommunitiesEdge(origin.id), origin.id, candidate.id)) ? 1 : 0;
		}
	};

	private static boolean isCommunityMember(Collection<Collection<Long>> userCommunities, long originId,
			long candidateId) {
		if (userCommunities != null)
			for (Collection<Long> community : userCommunities)
				if (community.contains(candidateId))
					return true;
		return false;
	}

	public abstract int eval(SocialNode origin, SocialNode candidate, SocialEdge incomingEdge, SocialEdge outgoingEdge,
			SocialGraphCache socialGraph, TimeInterval iv, Set<Long> explicitOutgoingDistanceTwo,
			Set<Long> commDistanceTwo);

	public boolean needsExplicitDistanceTwo() {
		return (this == EXPLICIT_DISTANCE_TWO);
	}

	public boolean needsCommunicationDistanceTwo() {
		return (this == COMMUNICATION_DISTANCE_TWO);
	}

	private static final IndicatorFunction[] requireExplicit = { EXPLICIT_DISTANCE_ONE, EXPLICIT_DISTANCE_TWO };

	public static EnumSet<IndicatorFunction> getApplicable(boolean isExplicit, boolean isSymmetric) {
		EnumSet<IndicatorFunction> indicators = EnumSet.allOf(IndicatorFunction.class);
		if (!isExplicit)
			indicators.removeAll(Arrays.asList(requireExplicit));
		return indicators;
	}

}
