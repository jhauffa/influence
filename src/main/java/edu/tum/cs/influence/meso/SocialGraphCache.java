package edu.tum.cs.influence.meso;

import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import edu.tum.cs.db.SocialGraphDao;
import edu.tum.cs.graph.SocialEdge;
import edu.tum.cs.graph.SocialNode;
import edu.tum.cs.graph.clique.VertexSubgraphStats;
import edu.tum.cs.time.TimeInterval;
import edu.tum.cs.util.io.Serializer;

public class SocialGraphCache {

	private final Map<Long, SocialNode> nodes;
	private final Map<Long, Map<Long, SocialEdge>> allOutgoingEdges;
	private final SocialNode globalNode;

	private final Map<Long, Set<Long>> cliques = new HashMap<Long, Set<Long>>();
	private final Map<Long, Collection<Collection<Long>>> communitiesPerc =
			new HashMap<Long, Collection<Collection<Long>>>();
	private final Map<Long, Collection<Collection<Long>>> communitiesEdge =
			new HashMap<Long, Collection<Collection<Long>>>();
	private final Map<Long, Long> numCliques = new HashMap<Long, Long>();
	private final Map<Long, Long> numCommunitiesPerc = new HashMap<Long, Long>();
	private final Map<Long, Long> numCommunitiesEdge = new HashMap<Long, Long>();

	private final boolean isExplicit;
	private final boolean isSymmetric;

	public SocialGraphCache(File cliquesFile, File commPercFile, File commEdgeFile) throws Exception {
		SocialGraphDao dao = new SocialGraphDao();
		Set<Long> userIds = dao.getUserIds();
		nodes = dao.getSocialNodes(userIds);
		globalNode = dao.getSocialNode(dao.getGlobalNodeId());
		isExplicit = dao.hasExplicitEdges();
		isSymmetric = dao.isExplicitSymmetric();

		allOutgoingEdges = new HashMap<Long, Map<Long, SocialEdge>>();
		for (SocialNode node : nodes.values()) {
			// construct the subgraph induced by the selected nodes
			Map<Long, SocialEdge> outgoingEdges = dao.getOutgoingSocialEdges(node.id);
			outgoingEdges.keySet().retainAll(userIds);
			allOutgoingEdges.put(node.id, outgoingEdges);

			// If there is no explicit social network graph, set the node incidence counts according to the
			// communication graph.
			if (!isExplicit) {
				node.reportedNumOutgoing = outgoingEdges.size();
				for (long targetId : outgoingEdges.keySet())
					nodes.get(targetId).reportedNumIncoming++;
			}
		}

		if (cliquesFile != null)
			prepareCliques(cliquesFile, userIds, cliques, numCliques);
		if (commPercFile != null)
			prepareCommunities(commPercFile, userIds, communitiesPerc, numCommunitiesPerc);
		if (commEdgeFile != null)
			prepareCommunities(commEdgeFile, userIds, communitiesEdge, numCommunitiesEdge);
	}

	private static void prepareCliques(File f, Collection<Long> userIds, Map<Long, Set<Long>> cliqueMap,
			Map<Long, Long> numCliqueMap) {
		Map<Long, VertexSubgraphStats<Long>> subgraphStats = Serializer.loadObjectFromFile(f);
		if (subgraphStats == null)
			throw new RuntimeException("Error loading cliques from file '" + f.getPath() + "'");
		for (Long userId : userIds) {
			VertexSubgraphStats<Long> pcl = subgraphStats.get(userId);
			if ((pcl != null) && (pcl.maxSizeSubgraph != null)) {
				cliqueMap.put(userId, new HashSet<Long>(pcl.maxSizeSubgraph));
				numCliqueMap.put(userId, pcl.numSubgraphs);
			}
		}
	}

	private static void prepareCommunities(File f, Collection<Long> userIds,
			Map<Long, Collection<Collection<Long>>> communityMap, Map<Long, Long> numCommunityMap) {
		Collection<Collection<Long>> communities = Serializer.loadObjectFromFile(f);
		if (communities == null)
			throw new RuntimeException("Error loading communities from file '" + f.getPath() + "'");
		for (Collection<Long> community : communities) {
			for (Long userId : community) {
				if (!userIds.contains(userId))
					continue;

				Collection<Collection<Long>> userCommunities = communityMap.get(userId);
				if (userCommunities == null) {
					userCommunities = new HashSet<Collection<Long>>();
					communityMap.put(userId, userCommunities);
				}

				userCommunities.add(community);
				Long numUserCommunities = numCommunityMap.get(userId);
				if (numUserCommunities == null)
					numUserCommunities = 0L;
				numCommunityMap.put(userId, numUserCommunities + 1);
			}
		}
	}

	public boolean isExplicit() {
		return isExplicit;
	}

	public boolean isSymmetric() {
		return isSymmetric;
	}

	public Collection<SocialNode> getNodes() {
		return nodes.values();
	}

	public SocialNode getNode(long userId) {
		return nodes.get(userId);
	}

	public SocialNode getGlobalNode() {
		return globalNode;
	}

	public Map<Long, Map<Long, SocialEdge>> getAllOutgoingEdges() {
		return allOutgoingEdges;
	}

	public Map<Long, SocialEdge> getOutgoingEdges(long userId) {
		return allOutgoingEdges.get(userId);
	}

	public long getNumCliques(long userId) {
		Long n = numCliques.get(userId);
		return (n != null) ? n : 0;
	}

	public long getNumCommunitiesPerc(long userId) {
		Long n = numCommunitiesPerc.get(userId);
		return (n != null) ? n : 0;
	}

	public long getNumCommunitiesEdge(long userId) {
		Long n = numCommunitiesEdge.get(userId);
		return (n != null) ? n : 0;
	}

	public Set<Long> getLargestMaximalClique(long userId) {
		return cliques.get(userId);
	}

	public Collection<Collection<Long>> getCommunitiesPerc(long userId) {
		return communitiesPerc.get(userId);
	}

	public Collection<Collection<Long>> getCommunitiesEdge(long userId) {
		return communitiesEdge.get(userId);
	}

	public Set<Long> getExplicitOutgoingDistanceTwo(Map<Long, SocialEdge> outgoingEdges) {
		Set<Long> explicitDistance2 = new HashSet<>();
		for (SocialEdge edge : outgoingEdges.values()) {
			if (edge.isExplicit) {
				explicitDistance2.add(edge.toUserId);
				Map<Long, SocialEdge> outgoingEdgesDistanceOne = allOutgoingEdges.get(edge.toUserId);
				for (SocialEdge edgeDistanceOne : outgoingEdgesDistanceOne.values()) {
					if (edgeDistanceOne.isExplicit)
						explicitDistance2.add(edgeDistanceOne.toUserId);
				}
			}
		}
		return explicitDistance2;
	}

	public Set<Long> getCommunicationDistanceTwo(Map<Long, SocialEdge> outgoingEdges, TimeInterval iv) {
		Set<Long> commDistTwo = new HashSet<>();
		for (SocialEdge edge : outgoingEdges.values()) {
			Map<Long, SocialEdge> outgoingEdgesDistanceOne = allOutgoingEdges.get(edge.toUserId);
			SocialEdge reverseEdge = outgoingEdgesDistanceOne.get(edge.fromUserId);
			if ((edge.getTheta(iv) != null) || ((reverseEdge != null) && (reverseEdge.getTheta(iv) != null))) {
				commDistTwo.add(edge.toUserId);
				for (SocialEdge edgeDistanceOne : outgoingEdgesDistanceOne.values()) {
					reverseEdge = allOutgoingEdges.get(edgeDistanceOne.toUserId).get(edgeDistanceOne.fromUserId);
					if ((edgeDistanceOne.getTheta(iv) != null) ||
						((reverseEdge != null) && (reverseEdge.getTheta(iv) != null)))
						commDistTwo.add(edgeDistanceOne.toUserId);
				}
			}
		}
		return commDistTwo;
	}

}
