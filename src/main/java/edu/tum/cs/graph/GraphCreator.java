package edu.tum.cs.graph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.tum.cs.db.SocialGraphDao;
import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.db.entities.SocialMediaUser;
import edu.tum.cs.db.loader.MessageLoader;
import edu.tum.cs.math.dist.DiscreteDistributionSampler;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

public class GraphCreator {

	public static DirectedGraph<Long, Integer> createGraph(boolean communicationNetwork) throws Exception {
		DirectedGraph<Long, Integer> graph = new DirectedSparseGraph<Long, Integer>();
		int edgeIdx = 0;
		SocialGraphDao dao = new SocialGraphDao();
		Set<Long> userIds = dao.getUserIds();
		for (Long userId : userIds) {
			Map<Long, SocialEdge> outgoingEdges = dao.getOutgoingSocialEdges(userId);
			edgeIdx = addUser(graph, edgeIdx, userId, userIds, outgoingEdges, communicationNetwork);
		}
		return graph;
	}

	public static DirectedGraph<Long, Integer> createGraph(Set<Long> userIds,
			Map<Long, Map<Long, SocialEdge>> outgoingEdgesMap, boolean communicationNetwork) {
		DirectedGraph<Long, Integer> graph = new DirectedSparseGraph<Long, Integer>();
		int edgeIdx = 0;
		for (Long userId : userIds) {
			Map<Long, SocialEdge> outgoingEdges = outgoingEdgesMap.get(userId);
			edgeIdx = addUser(graph, edgeIdx, userId, userIds, outgoingEdges, communicationNetwork);
		}
		return graph;
	}

	private static int addUser(DirectedGraph<Long, Integer> graph, int edgeIdx, long userId, Set<Long> userIds,
			Map<Long, SocialEdge> outgoingEdges, boolean communicationNetwork) {
		graph.addVertex(userId);
		for (SocialEdge edge : outgoingEdges.values()) {
			if (((communicationNetwork && edge.hasTheta()) || (!communicationNetwork && edge.isExplicit)) &&
					userIds.contains(edge.toUserId)) {
				graph.addVertex(edge.toUserId);
				graph.addEdge(edgeIdx++, userId, edge.toUserId, EdgeType.DIRECTED);
			}
		}
		return edgeIdx;
	}

	public static DirectedSparseGraph<Long, Integer> createNodeGraph(Set<Long> userIds) {
		DirectedSparseGraph<Long, Integer> graph = new DirectedSparseGraph<Long, Integer>();
		for (Long userId : userIds)
			graph.addVertex(userId);
		return graph;
	}

	// Note that in the directed graph generated by this method, edges are oriented in the opposite direction of
	// following, i.e., the followee points to the follower.
	public static DirectedSparseGraph<Long, Integer> createFollowerGraph(SocialMediaDao mediaDao, Set<Long> userIds)
			throws Exception {
		DirectedSparseGraph<Long, Integer> network = createNodeGraph(userIds);
		int idx = 0;
		for (Long userId : userIds) {
			SocialMediaUser user = mediaDao.getUser(userId);
			user.getOutgoingEdges().retainAll(userIds);	// outgoing edges == people I receive information from
			for (Long otherUserId : user.getOutgoingEdges())
				network.addEdge(idx++, otherUserId, userId, EdgeType.DIRECTED);
			user.getIncomingEdges().retainAll(userIds);
			for (Long otherUserId : user.getIncomingEdges())
				network.addEdge(idx++, userId, otherUserId, EdgeType.DIRECTED);
		}
		return network;
	}

	private static class MessageClassifier implements MessageLoader.MessageTypeCallback {
		private final DirectedSparseGraph<Long, Integer> networkGraph;
		private int idx = 0;

		public MessageClassifier(DirectedSparseGraph<Long, Integer> networkGraph) {
			this.networkGraph = networkGraph;
		}

		@Override
		public boolean foundNonAddressiveMessage(SocialMediaMessage message, long senderId,
				Collection<Long> exposedIds) {
			return true;
		}

		@Override
		public boolean foundAddressiveMessage(SocialMediaMessage message, long senderId,
				Collection<Long> recipientIds) {
			for (Long recipientId : recipientIds) {
				if (networkGraph.containsVertex(recipientId))	// do not add vertices for users not in the userId list
					networkGraph.addEdge(idx++, senderId, recipientId, EdgeType.DIRECTED);
			}
			return true;
		}

		@Override
		public boolean foundSharedMessage(SocialMediaMessage message, long senderId, long orginalSenderId) {
			return false;
		}
	}

	public static void applyClassifier(MessageLoader loader, DirectedSparseGraph<Long, Integer> networkGraph) {
		loader.setCallback(new MessageClassifier(networkGraph));
	}

	public static DirectedSparseGraph<Long, Integer> keepBidirectionalEdges(DirectedSparseGraph<Long, Integer> graph) {
		DirectedSparseGraph<Long, Integer> tfGraph = new DirectedSparseGraph<Long, Integer>();
		for (Long v : graph.getVertices())
			tfGraph.addVertex(v);
		for (Integer e : graph.getEdges()) {
			Pair<Long> v = graph.getEndpoints(e);
			if (graph.findEdge(v.getSecond(), v.getFirst()) != null)
				tfGraph.addEdge(e, v.getFirst(), v.getSecond(), EdgeType.DIRECTED);
		}
		return tfGraph;
	}

	/**
	 * Construct a graph by random sampling without replacement from the complement graph. The complement graph of G is
	 * constructed by removing the edges of G from a complete graph with the same vertices as G. Explicitly storing all
	 * edges of the complement graph, and then sampling from that set, takes a prohibitively large amount of memory, so
	 * we sample the indices of the complement graphs's edges.
	 */
	public static DirectedSparseGraph<Long, Integer> sampleComplementGraph(DirectedSparseGraph<Long, Integer> graph,
			int numEdges) {
		// sample edge indices of complement graph uniformly and without replacement
		int numVertices = graph.getVertexCount();
		int numComplementEdges = ((numVertices - 1) * numVertices) - graph.getEdgeCount();
		DiscreteDistributionSampler sampler = new DiscreteDistributionSampler();
		int[] compEdgeIndices = sampler.sample1DUniformWithoutReplacement(numComplementEdges, numEdges);  // asc. order

		// realize sample of complement graph by enumerating its edges
		List<Long> vertices = new ArrayList<Long>(graph.getVertices());	// ensure stable iteration order
		DirectedSparseGraph<Long, Integer> compGraph = new DirectedSparseGraph<Long, Integer>();
		int n = 0;
		int compEdgeIdx = 0;
		for (Long v1 : vertices) {
			if (n == numEdges)
				break;

			int maxIdx = compEdgeIdx + (numVertices - graph.outDegree(v1) - 1);
			for (Long v2 : vertices) {
				if ((v1 != v2) && (graph.findEdge(v1, v2) == null)) {
					if (compEdgeIndices[n] == compEdgeIdx) {
						compGraph.addEdge(n, new Pair<Long>(v1, v2), EdgeType.DIRECTED);
						n++;
						if (n == numEdges)
							break;
						if (compEdgeIndices[n] >= maxIdx) {
							compEdgeIdx = maxIdx;
							break;
						}
					}
					compEdgeIdx++;
				}
			}
		}
		return compGraph;
	}

}