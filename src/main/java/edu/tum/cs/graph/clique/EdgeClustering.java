package edu.tum.cs.graph.clique;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import edu.tum.cs.math.cluster.SingleLinkageClustering;
import edu.uci.ics.jung.graph.UndirectedGraph;
import edu.uci.ics.jung.graph.util.Pair;

public class EdgeClustering<V, E> extends SingleLinkageClustering<E> {

	private static final long serialVersionUID = -1725557819822706606L;

	private final UndirectedGraph<V, E> graph;

	public EdgeClustering(UndirectedGraph<V, E> graph, boolean parallel) {
		super(new ArrayList<E>(graph.getEdges()), parallel);
		this.graph = graph;
	}

	@Override
	protected double distance(E e1, E e2) {
		Pair<V> v1 = graph.getEndpoints(e1);
		Pair<V> v2 = graph.getEndpoints(e2);
		V i, j;
		if (v1.getFirst().equals(v2.getFirst())) {
			i = v1.getSecond();
			j = v2.getSecond();
		} else if (v1.getFirst().equals(v2.getSecond())) {
			i = v1.getSecond();
			j = v2.getFirst();
		} else if (v1.getSecond().equals(v2.getFirst())) {
			i = v1.getFirst();
			j = v2.getSecond();
		} else if (v1.getSecond().equals(v2.getSecond())) {
			i = v1.getFirst();
			j = v2.getFirst();
		} else	// no common vertex
			return 1.0;

		Set<V> n1 = new HashSet<V>(graph.getNeighbors(i));
		n1.add(i);
		Set<V> n2 = new HashSet<V>(graph.getNeighbors(j));
		n2.add(j);
		n1.retainAll(n2);
		n2.addAll(graph.getNeighbors(i));
		n2.add(i);

		return 1.0 - ((double) n1.size() / n2.size());
	}

	public Collection<Collection<V>> getVertexClusters(int level) {
		Collection<Collection<E>> edgeClusters = getClusters(level);
		Collection<Collection<V>> vertexClusters = new ArrayList<Collection<V>>(edgeClusters.size());
		for (Collection<E> cluster : edgeClusters) {
			Set<V> vertices = new HashSet<V>();
			for (E edge : cluster)
				vertices.addAll(graph.getIncidentVertices(edge));
			vertexClusters.add(vertices);
		}
		return vertexClusters;
	}

	public double getPartitionDensity(int level) {
		double d = 0.0;
		Collection<Collection<E>> clusters = getClusters(level);
		for (Collection<E> cluster : clusters) {
			int numEdges = cluster.size();
			if (numEdges > 1) {
				Set<V> vertices = new HashSet<V>();
				for (E edge : cluster)
					vertices.addAll(graph.getIncidentVertices(edge));
				int numVertices = vertices.size();
				d += numEdges * ((double) (numEdges - (numVertices - 1)) / ((numVertices - 2) * (numVertices - 1)));
			}
		}
		return (2.0 * d) / graph.getEdgeCount();
	}

}
