package edu.tum.cs.graph.clique;

import java.io.Serializable;
import java.util.Collection;

public class VertexSubgraphStats<V> implements Serializable {

	private static final long serialVersionUID = -7026984143440265140L;

	public long numSubgraphs;
	public Collection<V> maxSizeSubgraph;

	public VertexSubgraphStats() {
		this(0, null);
	}

	public VertexSubgraphStats(long numSubgraphs, Collection<V> maxSizeSubgraph) {
		this.numSubgraphs = numSubgraphs;
		this.maxSizeSubgraph = maxSizeSubgraph;
	}

}
