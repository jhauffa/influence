package edu.tum.cs.graph.clique;

import java.util.ArrayList;
import java.util.Collection;

import org.apache.commons.collections15.BidiMap;
import org.apache.commons.collections15.bidimap.DualHashBidiMap;

import edu.uci.ics.jung.algorithms.util.Indexer;
import edu.uci.ics.jung.graph.Graph;

public class VertexIndex<V> {

	private final BidiMap<V, Integer> index;
	private int nextId = 0;

	public VertexIndex() {
		index = new DualHashBidiMap<V, Integer>();
	}

	public <E> VertexIndex(Graph<V, E> graph) {
		index = Indexer.create(graph.getVertices());
		if (index.size() > Short.MAX_VALUE)
			throw new IllegalArgumentException("too many vertices");
	}

	public short get(V vertex) {
		return index.get(vertex).shortValue();
	}

	public short getOrInsert(V vertex) {
		Integer mappedVertex = index.get(vertex);
		if (mappedVertex == null) {
			mappedVertex = nextId++;
			index.put(vertex, mappedVertex);
		}
		return mappedVertex.shortValue();
	}

	public short[] packVertices(Collection<V> vertices) {
		short[] packedVertices = new short[vertices.size()];
		int idx = 0;
		for (V v : vertices)
			packedVertices[idx++] = index.get(v).shortValue();
		return packedVertices;
	}

	public Collection<V> unpackVertices(Collection<? extends Number> packedVertices) {
		Collection<V> vertices = new ArrayList<V>(packedVertices.size());
		for (Number v : packedVertices)
			vertices.add(index.getKey(v.intValue()));
		return vertices;
	}

	public Collection<V> unpackVertices(short[] packedVertices) {
		Collection<V> vertices = new ArrayList<V>(packedVertices.length);
		for (short v : packedVertices)
			vertices.add(index.getKey((int) v));
		return vertices;
	}

	public Collection<Collection<V>> unpackAllVertices(Collection<short[]> packedVertices) {
		Collection<Collection<V>> c = new ArrayList<Collection<V>>(packedVertices.size());
		for (short[] v : packedVertices)
			c.add(unpackVertices(v));
		return c;
	}

	public int size() {
		return nextId;
	}

}
