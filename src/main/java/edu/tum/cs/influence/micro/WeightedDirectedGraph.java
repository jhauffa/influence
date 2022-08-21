package edu.tum.cs.influence.micro;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.io.GraphMLMetadata;
import edu.uci.ics.jung.io.GraphMLWriter;
import edu.uci.ics.jung.io.graphml.EdgeMetadata;
import edu.uci.ics.jung.io.graphml.GraphMLReader2;
import edu.uci.ics.jung.io.graphml.GraphMetadata;
import edu.uci.ics.jung.io.graphml.HyperEdgeMetadata;
import edu.uci.ics.jung.io.graphml.NodeMetadata;

public class WeightedDirectedGraph<V, W> extends DirectedSparseGraph<V, WeightedEdge<W>> {

	private static final long serialVersionUID = 8443516872400472416L;

	public void writeGraphMl(File f) throws IOException {
		GraphMLWriter<V, WeightedEdge<W>> writer = new GraphMLWriter<V, WeightedEdge<W>>();
		Map<String, GraphMLMetadata<WeightedEdge<W>>> edgeMetadata =
				new HashMap<String, GraphMLMetadata<WeightedEdge<W>>>();
		edgeMetadata.put("weight", new GraphMLMetadata<WeightedEdge<W>>("", "",
				new Transformer<WeightedEdge<W>, String>() {
					@Override
					public String transform(WeightedEdge<W> e) {
						if (e.weight.getClass().isArray()) {
							// Convert array to string in a way that works for arrays of Object as well as for arrays of
							// primitive data types.
							int length = Array.getLength(e.weight);
							Object[] objArr = new Object[length];
							for (int i = 0; i < length; i++)
								objArr[i] = Array.get(e.weight, i);
							return Arrays.toString(objArr);
						}
						return e.weight.toString();
					}
		}));
		writer.setEdgeData(edgeMetadata);
		writer.save(this, new FileWriter(f));
	}

	public static <V, W> WeightedDirectedGraph<V, W> readGraphMl(File f,
			final Transformer<String, V> vertexTransformer,
			final Transformer<String, W> weightTransformer) throws Exception {
		Transformer<GraphMetadata, WeightedDirectedGraph<V, W>> graphTransformer =
				new Transformer<GraphMetadata, WeightedDirectedGraph<V, W>>() {
					@Override
					public WeightedDirectedGraph<V, W> transform(GraphMetadata m) {
						return new WeightedDirectedGraph<V, W>();
					}
		};
		Transformer<NodeMetadata, V> wrappedVertexTransformer = new Transformer<NodeMetadata, V>() {
					@Override
					public V transform(NodeMetadata m) {
						return vertexTransformer.transform(m.getId());
					}
		};
		Transformer<EdgeMetadata, WeightedEdge<W>> edgeTransformer =
				new Transformer<EdgeMetadata, WeightedEdge<W>>() {
					private int idx = 0;

					@Override
					public WeightedEdge<W> transform(EdgeMetadata m) {
						return new WeightedEdge<W>(idx++, weightTransformer.transform(m.getProperty("weight")));
					}
		};
		Transformer<HyperEdgeMetadata, WeightedEdge<W>> hyperEdgeTransformer =
				new Transformer<HyperEdgeMetadata, WeightedEdge<W>>() {
					@Override
					public WeightedEdge<W> transform(HyperEdgeMetadata m) {
						throw new RuntimeException("weighted graph cannot have hyper edges");
					}
		};

		GraphMLReader2<WeightedDirectedGraph<V, W>, V, WeightedEdge<W>> reader =
				new GraphMLReader2<WeightedDirectedGraph<V, W>, V, WeightedEdge<W>>(new FileReader(f),
						graphTransformer, wrappedVertexTransformer, edgeTransformer, hyperEdgeTransformer);
		try {
			reader.init();
			return reader.readGraph();
		} finally {
			reader.close();
		}
	}

	public static WeightedDirectedGraph<Long, Double> readInfluenceGraph(File f) throws Exception {
		return WeightedDirectedGraph.<Long, Double>readGraphMl(f,
				new Transformer<String, Long>() {
					@Override
					public Long transform(String id) {
						return Long.parseLong(id);
					}
				},
				new Transformer<String, Double>() {
					@Override
					public Double transform(String weightData) {
						return Double.valueOf(weightData);
					}
				});
	}

}
