package edu.tum.cs.influence.meso.eval;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import edu.tum.cs.graph.clique.MaximalClique;
import edu.tum.cs.graph.clique.VertexIndex;
import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.DirectedSparseGraph;

public class HomogeneousSubsets {

	private static final String sep = ";";
	private static final double alpha = 0.01;

	public static void main(String[] args) throws IOException {
		if (args.length < 2) {
			System.err.println("usage: " + HomogeneousSubsets.class.getSimpleName() + " multcomp margmean\n" +
					"\tmultcomp\tresults of multiple comparisons in CSV format\n" +
					"\tmargmean\tmarginal means in CSV format\n");
			return;
		}

		// read marginal means (expected format: "indicatorWeight;mean")
		Map<String, Double> margMeans = new HashMap<String, Double>();
		BufferedReader reader = new BufferedReader(new FileReader(args[1]));
		try {
			String line = reader.readLine();	// read and discard header
			while ((line = reader.readLine()) != null) {
				String[] data = line.split(sep);
				double mean = Double.parseDouble(data[1]);
				margMeans.put(data[0], mean);
			}
		} finally {
			reader.close();
		}

		// read multiple comparison results (expected format: "indicatorWeight1;indicatorWeight2;*;*;p")
		DirectedGraph<String, Integer> subsets = new DirectedSparseGraph<String, Integer>();
		reader = new BufferedReader(new FileReader(args[0]));
		try {
			String line = reader.readLine();	// read and discard header
			int i = 0;
			while ((line = reader.readLine()) != null) {
				String[] data = line.split(sep);
				subsets.addVertex(data[0]);
				subsets.addVertex(data[1]);

				// "p-value for a test that the difference between the corresponding two marginal means is 0"
				// "The small p-values [...] indicate that the estimated marginal means [...] significantly differ from
				//  each other."
				double p = Double.parseDouble(data[4]);
				if (p > alpha) {
					subsets.addEdge(i++, data[0], data[1]);
					// MATLAB explicitly lists both directions, R has chosen a more compact representation
					subsets.addEdge(i++, data[1], data[0]);
				}
			}
		} finally {
			reader.close();
		}

		System.out.println("homogeneous groups, alpha = " + alpha);
		System.out.println("group;min;max;mean;members");
		VertexIndex<String> index = new VertexIndex<String>(subsets);
		MaximalClique<String, Integer> mc = new MaximalClique<String, Integer>(subsets, index, true, 1);
		int i = 1;
		for (Collection<String> cluster : index.unpackAllVertices(mc.getAllMaximalSubgraphs())) {
			SummaryStatistics stats = new SummaryStatistics();
			for (String s : cluster)
				stats.addValue(margMeans.get(s));

			System.out.print((i++) + ";" + stats.getMin() + ";" + stats.getMax() + ";" + stats.getMean());
			for (String s : cluster)
				System.out.print(";" + s + ";" + margMeans.get(s));
			System.out.println();
		}
		System.out.println();

		System.out.println("transitive closure, alpha = " + alpha);
		System.out.println("group;min;max;mean;members");
		WeakComponentClusterer<String, Integer> wcc = new WeakComponentClusterer<String, Integer>();
		i = 1;
		for (Collection<String> cluster : wcc.transform(subsets)) {
			SummaryStatistics stats = new SummaryStatistics();
			for (String s : cluster)
				stats.addValue(margMeans.get(s));

			System.out.print((i++) + ";" + stats.getMin() + ";" + stats.getMax() + ";" + stats.getMean());
			for (String s : cluster)
				System.out.print(";" + s + ";" + margMeans.get(s));
			System.out.println();
		}
	}

}
