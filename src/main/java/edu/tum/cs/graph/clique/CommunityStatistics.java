package edu.tum.cs.graph.clique;

import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import edu.tum.cs.math.dist.HistogramImpl.IdentityHistogram;
import edu.tum.cs.util.io.Serializer;

public class CommunityStatistics {

	public static void main(String[] args) throws Exception {
		Collection<Collection<Long>> communities = Serializer.loadObjectFromFile(new File(args[0]));
		int numNodes = Integer.parseInt(args[1]);

		Set<Long> coveredNodes = new HashSet<Long>();
		IdentityHistogram histo = new IdentityHistogram();
		DescriptiveStatistics stats = new DescriptiveStatistics();
		for (Collection<Long> c : communities) {
			int n = c.size();
			histo.addValue(n);
			stats.addValue(n);
			for (Long v : c)
				coveredNodes.add(v);
		}

		System.out.println(communities.size() + " communities, " + (((double) coveredNodes.size() / numNodes) * 100.0) +
				"% coverage (" + coveredNodes.size() + " of " + numNodes + " nodes)");
		System.out.println(stats);
		System.out.println("largest community includes " + ((stats.getMax() / coveredNodes.size()) * 100.0) +
				"% of covered nodes");
		System.out.println();
		System.out.println(histo);
	}

}
