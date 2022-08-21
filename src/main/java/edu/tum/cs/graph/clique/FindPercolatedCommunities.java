package edu.tum.cs.graph.clique;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Writer;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;

import edu.tum.cs.math.dist.HistogramImpl.IdentityHistogram;
import edu.tum.cs.util.io.Serializer;

public class FindPercolatedCommunities {

	private static final Logger logger = Logger.getLogger(FindPercolatedCommunities.class.getName());
	private static int k = 3;
	private static boolean useLifoOrder = false;
	private static String pathPrefix = ".";

	private static int loadCliques(File cliqueFile, CliquePercolationParallel cf, VertexIndex<Long> index)
			throws IOException {
		IdentityHistogram histo = new IdentityHistogram();

		int n = 0;
		BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(
				new FileInputStream(cliqueFile))), 32 * 1024 * 1024);
		try {
			String line;
			while ((line = reader.readLine()) != null) {
				String[] part = line.split(";");
				histo.addValue(part.length);

				short[] clique = new short[part.length];
				for (int i = 0; i < part.length; i++) {
					Long node = Long.parseLong(part[i]);
					clique[i] = index.getOrInsert(node);
				}
				cf.addClique(clique);
				if ((++n % 10000000) == 0)
					logger.info("loaded " + n + " cliques");
			}
		} finally {
			reader.close();
		}

		System.out.println(histo);
		return n;
	}

	public static void saveCommunitiesAsText(Collection<Collection<Long>> communities, File f)
			throws IOException {
		Writer writer = new BufferedWriter(new FileWriter(f));
		try {
			for (Collection<Long> community : communities) {
				int idx = 0;
				for (long v : community) {
					if (idx++ > 0)
						writer.write(" ");
					writer.write(Long.toString(v));
				}
				writer.write("\n");
			}
		} finally {
			writer.close();
		}
	}

	public static void printCoverage(Collection<Collection<Long>> communities, int numInputNodes) {
		Set<Long> coveredNodes = new HashSet<Long>();
		for (Collection<Long> community : communities)
			for (Long node : community)
				coveredNodes.add(node);
		logger.info((((double) coveredNodes.size() / numInputNodes) * 100.0) + "% coverage (of " + numInputNodes +
				" nodes)");
	}

	public static void main(String[] args) throws IOException {
		if (args.length == 0) {
			logger.severe("usage: " + FindPercolatedCommunities.class.getSimpleName() + " cliques [clique size] " +
					"[lifo/pq] [output directory]");
			return;
		}
		String cliqueFileName = args[0];
		if (args.length > 1)
			k = Integer.parseInt(args[1]);
		if (args.length > 2)
			useLifoOrder = args[2].equalsIgnoreCase("lifo");
		if (args.length > 3)
			pathPrefix = args[3];

		CliquePercolationParallel cf = new CliquePercolationParallel(k, useLifoOrder);

		logger.info("loading cliques...");
		VertexIndex<Long> index = new VertexIndex<Long>();
		int numCliques = loadCliques(new File(cliqueFileName), cf, index);
		logger.info("loaded " + numCliques + " maximal cliques");

		logger.info("extracting communities for k = " + k);
		long startTime = System.currentTimeMillis();
		Collection<Collection<Long>> communities = index.unpackAllVertices(cf.getCommunities());
		long endTime = System.currentTimeMillis();
		logger.info("found " + communities.size() + " communities in " + (endTime - startTime) + " milliseconds");

		File f = new File(pathPrefix, "communities-" + k + ".ser");
		logger.info("writing communities to '" + f.getPath() + "'...");
		Serializer.saveObjectToFile(communities, f);
		saveCommunitiesAsText(communities, new File(pathPrefix, "onmi-" + k + ".txt"));

		// This specifically computes the coverage of nodes that are member of at least one clique!
		printCoverage(communities, index.size());
		logger.info("done");
	}

}
