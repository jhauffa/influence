package edu.tum.cs.graph.clique;

import java.io.File;
import java.util.logging.Logger;

import edu.tum.cs.graph.GraphCreator;
import edu.tum.cs.util.io.Serializer;
import edu.uci.ics.jung.graph.DirectedGraph;

public class FindSubgraphs {

	private static final Logger logger = Logger.getLogger(FindSubgraphs.class.getName());

	public static void main(String[] args) throws Exception {
		boolean requireComm = false;
		boolean requireBidir = false;
		boolean plexes = false;
		String pathPrefix = ".";
		int minSize = 3;
		for (String arg : args) {
			if (arg.equals("comm"))
				requireComm = true;
			else if (arg.equals("bidir"))
				requireBidir = true;
			else if (arg.startsWith("plex"))
				plexes = true;
			else if (arg.startsWith("m"))
				minSize = Integer.parseInt(arg.substring(1));
			else
				pathPrefix = arg;
		}
		String subgraphType = plexes ? "plexes" : "cliques";

		logger.info("Generating social network graph");
		DirectedGraph<Long, Integer> graph = GraphCreator.createGraph(requireComm);

		String baseName = subgraphType + "-" + graph.getVertexCount() + (requireComm ? "-comm" : "") +
				"-" + requireBidir;

		logger.info("Calculating " + subgraphType);
		long startTime = System.currentTimeMillis();

		VertexIndex<Long> index = new VertexIndex<Long>(graph);
		MaximalClique<Long, Integer> subgraphFinder;
		if (plexes)
			subgraphFinder = new MaximalPlex<Long, Integer>(graph, index, 1, requireBidir, minSize);
		else
			subgraphFinder = new MaximalClique<Long, Integer>(graph, index, requireBidir, minSize);
		subgraphFinder.saveAllMaximalSubgraphs(new File(pathPrefix, baseName + ".txt.gz"));

		long endTime = System.currentTimeMillis();
		logger.info("Finished, calculation took " + (endTime - startTime) + " milliseconds.");

		File f = new File(pathPrefix, baseName + "-node.ser.gz");
		logger.info("Writing per-node results to '" + f.getPath() + "'");
		Serializer.saveObjectToFile(subgraphFinder.getLargestMaximalSubgraphForNodes(), f);
		logger.info("Done.");
	}

}
