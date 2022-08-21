package edu.tum.cs.influence.micro;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

import edu.tum.cs.math.dist.HistogramImpl.QuantizedProbabilityHistogram;
import edu.tum.cs.util.ExperimentConfiguration;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

public class EvaluateNetworks {

	public static class EvaluationResult {
		public final TimeSeriesVariant dataVariant;
		public final ModelVariant modelVariant;
		public final String graphFileName;
		public WeightedDirectedGraph<Long, Integer> annotPositive;
		public double precision, specificity, recallLowerBound, botMisclassificationRate;

		public EvaluationResult(TimeSeriesVariant dataVariant, ModelVariant modelVariant, String graphFileName) {
			this.dataVariant = dataVariant;
			this.modelVariant = modelVariant;
			this.graphFileName = graphFileName;
		}

		@Override
		public String toString() {
			return dataVariant + ";" + modelVariant + ";" + graphFileName + ";" + precision + ";" + specificity + ";" +
					recallLowerBound + ";" + botMisclassificationRate;
		}

		private static String getCsvHeader() {
			return TimeSeriesVariant.getCsvHeader() + ";" + ModelVariant.getCsvHeader() +
					";graph;precision;specificity;recallLowerBound;botMisclassificationRate";
		}

		public static ArrayList<EvaluationResult> readResults(String fileName) throws IOException {
			ArrayList<EvaluationResult> results = new ArrayList<EvaluationResult>();
			BufferedReader r = new BufferedReader(new FileReader(fileName));
			try {
				String line = r.readLine();	// skip column header
				while ((line = r.readLine()) != null) {
					String[] parts = line.split(";");
					Queue<String> remainingParts = new LinkedList<String>(Arrays.asList(parts));
					TimeSeriesVariant dataVariant = TimeSeriesVariant.readCsv(remainingParts);
					ModelVariant modelVariant = ModelVariant.readCsv(remainingParts);
					String graphFileName = remainingParts.poll();
					results.add(new EvaluationResult(dataVariant, modelVariant, graphFileName));
				}
			} finally {
				r.close();
			}
			return results;
		}

		public static void writeResults(String fileName, ArrayList<EvaluationResult> results) throws IOException {
			PrintWriter w = new PrintWriter(new FileWriter(fileName));
			try {
				w.println(getCsvHeader());
				for (EvaluationResult r : results)
					w.println(r);
			} finally {
				w.close();
			}
		}
	}

	private static <T> Map<T, WeightedDirectedGraph<Long, Integer>> readAnnotations(String fileName, List<T> categories)
			throws IOException {
		Map<String, T> colMap = new HashMap<String, T>();
		for (T c : categories)
			colMap.put(c.toString(), c);

		Map<T, WeightedDirectedGraph<Long, Integer>> annot = new HashMap<T, WeightedDirectedGraph<Long, Integer>>();
		int edgeIdx = 0;
		BufferedReader r = new BufferedReader(new FileReader(fileName));
		try {
			String line = r.readLine();	// first line is column header
			String[] header = line.split(";");
			for (int i = 0; i < header.length - 5; i++)
				annot.put(colMap.get(header[4 + i]), new WeightedDirectedGraph<Long, Integer>());

			while ((line = r.readLine()) != null) {
				String[] parts = line.split(";");
				Long v1 = Long.parseLong(parts[0]);
				Long v2 = Long.parseLong(parts[1]);
				Integer annotation = Integer.parseInt(parts[parts.length - 1]);
				for (int i = 0; i < parts.length - 5; i++) {
					if (Integer.parseInt(parts[4 + i]) != 0) {
						annot.get(colMap.get(header[4 + i])).addEdge(new WeightedEdge<Integer>(edgeIdx++, annotation),
								v1, v2, EdgeType.DIRECTED);
					}
				}
			}
		} finally {
			r.close();
		}
		return annot;
	}

	private static void readPositiveAnnotations(ArrayList<EvaluationResult> experiments, String fileName)
			throws IOException {
		List<Integer> columns = new ArrayList<Integer>(experiments.size());
		for (int i = 0; i < experiments.size(); i++)
			columns.add(i);
		Map<Integer, WeightedDirectedGraph<Long, Integer>> annot = readAnnotations(fileName, columns);
		for (Map.Entry<Integer, WeightedDirectedGraph<Long, Integer>> e : annot.entrySet())
			experiments.get(e.getKey()).annotPositive = e.getValue();
	}

	private static final int numBins = 25;

	private static void computeEvaluationMetrics(EvaluationResult r, int variantIdx, int maxIdx,
			DirectedSparseGraph<Long, boolean[]> unionGraph, QuantizedProbabilityHistogram[] histograms,
			WeightedDirectedGraph<Long, Integer> annotNegative, Set<Long> botIds, File resultsDir) throws Exception {
		if (r.annotPositive != null) {
			// check consistency of annotations: sets of positive and negative edges must be disjoint, positive set must
			// not contain any edges where the target node is on the list of bot accounts
			for (WeightedEdge<Integer> e : r.annotPositive.getEdges()) {
				Pair<Long> v = r.annotPositive.getEndpoints(e);
				if (annotNegative.findEdge(v.getFirst(), v.getSecond()) != null) {
					throw new RuntimeException("annotation mismatch for experiment " + r.dataVariant + ";" +
							r.modelVariant + ": edge " + v.getFirst() + "-" + v.getSecond() +
							" has both positive and negative annotations");
				}
				if (botIds.contains(v.getSecond())) {
					throw new RuntimeException("annotation mismatch for experiment " + r.dataVariant + ";" +
							r.modelVariant + ": bot account " + v.getSecond() +
							" has positive annotation for influence from" + v.getFirst());
				}
			}
		}

		// load influence graph, then update edge presence information and influence magnitude histogram
		WeightedDirectedGraph<Long, Double> graph = WeightedDirectedGraph.readInfluenceGraph(
				new File(resultsDir, r.graphFileName));
		histograms[variantIdx] = new QuantizedProbabilityHistogram(numBins);
		for (WeightedEdge<Double> e : graph.getEdges()) {
			Pair<Long> v = graph.getEndpoints(e);
			boolean[] isPresent = unionGraph.findEdge(v.getFirst(), v.getSecond());
			if (isPresent == null) {
				isPresent = new boolean[maxIdx];
				unionGraph.addEdge(isPresent, v, EdgeType.DIRECTED);
			}
			isPresent[variantIdx] = true;
			histograms[variantIdx].addValue(e.weight);
		}

		// compute precision, specificity, a lower bound on recall, and the bot misclassification rate
		int numTruePositive = 0;
		int numEdgeToBot = 0;
		for (WeightedEdge<Double> e : graph.getEdges()) {
			Pair<Long> v = graph.getEndpoints(e);
			if ((r.annotPositive != null) && (r.annotPositive.findEdge(v.getFirst(), v.getSecond()) != null))
				numTruePositive++;
			if (botIds.contains(v.getSecond()))
				numEdgeToBot++;
		}
		r.precision = (double) numTruePositive / graph.getEdgeCount();
		r.botMisclassificationRate = (double) numEdgeToBot / graph.getEdgeCount();

		int numTrueNegative = 0;
		for (WeightedEdge<Integer> e : annotNegative.getEdges()) {
			Pair<Long> v = annotNegative.getEndpoints(e);
			if (graph.findEdge(v.getFirst(), v.getSecond()) == null)
				numTrueNegative++;
		}
		r.specificity = (double) numTrueNegative / annotNegative.getEdgeCount();
		r.recallLowerBound = (1.0 - r.specificity) * (annotNegative.getEdgeCount() - 1) / ((1.0 / r.precision) - 1.0);
	}

	private static double computeJaccardIndex(int[][] numCommonEdges, int i, int j) {
		double denom = numCommonEdges[i][i] + numCommonEdges[j][j] - numCommonEdges[i][j];
		double jd = 1.0;
		if (denom != 0.0)
			jd = (double) numCommonEdges[i][j] / denom;
		return jd;
	}

	private static double computeOverlapCoeff(int[][] numCommonEdges, int i, int j) {
		double denom = Math.min(numCommonEdges[i][i], numCommonEdges[j][j]);
		double oc = 1.0;
		if (denom != 0.0)
			oc = (double) numCommonEdges[i][j] / denom;
		return oc;
	}

	private static void writeOverlap(String fileName, DirectedSparseGraph<Long, boolean[]> unionGraph,
			String[] variantIds, boolean useJaccard) throws IOException {
		int maxIdx = variantIds.length;
		int[][] numCommonEdges = new int[maxIdx][maxIdx];
		for (boolean[] isPresent : unionGraph.getEdges()) {
			for (int i = 0; i < isPresent.length; i++) {
				for (int j = 0; j <= i; j++) {
					if (isPresent[i] && isPresent[j]) {
						numCommonEdges[i][j]++;
						if (i != j)
							numCommonEdges[j][i]++;
					}
				}
			}
		}

		// compute Jaccard distance / overlap coefficient and write to CSV file
		PrintWriter w = new PrintWriter(new FileWriter(fileName));
		try {
			for (String id : variantIds)
				w.print(";" + id);
			w.println();
			for (int i = 0; i < maxIdx; i++) {
				w.print(variantIds[i]);
				for (int j = 0; j < maxIdx; j++) {
					double v;
					if (useJaccard)
						v = computeJaccardIndex(numCommonEdges, i, j);
					else
						v = computeOverlapCoeff(numCommonEdges, i, j);
					w.print(";" + v);
				}
				w.println();
			}
		} finally {
			w.close();
		}
	}

	private static void writeHistograms(String fileName, QuantizedProbabilityHistogram[] histograms,
			String[] variantIds) throws IOException {
		PrintWriter w = new PrintWriter(new FileWriter(fileName));
		try {
			for (int i = 0; i < variantIds.length; i++) {
				w.println(variantIds[i] + ";;;");
				w.println(histograms[i].toString());
				w.println(";;;");
			}
		} finally {
			w.close();
		}
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 2) {
			System.err.println("usage: " + EvaluateNetworks.class.getSimpleName() +
					" results negative [positive [bots]]");
		}

		// load results of individual influence measurement experiments
		ArrayList<EvaluationResult> results = EvaluationResult.readResults(args[0]);
		File resultsDir = (new File(args[0])).getAbsoluteFile().getParentFile();

		// load positive/negative/bot network annotations
		Map<ModelVariant.NetworkType, WeightedDirectedGraph<Long, Integer>> annotNegative = readAnnotations(args[1],
				Arrays.asList(ModelVariant.NetworkType.values()));
		Set<Long> botIds = Collections.<Long>emptySet();
		if (args.length > 2) {
			readPositiveAnnotations(results, args[2]);
			if (args.length > 3)
				botIds = ExperimentConfiguration.loadUserIds(args[3], Integer.MAX_VALUE);
		}

		// process GraphML file of each experiment and write results to CSV file
		int idx = 0, maxIdx = results.size();
		String[] variantIds = new String[maxIdx];
		DirectedSparseGraph<Long, boolean[]> unionGraph = new DirectedSparseGraph<Long, boolean[]>();
		QuantizedProbabilityHistogram[] histograms = new QuantizedProbabilityHistogram[maxIdx];
		for (EvaluationResult r : results) {
			System.err.println("processing " + r.graphFileName + " (" + (idx + 1) + "/" + results.size() + ")");
			variantIds[idx] = r.dataVariant.toString() + "-" + r.modelVariant.getFilePrefix();
			computeEvaluationMetrics(r, idx, maxIdx, unionGraph, histograms,
					annotNegative.get(r.modelVariant.getNetwork()), botIds, resultsDir);
			idx++;
		}
		EvaluationResult.writeResults("eval-post.csv", results);
		writeOverlap("eval-jaccard.csv", unionGraph, variantIds, true);
		writeOverlap("eval-overlap.csv", unionGraph, variantIds, false);
		writeHistograms("eval-histo.csv", histograms, variantIds);
		System.err.println("done");
	}

}
