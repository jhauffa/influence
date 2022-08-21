package edu.tum.cs.nlp.topic;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Locale;

import edu.tum.cs.nlp.corpus.Index;

public class TopicWordDistribution {

	private static final int numTopWordsCsv = 50;
	private static final int numTopWordsLatex = 5;

	public static int[] topWords(double[] topicWordDistr, int n) {
		int[] topWordIdx = new int[n];
		Arrays.fill(topWordIdx, -1);
		for (int i = 1; i < topicWordDistr.length; i++) {
			if ((topWordIdx[0] == -1) || (topicWordDistr[i] > topicWordDistr[topWordIdx[0]])) {
				topWordIdx[0] = i;
				int j = 1;
				while ((j < n) &&
					   ((topWordIdx[j] == -1) ||
						(topicWordDistr[topWordIdx[j - 1]] > topicWordDistr[topWordIdx[j]]))) {
					int tmp = topWordIdx[j];
					topWordIdx[j] = topWordIdx[j - 1];
					topWordIdx[j - 1] = tmp;
					j++;
				}
			}
		}
		return topWordIdx;
	}

	public static void writeTopicsCsv(double[][] topicWordDistr, Index<String> bagOfWords, PrintWriter writer)
			throws IOException {
		for (int i = 0; i < topicWordDistr.length; i++) {
			writer.println(i + 1);

			int[] topWords = topWords(topicWordDistr[i], numTopWordsCsv);
			for (int j = numTopWordsCsv - 1; j >= 0; j--)
				writer.println(bagOfWords.getElement(topWords[j]) + "\t" + topicWordDistr[i][topWords[j]]);
			writer.println();
		}
		writer.flush();
	}

	public static void saveTopicsCsv(double[][] topicWordDistr, Index<String> bagOfWords, File f) throws IOException {
		PrintWriter writer = new PrintWriter(new FileWriter(f));
		try {
			writeTopicsCsv(topicWordDistr, bagOfWords, writer);
		} finally {
			writer.close();
		}
	}

	public static void writeTopicsLatex(double[][] topicWordDistr, Index<String> bagOfWords, PrintWriter writer)
			throws IOException {
		DecimalFormat fmt = (DecimalFormat) NumberFormat.getInstance(Locale.US);
		fmt.applyPattern("0.0000");

		for (int i = 0; i < topicWordDistr.length; i++) {
			writer.println("\n\\centering{\\textbf{Topic " + (i + 1) + "}}\\\\");

			int[] topWords = topWords(topicWordDistr[i], numTopWordsLatex);
			for (int j = numTopWordsLatex - 1; j >= 0; j--) {
				String word = bagOfWords.getElement(topWords[j]);
				word = word.replaceAll("_", "\\_").replaceAll("%", "\\%");
				writer.println(word + "\\hfill" + fmt.format(topicWordDistr[i][topWords[j]]) + "\\\\");
			}

			if (((i + 1) % 5) == 0)
				writer.println("\\vfill\n\\columnbreak");
			writer.println();
		}
		writer.flush();
	}

	public static void saveTopicsLatex(double[][] topicWordDistr, Index<String> bagOfWords, File f) throws IOException {
		PrintWriter writer = new PrintWriter(new FileWriter(f));
		try {
			writeTopicsLatex(topicWordDistr, bagOfWords, writer);
		} finally {
			writer.close();
		}
	}

}
