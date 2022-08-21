package edu.tum.cs.nlp.topic;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ml.clustering.CentroidCluster;
import org.apache.commons.math3.ml.clustering.Clusterable;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;
import org.apache.commons.math3.ml.distance.EuclideanDistance;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.HmmBuilder;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import edu.tum.cs.math.MultiDimensionalScaling;
import edu.tum.cs.math.dist.DiscreteDistributionSampler;
import edu.tum.cs.nlp.corpus.Index;
import edu.tum.cs.nlp.topic.FitSocialTopicModel.SocialTopicModelParameters;
import edu.tum.cs.time.hmm.HmmDivergence;
import edu.tum.cs.time.hmm.HmmPlotter;
import edu.tum.cs.time.hmm.InteractionHmm;
import edu.tum.cs.util.arrays.Sparse2DIterator;
import edu.tum.cs.util.io.Serializer;

public class AnalyzeSocialTopicModel {

	private static class HmmParameters implements Clusterable, Serializable {
		private static final long serialVersionUID = -4630385431288238478L;

		public final int idx;
		public final double[][] pTrans, pOut;
		public final int seqLength;
		private double[] embedding;
		public List<Integer> clusterAssignments = new ArrayList<Integer>();

		public HmmParameters(int idx, double[][] pTrans, double[][] pOut, int seqLength) {
			this.idx = idx;
			this.pTrans = pTrans;
			this.pOut = pOut;
			this.seqLength = seqLength;
		}

		public void setPoint(double[] embedding) {
			this.embedding = embedding;
		}

		@Override
		public double[] getPoint() {
			return embedding;
		}

		public boolean equals(HmmParameters other, boolean transitionMatrixOnly) {
			boolean isEqual = Arrays.deepEquals(pTrans, other.pTrans);
			if (isEqual && !transitionMatrixOnly)
				isEqual = Arrays.deepEquals(pOut, other.pOut);
			return isEqual;
		}
	}

	private static final int maxN = 1000;
	private static final int[] numClusters = { 2, 3, 4, 5, 6 };

	private static void printMatrix(double[][] m) {
		for (double[] p : m) {
			System.out.print("\t");
			for (double v : p)
				System.out.print(v + ";");
			System.out.println();
		}
	}

	private static void writeTransitionPlot(String id, double[][] pTrans) throws Exception {
		Hmm<ObservationReal> hmm = (new HmmBuilder<ObservationReal>(pTrans.length)).withA(pTrans).build();
		InteractionHmm ihmm = new InteractionHmm(hmm, false);
		FileOutputStream os = new FileOutputStream(id + ".svg");
		try {
			HmmPlotter.plotTransitions(os, ihmm, 0.01, true, true);
		} finally {
			os.close();
		}
	}

	private static void embedHmmParameters(ArrayList<HmmParameters> hmms, boolean transitionDistOnly) {
		// construct the dissimilarity matrix
		int n = hmms.size();
		double[][] dist = new double[n][n];
		for (int i = 0; i < n; i++) {
			HmmParameters paramI = hmms.get(i);
			for (int j = i + 1; j < n; j++) {
				HmmParameters paramJ = hmms.get(j);
				double dij;
				if (transitionDistOnly)
					dij = HmmDivergence.distTransitionMatrix(paramI.pTrans, paramJ.pTrans);
				else
					dij = HmmDivergence.divergenceRate(paramI.pTrans, paramI.pOut, paramJ.pTrans, paramJ.pOut);
				if ((dij == 0.0) && (i != j))
					throw new RuntimeException("i = " + i + ", j = " + j);
				dist[i][j] = dist[j][i] = dij;
			}
		}

		// embed each set of parameters in a 2D Euclidean space
		MultiDimensionalScaling mds = new MultiDimensionalScaling(dist);
		RealMatrix emb = mds.getEmbedding(2);
		for (int i = 0; i < n; i++)
			hmms.get(i).setPoint(new double[] { emb.getEntry(i, 0), emb.getEntry(i, 1) });
	}

	private static void analyzeSocialTopicModel(ArrayList<HmmParameters> hmms, boolean embedOnly,
			boolean transitionDistOnly, int minSeqLength) throws Exception {
		System.out.println("=== embedding" + (embedOnly ? "" : "+clustering") + ", " +
			(transitionDistOnly ? "JSD of transition matrix" : "divergence rate of Yang et al.") +
			", min. sequence length = " + minSeqLength);

		if (minSeqLength > 0) {
			ArrayList<HmmParameters> filteredHmms = new ArrayList<HmmParameters>();
			for (HmmParameters hmm : hmms) {
				if (hmm.seqLength >= minSeqLength)
					filteredHmms.add(hmm);
			}
			hmms = filteredHmms;
		}

		ArrayList<HmmParameters> uniqueHmms = new ArrayList<HmmParameters>();
		for (HmmParameters hmm : hmms) {
			boolean isDuplicate = false;
			for (HmmParameters otherHmm : uniqueHmms) {
				isDuplicate = otherHmm.equals(hmm, transitionDistOnly);
				if (isDuplicate)
					break;
			}
			if (!isDuplicate)
				uniqueHmms.add(hmm);
		}
		System.out.println("before deduplication: " + hmms.size() + ", after: " + uniqueHmms.size());
		hmms = uniqueHmms;

		if (hmms.size() > maxN) {
			ArrayList<HmmParameters> subSample = new ArrayList<HmmParameters>(maxN);
			DiscreteDistributionSampler sampler = new DiscreteDistributionSampler();
			for (int idx : sampler.sample1DUniformWithoutReplacement(hmms.size(), maxN))
				subSample.add(hmms.get(idx));
			hmms = subSample;
		}

		embedHmmParameters(hmms, transitionDistOnly);

		// perform clustering (k-means with different numbers of clusters)
		for (HmmParameters hmm : hmms)
			hmm.clusterAssignments.clear();
		List<String> clusteringOrder = new ArrayList<String>();
		if (!embedOnly) {
			EuclideanDistance euclideanDist = new EuclideanDistance();
			for (int c : numClusters) {
				KMeansPlusPlusClusterer<HmmParameters> clustKm = new KMeansPlusPlusClusterer<HmmParameters>(c);
				List<CentroidCluster<HmmParameters>> clustersKm = clustKm.cluster(hmms);
				clusteringOrder.add("km" + c);
				System.out.println("k-means with " + c + " clusters");
				int clusterIdx = 0;
				for (CentroidCluster<HmmParameters> cluster : clustersKm) {
					HmmParameters exemplar = null;
					double minDist = Double.MAX_VALUE;
					for (HmmParameters param : cluster.getPoints()) {
						double dist = euclideanDist.compute(param.getPoint(), cluster.getCenter().getPoint());
						if (dist < minDist) {
							exemplar = param;
							minDist = dist;
						}
						param.clusterAssignments.add(clusterIdx);
					}

					System.out.println("\texemplar for cluster " + clusterIdx + " (index " + exemplar.idx + ")");
					System.out.println("\tpi:");
					printMatrix(exemplar.pTrans);
					System.out.println("\ttheta:");
					printMatrix(exemplar.pOut);
					writeTransitionPlot("km-" + (transitionDistOnly ? "trans-" : "yang-") + c + "-cluster" + clusterIdx,
							exemplar.pTrans);
					clusterIdx++;
				}
			}
			System.out.println();
		}

		// print embedding coordinates and cluster assignments in CSV format
		System.out.print("idx;x;y;seqlen");
		for (String clusteringId : clusteringOrder)
			System.out.print(";" + clusteringId);
		System.out.println();
		for (HmmParameters param : hmms) {
			System.out.print(param.idx + ";" + param.embedding[0] + ";" + param.embedding[1] + ";" + param.seqLength);
			for (Integer c : param.clusterAssignments)
				System.out.print(";" + c);
			System.out.println();
		}
		System.out.println();
	}

	private static ArrayList<HmmParameters> filterHmmParameters(SocialTopicModelParameters soc,
			boolean excludeNonAddressive) {
		// collect parameters pi and theta for each author-recipient pair
		ArrayList<HmmParameters> hmms = new ArrayList<HmmParameters>();
		int idx = 0;
		Sparse2DIterator<double[][]> it = soc.pi.sparseIterator();
		while (it.hasNext()) {
			it.advance();
			int author = it.getRow();
			int recipient = it.getColumn();
			if (excludeNonAddressive && (author == recipient))
				continue;

			double[][] pi = it.getValue();
			double[][] theta = soc.theta.get(author, recipient);
			int seqLength = soc.seqLengths.get((author << 16) | recipient);
			hmms.add(new HmmParameters(idx++, pi, theta, seqLength));
		}
		return hmms;
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 2) {
			System.err.println("usage: " + AnalyzeSocialTopicModel.class.getSimpleName() + " parameters bag-of-words " +
					"[min. seq. length]");
			return;
		}

		SocialTopicModelParameters soc = Serializer.loadObjectFromFile(new File(args[0]));
		Index<String> bagOfWords = Serializer.loadObjectFromFile(new File(args[1]));
		TopicWordDistribution.writeTopicsCsv(soc.phi, bagOfWords, new PrintWriter(System.out));
		System.out.println();

		int minSeqLength = 2;
		if (args.length >= 3)
			minSeqLength = Integer.parseInt(args[2]);

		ArrayList<HmmParameters> hmms = filterHmmParameters(soc, true);
		analyzeSocialTopicModel(hmms, true, false, 0);
		analyzeSocialTopicModel(hmms, false, false, minSeqLength);
		analyzeSocialTopicModel(hmms, true, true, 0);
		analyzeSocialTopicModel(hmms, false, true, minSeqLength);
	}

}
