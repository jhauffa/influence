package edu.tum.cs.time.hmm;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import com.google.common.primitives.Doubles;

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.graphics.Insets2D;
import de.erichseifert.gral.io.plots.DrawableWriter;
import de.erichseifert.gral.io.plots.DrawableWriterFactory;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import edu.tum.cs.math.dist.DiscreteDistribution;

public class TestGeneralization {

	private static int maxQuantiles = 10000;

	private static void addQuantile(DataTable qqData, double[] sortedData, double q, Percentile pct,
			RealDistribution dist, double[] limits) {
		double x = pct.evaluate(sortedData, q * 100.0);
		double y = dist.inverseCumulativeProbability(q);
		if (!Double.isInfinite(y)) {
			qqData.add(x, y);
			limits[0] = Math.min(limits[0], Math.min(x, y));
			limits[1] = Math.max(limits[1], Math.max(x, y));
		}
	}

	private static void writeQQPlot(OutputStream os, double[] sortedData, RealDistribution dist) throws Exception {
		DataTable qqData = new DataTable(2, Double.class);
		double[] limits = new double[] { 0.0, 1.0 };
		Percentile pct = (new Percentile()).withEstimationType(Percentile.EstimationType.R_7);	// mimic R's default
		if (sortedData.length <= maxQuantiles) {
			for (int i = 1; i <= sortedData.length; i++)
				addQuantile(qqData, sortedData, (double) i / sortedData.length, pct, dist, limits);
		} else {
			// Use the simple subsampling strategy described here (similar to Chebyshev nodes):
			// https://stats.stackexchange.com/questions/35220/removing-extraneous-points-near-the-centre-of-a-qq-plot
			// Spend half of the "quantile budget" on the center of the distribution, half on the tails.
			int numTail = maxQuantiles / 4;
			double centerShift = (double) numTail / sortedData.length;
			double centerScale = (double) (sortedData.length - (2 * numTail)) / sortedData.length;
			for (int i = 1; i <= maxQuantiles; i++) {
				double q;
				if (i < numTail) {
					q = (double) i / sortedData.length;
				} else if (i < (3 * numTail)) {
					int idx = i - numTail;
					q = (1.0 + Math.sin((double) idx / ((2 * numTail) + 1) * Math.PI - Math.PI / 2.0)) / 2.0;
					q = centerShift + (q * centerScale);
				} else {
					q = (double) (sortedData.length - maxQuantiles + i) / sortedData.length;
				}
				addQuantile(qqData, sortedData, q, pct, dist, limits);
			}
		}

		XYPlot plot = new XYPlot(qqData);
		DecimalFormat tickFmt = new DecimalFormat("0.##E0", DecimalFormatSymbols.getInstance(Locale.US));
		plot.getAxisRenderer(XYPlot.AXIS_X).setTickLabelFormat(tickFmt);
		plot.getAxisRenderer(XYPlot.AXIS_Y).setTickLabelFormat(tickFmt);
		plot.setInsets(new Insets2D.Double(10.0, 55.0, 35.0, 15.0));
		DataTable refData = new DataTable(2, Double.class);
		refData.add(limits[0], limits[0]);
		refData.add(limits[1], limits[1]);
		plot.add(refData);
		plot.setPointRenderers(refData, null);
		plot.setLineRenderers(refData, new DefaultLineRenderer2D());
		DrawableWriter dw = DrawableWriterFactory.getInstance().get("image/svg+xml");
		dw.write(plot, os, 1000, 1000);
	}

	private static void evaluateDist(String id, double[] sortedObs, AbstractRealDistribution dist, int numParam,
			String outputDir) throws Exception {
		double ll = 0.0;
		for (double v : sortedObs)
			ll += dist.logDensity(v);
		double bic = (-2 * ll) + (numParam * Math.log(sortedObs.length));
		double pp = DiscreteDistribution.perplexity(ll, sortedObs.length);

		System.out.println(id + ": LL = " + ll + ", BIC = " + bic + ", perplexity = " + pp);
		OutputStream os = new FileOutputStream(new File(outputDir, "qq-" + id + ".svg"));
		try {
			writeQQPlot(os, sortedObs, dist);
		} finally {
			os.close();
		}
	}

	private static double[] loadObservations(String fileName, SummaryStatistics deltaStats,
			SummaryStatistics logDeltaStats) throws Exception {
		ArrayList<Double> observations = new ArrayList<Double>();
		for (TimeSequence seq : TimeSequence.readSequences(new File(fileName))) {
			for (ObservationReal obs : seq.data) {
				if (obs.value > 0.0) {
					if (deltaStats != null)
						deltaStats.addValue(obs.value);
					if (logDeltaStats != null)
						logDeltaStats.addValue(Math.log(obs.value));
					observations.add(obs.value);
				}
			}
		}
		double[] sortedObs = Doubles.toArray(observations);
		Arrays.sort(sortedObs);
		return sortedObs;
	}

	private static final int maxObs = 1000000;

	public static void main(String[] args) throws Exception {
		if (args.length < 2) {
			System.err.println("usage: " + TestGeneralization.class.getSimpleName() +
					" chains-train chains-test [out]");
			return;
		}
		String outputDir = ".";
		if (args.length > 2)
			outputDir = args[2];

		SummaryStatistics deltaStats = new SummaryStatistics();
		SummaryStatistics logDeltaStats = new SummaryStatistics();
		double[] sortedObsTrain = loadObservations(args[0], deltaStats, logDeltaStats);
		double[] sortedObsTest = loadObservations(args[1], null, null);

		double[] sampleObs = null;
		if (sortedObsTrain.length > maxObs) {
			System.out.println("subsampling " + sortedObsTrain.length + " observations to " + maxObs +
					" for fitting of power-law distribution");
			sampleObs = new double[maxObs];
			int s = sortedObsTrain.length / maxObs;
			for (int i = 0; i < maxObs; i++)
				sampleObs[i] = sortedObsTrain[i * s];
		}

		NormalDistribution distNormal = new NormalDistribution(deltaStats.getMean(), deltaStats.getStandardDeviation());
		System.out.println("normal distribution: mean = " + distNormal.getMean() + ", sd = " +
				distNormal.getStandardDeviation());
		evaluateDist("normal-train", sortedObsTrain, distNormal, 2, outputDir);
		evaluateDist("normal-test", sortedObsTest, distNormal, 2, outputDir);

		LogNormalDistribution distLogNormal = new LogNormalDistribution(logDeltaStats.getMean(),
				logDeltaStats.getStandardDeviation());
		System.out.println("log-normal distribution: mean = " + distLogNormal.getScale() + ", sd = " +
				distLogNormal.getShape());
		evaluateDist("log-normal-train", sortedObsTrain, distLogNormal, 2, outputDir);
		evaluateDist("log-normal-test", sortedObsTest, distLogNormal, 2, outputDir);

		ExponentialDistribution distExp = new ExponentialDistribution(deltaStats.getMean());
		System.out.println("exponential distribution: mean = " + distExp.getMean());
		evaluateDist("exp-train", sortedObsTrain, distExp, 1, outputDir);
		evaluateDist("exp-test", sortedObsTest, distExp, 1, outputDir);
	}

}
