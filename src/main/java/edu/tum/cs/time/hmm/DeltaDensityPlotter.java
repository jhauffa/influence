package edu.tum.cs.time.hmm;

import java.awt.Color;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import com.google.common.primitives.Doubles;

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.graphics.Insets2D;
import de.erichseifert.gral.io.plots.DrawableWriter;
import de.erichseifert.gral.io.plots.DrawableWriterFactory;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import de.erichseifert.gral.plots.lines.LineRenderer;

public class DeltaDensityPlotter {

	/* Kernel density estimator adapted from "Smile" by Haifeng Li. */
	private static class KernelDensity {
		private final double[] sortedData;
		private final double h;
		private final NormalDistribution gaussian;

		public KernelDensity(double[] sortedData, double sd) {
			this.sortedData = sortedData;
			int n = sortedData.length;
			double iqr = sortedData[n * 3 / 4] - sortedData[n / 4];
			h = 1.06 * Math.min(sd, iqr / 1.34) / Math.pow(n, 0.2);
			gaussian = new NormalDistribution(0.0, h);
		}

		public double density(double x) {
			int start = Arrays.binarySearch(sortedData, x - 5 * h);
			if (start < 0)
				start = -start - 1;
			int end = Arrays.binarySearch(sortedData, x + 5 * h);
			if (end < 0)
				end = -end - 1;

			double p = 0.0;
			for (int i = start; i < end; i++)
				p += gaussian.density(sortedData[i] - x);
			return p / sortedData.length;
		}
	}

	private static final int numPointsDensity = 1000;
	private static final double globalMinX = 1.0;
	private static final double globalMinY = 1E-11;

	private static double[][] buildDensityCurve(double[] sortedData, double sd, double[] limits) {
		KernelDensity dist = new KernelDensity(sortedData, sd);
		double[][] points = new double[numPointsDensity + 1][2];
		double minY = 1.0, maxY = 0.0;
		double maxExp = Math.log(sortedData[sortedData.length - 1] - sortedData[0] + 1.0);
		for (int i = 0; i <= numPointsDensity; i++) {
			double x = sortedData[0] + Math.exp(((double) i / numPointsDensity) * maxExp) - 1.0;
			double y = dist.density(x);
			points[i][0] = x;
			points[i][1] = y;
			if (y != 0.0) {
				if (y < minY)
					minY = y;
				if (y > maxY)
					maxY = y;
			}
		}
		minY = Math.pow(10.0, Math.floor(Math.log10(minY)));
		maxY = Math.pow(10.0, Math.ceil(Math.log10(maxY)));
		if (minY < limits[0])
			limits[0] = minY;
		if (maxY > limits[1])
			limits[1] = maxY;
		return points;
	}

	private static DataTable filterDensityCurve(double[][] curve, double[] limits) {
		DataTable points = new DataTable(2, Double.class);
		for (double[] point : curve) {
			if (point[1] != 0.0)
				points.add(point[0], point[1]);
			else
				points.add(point[0], limits[0]);
		}
		return points;
	}

	private static XYPlot buildDensityPlot(DataTable[] curves, Color[] colors, double[] limits) {
		XYPlot plot = new XYPlot();
		plot.setInsets(new Insets2D.Double(10.0, 70.0, 30.0, 10.0));
		for (int i = 0; i < curves.length; i++) {
			plot.add(curves[i]);
			plot.setPointRenderers(curves[i], null);
			LineRenderer lr = new DefaultLineRenderer2D();
			lr.setColor(colors[i]);
			plot.setLineRenderers(curves[i], lr);
		}
		if (plot.getAxis(XYPlot.AXIS_X).getMin().doubleValue() < globalMinX)
			plot.getAxis(XYPlot.AXIS_X).setMin(globalMinX);
		plot.getAxis(XYPlot.AXIS_Y).setRange(limits[0], limits[1] + (limits[1] / 100));
		plot.setAxisRenderer(XYPlot.AXIS_X, new LogTimeAxisRenderer());
		plot.setAxisRenderer(XYPlot.AXIS_Y, new LogProbabilityAxisRenderer(limits[0]));
		return plot;
	}

	private static double[] loadObservations(String fileName, SummaryStatistics stats) throws Exception {
		ArrayList<Double> observations = new ArrayList<Double>();
		for (TimeSequence seq : TimeSequence.readSequences(new File(fileName))) {
			for (ObservationReal obs : seq.data) {
				if (obs.value > 0.0) {
					observations.add(obs.value);
					stats.addValue(obs.value);
				}
			}
		}
		double[] sortedObs = Doubles.toArray(observations);
		observations = null;
		Arrays.sort(sortedObs);
		return sortedObs;
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 2) {
			System.err.println("usage: " + DeltaDensityPlotter.class.getSimpleName() +
					" chains color [chains color ...]");
			return;
		}

		int numCurves = args.length / 2;
		double[][][] densityCurveData = new double[numCurves][][];
		Color[] colors = new Color[numCurves];
		double[] limits = { 1.0, 0.0 };
		for (int i = 0; i < numCurves; i++) {
			SummaryStatistics stats = new SummaryStatistics();
			double[] sortedObs = loadObservations(args[i * 2], stats);
			colors[i] = Color.decode(args[(i * 2) + 1]);
			densityCurveData[i] = buildDensityCurve(sortedObs, stats.getStandardDeviation(), limits);
		}
		limits[0] = Math.max(limits[0], globalMinY);
		DataTable[] densityCurves = new DataTable[numCurves];
		for (int i = 0; i < numCurves; i++)
			densityCurves[i] = filterDensityCurve(densityCurveData[i], limits);
		densityCurveData = null;

		// generate density plot
		XYPlot densityPlot = buildDensityPlot(densityCurves, colors, limits);
		DrawableWriter dw = DrawableWriterFactory.getInstance().get("image/svg+xml");
		OutputStream os = new FileOutputStream(new File("density-est.svg"));
		try {
			dw.write(densityPlot, os, 600, 400);
		} finally {
			os.close();
		}
	}

}
