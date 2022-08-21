package edu.tum.cs.time.hmm;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Locale;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import de.erichseifert.gral.data.DataSeries;
import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.graphics.Insets2D;
import de.erichseifert.gral.io.plots.DrawableWriter;
import de.erichseifert.gral.io.plots.DrawableWriterFactory;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.plots.axes.LinearRenderer2D;
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import be.ac.ulg.montefiore.run.jahmm.OpdfLogNormal;
import be.ac.ulg.montefiore.run.jahmm.io.HmmWriter;
import be.ac.ulg.montefiore.run.jahmm.io.OpdfGaussianWriter;
import be.ac.ulg.montefiore.run.jahmm.io.OpdfLogNormalWriter;

public class HmmPlotter {

	private static final String svgHeader =
			"<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n" +
			"<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" viewBox=\"0 0 1200 300\" height=\"63mm\" width=\"210mm\">\n" +
			"  <defs>\n" +
			"    <marker style=\"overflow:visible\" id=\"arrow\" refX=\"0.0\" refY=\"0.0\" orient=\"auto\" markerUnits=\"userSpaceOnUse\">\n" +
			"      <path style=\"stroke-width:0.625;stroke:#000000;fill:#000000\" transform=\"scale(1.1) rotate(180) translate(1,0)\"\n" +
			"        d=\"M 8.7185878,4.0337352 L -2.2072895,0.016013256 L 8.7185884,-4.0017078 C 6.9730900,-1.6296469 6.9831476,1.6157441 8.7185878,4.0337352 z\" />\n" +
			"    </marker>\n" +
			"  </defs>\n";
	private static final String svgFooter =
			"</svg>\n";

	private static final long[] durations = new long[] { 365 * 24 * 60 * 60, 24 * 60 * 60, 60 * 60, 60, 1 };
	private static final String[] suffixes = new String[] { "y", "d", "h", "m", "s" };

	public static String formatDuration(double d, int maxElements) {
		d = Math.round(d);
		if (d == 0.0)
			return "< 1s";

		StringBuilder sb = new StringBuilder();
		int numElements = 0;
		for (int i = 0; i < durations.length; i++) {
			int curUnit = (int) (d / durations[i]);
			if (curUnit > 0) {
				sb.append(curUnit).append(suffixes[i]).append(' ');
				if ((maxElements > 0) && (++numElements >= maxElements))
					break;
				d -= curUnit * durations[i];
			}
		}
		return sb.toString().trim();
	}

	private static double evalCubicBezier(double t, double start, double ctrl1, double ctrl2, double end) {
		double tInv = 1.0 - t;
		return (tInv * tInv * tInv) * start +
				3 * (tInv * tInv) * t * ctrl1 +
				3 * tInv * (t * t) * ctrl2 +
				(t * t * t) * end;
	}

	private static double[] getCubicBezierBounds(double start, double ctrl1, double ctrl2, double end) {
		double denom = -start + 3 * ctrl1 - 3 * ctrl2 + end;
		double p2 = (start - 2 * ctrl1 + ctrl2) / denom;
		double q = (ctrl1 - start) / denom;
		double t1 = -p2 + Math.sqrt((p2 * p2) - q);
		double t2 = -p2 - Math.sqrt((p2 * p2) - q);
		return new double[] {
				evalCubicBezier(t1, start, ctrl1, ctrl2, end), evalCubicBezier(t2, start, ctrl1, ctrl2, end) };
	}

	private static final double minInitialProb = 0.1;
	private static final double minTransitionProb = 0.1;
	private static final double minPrintProb = 0.01;

	private static final double minStrokeWidth = 1.0, maxStrokeWidth = 5.0;
	private static final double baselineY = 180.0;
	private static final double singleRadius = 15.0, doubleRadius = singleRadius + 3 * minStrokeWidth;
	private static final double heightScale = 20.0;
	private static final double radiusScaleFrom = 0.92;
	private static final double radiusScaleTo = 1.1;
	private static final int labelFontSize = 10;

	public static void plotTransitions(OutputStream os, InteractionHmm hmm, double minTransitionProb,
			boolean ignoreInitialProbabilities, boolean ignoreOutputProbabilities) throws IOException {
		Writer w = new OutputStreamWriter(os, Charset.forName("UTF-8"));
		DecimalFormat svgFmt = (DecimalFormat) NumberFormat.getInstance(Locale.US);
		svgFmt.applyPattern("0.0######");

		w.write(svgHeader);

		int numStates = hmm.nbStates();
		boolean[] isInitialState = new boolean[numStates];
		if (!ignoreInitialProbabilities)
			for (int i = 0; i < isInitialState.length; i++)
				isInitialState[i] = (hmm.getPi(i) >= minInitialProb);

		// compute order and position of states
		final double[] statePos = new double[numStates];
		double maxStatePos = Double.NEGATIVE_INFINITY;
		if (!ignoreOutputProbabilities) {
			for (int i = 0; i < statePos.length; i++) {
				if (hmm.hasLogNormalEmissions())
					statePos[i] = Math.log10(((OpdfLogNormal) hmm.getOpdf(i)).arithmeticMean());
				else
					statePos[i] = Math.log10(((OpdfGaussian) hmm.getOpdf(i)).mean());
				if (statePos[i] > maxStatePos)
					maxStatePos = statePos[i];
			}
		} else {
			for (int i = 0; i < numStates; i++)
				statePos[i] = i + 1.0;
			maxStatePos = statePos[numStates - 1];
		}
		Integer[] posToStateIdx = new Integer[numStates];
		for (int i = 0; i < posToStateIdx.length; i++)
			posToStateIdx[i] = i;
		Arrays.sort(posToStateIdx, new Comparator<Integer>() {
			@Override
			public int compare(Integer o1, Integer o2) {
				return Double.compare(statePos[o1], statePos[o2]);
			}
		});

		double effectiveBorder = 5 * doubleRadius;
		double effectiveWidth = 1200.0 - (2 * effectiveBorder);
		double[] statePosX = new double[numStates];
		for (int i = 0; i < numStates; i++)
			statePosX[i] = effectiveBorder + ((statePos[posToStateIdx[i]] / maxStatePos) * effectiveWidth);

		// transition edge layout
		int[][] outgoingIdx = new int[numStates][numStates];
		int[][] numOutgoingEdges = new int[2][numStates];
		int[][] incomingIdx = new int[numStates][numStates];
		int[][] numIncomingEdges = new int[2][numStates];
		// visit target states in order of increasing distance
		for (int d = 1; d < numStates; d++) {
			for (int i = 0; i < numStates; i++) {
				int j = i + d;
				if ((j < numStates) && (hmm.getAij(posToStateIdx[i], posToStateIdx[j]) > minTransitionProb)) {
					outgoingIdx[i][j] = numOutgoingEdges[0][i]++;
					incomingIdx[i][j] = numIncomingEdges[0][j]++;
				}
				j = i - d;
				if ((j >= 0) && (hmm.getAij(posToStateIdx[i], posToStateIdx[j]) > minTransitionProb)) {
					outgoingIdx[i][j] = numOutgoingEdges[1][i]++;
					incomingIdx[i][j] = numIncomingEdges[1][j]++;
				}
			}
		}

		// draw transitions
		int maxHeightIdx = 0;
		StringBuilder labelDefs = new StringBuilder();
		for (int i = 0; i < numStates; i++) {
			int stateI = posToStateIdx[i];
			double fromR = (isInitialState[stateI] ? doubleRadius : singleRadius) * radiusScaleFrom;
			for (int j = 0; j < numStates; j++) {
				int stateJ = posToStateIdx[j];
				if (hmm.getAij(stateI, stateJ) > minTransitionProb) {
					double toR = (isInitialState[stateJ] ? doubleRadius : singleRadius) * radiusScaleTo;
					double strokeWidth = minStrokeWidth + (hmm.getAij(stateI, stateJ) *
							(maxStrokeWidth - minStrokeWidth));

					double fromPosX, fromPosY;
					double toPosX, toPosY;
					double ctrlPos1X, ctrlPos1Y, ctrlPos2X, ctrlPos2Y;
					double labelX, labelY;
					String labelAlign, labelBaseline;
					if (i != j) {
						// regular transitions
						int fromIdx = outgoingIdx[i][j];
						int toIdx = incomingIdx[i][j];
						double fromRad, toRad, midPosY;
						if (i < j) {	// above
							fromRad = (2.0 * Math.PI) -
									((fromIdx + 1) * ((0.5 * Math.PI) / (numOutgoingEdges[0][i] + 2)));
							toRad = Math.PI + ((toIdx + 1) * ((0.5 * Math.PI) / (numIncomingEdges[0][j] + 2)));
							midPosY = -1.0;
						} else {	// below
							fromRad = Math.PI -
									((fromIdx + 1) * ((0.5 * Math.PI) / (numOutgoingEdges[1][i] + 2)));
							toRad = ((toIdx + 1) * ((0.5 * Math.PI) / (numIncomingEdges[1][j] + 2)));
							midPosY = 1.0;
						}
						fromPosX = statePosX[i] + (fromR * Math.cos(fromRad));
						fromPosY = baselineY + (fromR * Math.sin(fromRad));
						toPosX = statePosX[j] + (toR * Math.cos(toRad));
						toPosY = baselineY + (toR * Math.sin(toRad));

						int heightIdx = Math.max(fromIdx, toIdx) + 1;
						if (heightIdx > maxHeightIdx)
							maxHeightIdx = heightIdx;
						double midOffsetX = (toPosX - fromPosX) * Math.pow(0.5, heightIdx);
						midPosY = fromPosY + (midPosY * (heightIdx * heightScale));
						ctrlPos1X = fromPosX + midOffsetX;
						ctrlPos2X = toPosX - midOffsetX;
						ctrlPos1Y = ctrlPos2Y = midPosY;

						double[] bounds = getCubicBezierBounds(fromPosY, ctrlPos1Y, ctrlPos2Y, toPosY);
						labelX = fromPosX + ((toPosX - fromPosX) / 2);
						if (i < j) {
							labelY = Math.min(bounds[0], bounds[1]);
							labelBaseline = "top";
						} else {
							labelY = Math.max(bounds[0], bounds[1]);
							labelBaseline = "bottom";
						}
						labelAlign = "middle";
					} else {
						// self transitions
						double fromRad = -0.05 * Math.PI;
						double toRad = 0.05 * Math.PI;
						fromPosX = statePosX[i] + (fromR * Math.cos(fromRad));
						fromPosY = baselineY + (fromR * Math.sin(fromRad));
						toPosX = statePosX[j] + (toR * Math.cos(toRad));
						toPosY = baselineY + (toR * Math.sin(toRad));
						ctrlPos1X = ctrlPos2X = fromPosX + toR + heightScale;
						ctrlPos1Y = fromPosY - heightScale;
						ctrlPos2Y = toPosY + heightScale;

						labelX = fromPosX + fromR + (2 * strokeWidth);
						labelY = fromPosY + ((toPosY - fromPosY) / 2);
						labelBaseline = "middle";
						labelAlign = "left";
					}

 					w.write("  <path style=\"fill:none;stroke:#000000;marker-end:url(#arrow);stroke-width:" +
							svgFmt.format(strokeWidth) +
							"\" d=\"M " + svgFmt.format(fromPosX) + "," + svgFmt.format(fromPosY) +
							" C " + svgFmt.format(ctrlPos1X) + "," + svgFmt.format(ctrlPos1Y) +
							" " + svgFmt.format(ctrlPos2X) + "," + svgFmt.format(ctrlPos2Y) +
							" " + svgFmt.format(toPosX) + "," + svgFmt.format(toPosY) + "\" />\n");

					// transition probability label
					labelDefs.append("  <text x=\"" + svgFmt.format(labelX) + "\" y=\"" + svgFmt.format(labelY) +
							"\" text-anchor=\"" + labelAlign + "\" dominant-baseline=\"" + labelBaseline + "\" " +
							"fill=\"#000000\" stroke=\"#ffffff\" " +
							"stroke-width=\"3\" paint-order=\"stroke\" font-family=\"Arial\" font-size=\"" +
							labelFontSize + "pt\">" +
							LogProbabilityAxisRenderer.formatProbability(hmm.getAij(stateI, stateJ), minPrintProb,
									false) + "</text>\n");
				}
			}
		}
		// labels should not be covered by paths, so they need to be drawn afterwards
		w.write(labelDefs.toString());

		// draw states
		boolean prevLabelShifted = false;
		for (int i = 0; i < numStates; i++) {
			int stateI = posToStateIdx[i];
			String circlePrefix = "  <circle cx=\"" + svgFmt.format(statePosX[i]) +
					"\" cy=\"" + svgFmt.format(baselineY) +
					"\" style=\"fill:#ffffff;stroke:#000000;stroke-width:" + svgFmt.format(minStrokeWidth) + "\" r=\"";
			String circleSuffix = "\" />\n";
			if (isInitialState[stateI])
				w.write(circlePrefix + svgFmt.format(doubleRadius) + circleSuffix);
			w.write(circlePrefix + svgFmt.format(singleRadius) + circleSuffix);
			w.write("  <text x=\"" + svgFmt.format(statePosX[i]) + "\" y=\"" + svgFmt.format(baselineY) + "\" " +
					"text-anchor=\"middle\" dominant-baseline=\"central\" " +
					"font-family=\"Arial\" font-weight=\"bold\" font-size=\"12pt\">" + (i + 1) + "</text>\n");

			// initial state probability
			double topBaselineY = baselineY - ((maxHeightIdx + 4) * heightScale);
			if (!prevLabelShifted && (i > 0) && ((statePosX[i] - statePosX[i - 1]) < (3 * doubleRadius))) {
				topBaselineY += (3 * (labelFontSize + 3)) + 2;
				prevLabelShifted = true;
			} else
				prevLabelShifted = false;
			if (!ignoreInitialProbabilities) {
				w.write("  <text x=\"" + svgFmt.format(statePosX[i]) + "\" y=\"" + svgFmt.format(topBaselineY) +
						"\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"" + labelFontSize + "pt\">");
				if (i == 0)
					w.write("Pi = ");
				w.write(LogProbabilityAxisRenderer.formatProbability(hmm.getPi(stateI), minPrintProb, false) +
						"</text>\n");
				topBaselineY += labelFontSize + 3;
			}

			// output mean/std.dev
			if (!ignoreOutputProbabilities) {
				double mu, sigma;
				if (hmm.hasLogNormalEmissions()) {
					OpdfLogNormal opdf = (OpdfLogNormal) hmm.getOpdf(stateI);
					mu = opdf.arithmeticMean();
					sigma = Math.sqrt(opdf.arithmeticVariance());
				} else {
					OpdfGaussian opdf = (OpdfGaussian) hmm.getOpdf(stateI);
					mu = opdf.mean();
					sigma = Math.sqrt(opdf.variance());
				}
				w.write("  <text x=\"" + svgFmt.format(statePosX[i]) + "\" y=\"" + svgFmt.format(topBaselineY) +
						"\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"" + labelFontSize + "pt\">");
				if (i == 0)
					w.write("\u03bc = ");
				w.write(formatDuration(mu, 2) + "</text>\n");
				topBaselineY += labelFontSize + 3;
				w.write("  <text x=\"" + svgFmt.format(statePosX[i]) + "\" y=\"" + svgFmt.format(topBaselineY) +
						"\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"" + labelFontSize + "pt\">");
				if (i == 0)
					w.write("\u03c3 = ");
				w.write(formatDuration(sigma, 2) + "</text>\n");
			}
		}

		w.write(svgFooter);
		w.close();
	}

	private static double evalGaussian(double x, OpdfGaussian g) {
		return (1.0 / Math.sqrt(2.0 * g.variance() * Math.PI)) *
				Math.exp(-((x - g.mean()) * (x - g.mean())) / (2.0 * g.variance()));
	}

	private static double evalLogNormal(double x, OpdfLogNormal d) {
		if (x <= 0.0)
			return 0.0;
		double x1 = (Math.log(x) - d.mean()) / Math.sqrt(d.variance());
		return Math.exp(-0.5 * x1 * x1) / (Math.sqrt(2.0 * d.variance() * Math.PI) * x);
	}

	private static double getGaussianBound(double minY, double mean, double variance) {
		return mean + Math.sqrt(-Math.log(minY * Math.sqrt(2.0 * variance * Math.PI)) * (2.0 * variance));
	}

	private static double getGaussianBound(double minY, OpdfGaussian g) {
		return getGaussianBound(minY, g.mean(), g.variance());
	}

	private static double getLogNormalBound(double minY, OpdfLogNormal d) {
		return Math.exp(getGaussianBound(minY, d.mean(), d.variance()));
	}

	private static final int numPoints = 500;

	private static XYPlot prepareGaussianPlot(InteractionHmm hmm) {
		int numStates = hmm.nbStates();

		// compute range of y axis
		DataTable modePoints = new DataTable(2, Double.class);
		int minYExp = 0, maxYExp = Integer.MIN_VALUE;
		for (int i = 0; i < numStates; i++) {
			double mu = ((OpdfGaussian) hmm.getOpdf(i)).mean();		// mean == mode
			double y = evalGaussian(mu, (OpdfGaussian) hmm.getOpdf(i));
			modePoints.add(mu, y);
			double logY = Math.log10(y);
			int yExp = (int) Math.floor(logY);
			if (yExp < minYExp)
				minYExp = yExp;
			yExp = (int) Math.ceil(logY);
			if (yExp > maxYExp)
				maxYExp = yExp;
		}
		minYExp--;
		double minY = Math.pow(10.0, minYExp);
		double maxY = Math.pow(10.0, maxYExp);

		// compute range of x axis
		int maxXExp = 0;
		for (int i = 0; i < numStates; i++) {
			int xExp = (int) Math.ceil(Math.log10(getGaussianBound(minY, (OpdfGaussian) hmm.getOpdf(i))));
			if (xExp > maxXExp)
				maxXExp = xExp;
		}

		// fill data table with samples from PDF
		DataTable data = new DataTable(numStates + 1, Double.class);
		double xExp = 0.0;
		double inc = (double) maxXExp / numPoints;
		for (int i = 0; i < numPoints; i++) {
			Double[] v = new Double[numStates + 1];
			v[0] = Math.pow(10.0, xExp);
			for (int j = 0; j < numStates; j++)
				v[j + 1] = Math.max(minY, evalGaussian(v[0], (OpdfGaussian) hmm.getOpdf(j)));
			data.add(v);
			xExp += inc;
		}

		XYPlot plot = new XYPlot();
		for (int i = 0; i < numStates; i++) {
			DataSeries s = new DataSeries(data, 0, i + 1);
			plot.add(s);
			plot.setPointRenderers(s, null);
			plot.setLineRenderers(s, new DefaultLineRenderer2D());
		}
		plot.add(modePoints);
		plot.setLineRenderers(modePoints, null);

		// setup log-log plot
		plot.getAxis(XYPlot.AXIS_X).setRange(1.0, Math.pow(10.0, maxXExp));
		plot.getAxis(XYPlot.AXIS_Y).setRange(minY, maxY + (maxY / 100));
		plot.setAxisRenderer(XYPlot.AXIS_Y, new LogProbabilityAxisRenderer(minY));
		return plot;
	}

	private static XYPlot prepareLogNormalPlot(InteractionHmm hmm) {
		int numStates = hmm.nbStates();

		// compute scale factor for each distribution
		double[] scale = new double[numStates];
		for (int i = 0; i < numStates; i++) {
			OpdfLogNormal dist = (OpdfLogNormal) hmm.getOpdf(i);
			double mu = Math.exp(dist.mean() - dist.variance());	// mode of log-normal distribution
			scale[i] = 1.0 / evalLogNormal(mu, dist);
		}

		// compute range of x axis
		int maxXExp = 0;
		for (int i = 0; i < numStates; i++) {
			int xExp = (int) Math.ceil(Math.log10(getLogNormalBound(1E-1, (OpdfLogNormal) hmm.getOpdf(i))));
			if (xExp > maxXExp)
				maxXExp = xExp;
		}

		// fill data table with samples from PDF
		DataTable data = new DataTable(numStates + 1, Double.class);
		double xExp = 0.0;
		double inc = (double) maxXExp / numPoints;
		for (int i = 0; i < numPoints; i++) {
			Double[] v = new Double[numStates + 1];
			v[0] = Math.pow(10.0, xExp);
			for (int j = 0; j < numStates; j++)
				v[j + 1] = Math.max(1E-7, scale[j] * evalLogNormal(v[0], (OpdfLogNormal) hmm.getOpdf(j)));
			data.add(v);
			xExp += inc;
		}

		XYPlot plot = new XYPlot();
		for (int i = 0; i < numStates; i++) {
			DataSeries s = new DataSeries(data, 0, i + 1);
			plot.add(s);
			plot.setPointRenderers(s, null);
			plot.setLineRenderers(s, new DefaultLineRenderer2D());
		}

		// setup log-linear plot
		plot.getAxis(XYPlot.AXIS_X).setRange(1.0, Math.pow(10.0, maxXExp));
		plot.setAxisRenderer(XYPlot.AXIS_Y, new LinearRenderer2D());
		plot.getAxis(XYPlot.AXIS_Y).setRange(0.0, 1.02);
		plot.getAxisRenderer(XYPlot.AXIS_Y).setTickLabelsVisible(false);
		return plot;
	}

	public static void plotOutput(OutputStream os, InteractionHmm hmm) throws IOException {
		XYPlot plot;
		if (hmm.hasLogNormalEmissions())
			plot = prepareLogNormalPlot(hmm);
		else
			plot = prepareGaussianPlot(hmm);

		plot.setInsets(new Insets2D.Double(10.0, 70.0, 30.0, 10.0));
		plot.setAxisRenderer(XYPlot.AXIS_X, new LogTimeAxisRenderer());

		// render chart as SVG
		DrawableWriter dw = DrawableWriterFactory.getInstance().get("image/svg+xml");
		dw.write(plot, os, 1200, hmm.hasLogNormalEmissions() ? 250 : 400);
	}

	private static void computeExpectedBurstLengths(final InteractionHmm hmm, boolean ignoreInitialProbabilities) {
		int n = hmm.nbStates();

		// sort states by mean of output distribution and build permuted transition matrix
		Integer[] posToStateIdx = new Integer[n];
		for (int i = 0; i < posToStateIdx.length; i++)
			posToStateIdx[i] = i;
		Arrays.sort(posToStateIdx, new Comparator<Integer>() {
			@Override
			public int compare(Integer o1, Integer o2) {
				if (hmm.hasLogNormalEmissions()) {
					return Double.compare(((OpdfLogNormal) hmm.getOpdf(o1)).arithmeticMean(),
							((OpdfLogNormal) hmm.getOpdf(o2)).arithmeticMean());
				}
				return Double.compare(((OpdfGaussian) hmm.getOpdf(o1)).mean(),
						((OpdfGaussian) hmm.getOpdf(o2)).mean());
			}
		});

		RealMatrix a = MatrixUtils.createRealMatrix(n, n);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				a.setEntry(i, j, hmm.getAij(posToStateIdx[i], posToStateIdx[j]));

		for (int i = 0; i < (n - 1); i++) {
			// Construct an absorbing Markov chain by merging all states > i into a single absorbing state and compute
			// the expectation of the sum of outputs associated with the non-absorbing states.
			String comment = "";

			// Computing the expected counts of visiting the non-absorbing states is equivalent to computing
			// Pi^T * (I - Q)^-1, where Pi is the column vector of initial probabilities (re-normalized), I is the
			// identity matrix and Q is the quadratic sub-matrix of size (i + 1) of the state transition matrix. This in
			// turn is equivalent to solving (I - Q)^T * x = Pi.
			RealMatrix aa = a.getSubMatrix(0, i, 0, i);
			RealMatrix pi = MatrixUtils.createRealMatrix((i + 1), 1);
			if (!ignoreInitialProbabilities) {
				double sum = 0.0;
				for (int j = 0; j < (i + 1); j++) {
					pi.setEntry(j, 0, hmm.getPi(posToStateIdx[j]));
					sum += hmm.getPi(posToStateIdx[j]);
				}
				if (sum == 0.0) {
					comment += "all states have 0 initial probability";
					sum = 1.0;
				}
				for (int j = 0; j < (i + 1); j++)
					pi.setEntry(j, 0, pi.getEntry(j, 0) / sum);
			} else {
				for (int j = 0; j < (i + 1); j++)
					pi.setEntry(j, 0, 1.0 / (i + 1));
			}

			RealMatrix m = MatrixUtils.createRealIdentityMatrix(i + 1).subtract(aa).transpose();
			LUDecomposition lu = new LUDecomposition(m);
			RealMatrix x = lu.getSolver().solve(pi);

			// The expected visit counts are multiplied by the expectation of the output distribution, i.e. its mean,
			// and summed up (linearity of expectation).
			double c = 0.0;	// expected number of events
			double t = 0.0;	// expected duration
			for (int j = 0; j < (i + 1); j++) {
				c += x.getEntry(j, 0);
				double out;
				if (hmm.hasLogNormalEmissions()) {
					out = ((OpdfLogNormal) hmm.getOpdf(posToStateIdx[j])).arithmeticMean();
				} else
					out = ((OpdfGaussian) hmm.getOpdf(posToStateIdx[j])).mean();
				t += out * x.getEntry(j, 0);
			}

			System.out.println(i + ";" + t + ";" + HmmPlotter.formatDuration(t, 0) + ";" + c + ";" + comment);
		}
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 1) {
			System.err.println("usage: " + HmmPlotter.class.getSimpleName() + " hmm [out] [ignore p_i?]");
			return;
		}
		String outputDir = ".";
		if (args.length > 1)
			outputDir = args[1];
		boolean ignoreInitialProbabilities = false;
		if (args.length > 2)
			ignoreInitialProbabilities = Boolean.parseBoolean(args[2]);

		InteractionHmm hmm = InteractionHmm.read(new File(args[0]));

		computeExpectedBurstLengths(hmm, ignoreInitialProbabilities);

		Writer w = new FileWriter(new File(outputDir, "hmm-param.txt"));
		try {
			if (hmm.hasLogNormalEmissions())
				HmmWriter.write(w, new OpdfLogNormalWriter(), hmm);
			else
				HmmWriter.write(w, new OpdfGaussianWriter(), hmm);
		} finally {
			w.close();
		}

		OutputStream os = new FileOutputStream(new File(outputDir, "hmm-trans.svg"));
		try {
			plotTransitions(os, hmm, minTransitionProb, ignoreInitialProbabilities, false);
		} finally {
			os.close();
		}

		os = new FileOutputStream(new File(outputDir, "hmm-out.svg"));
		try {
			plotOutput(os, hmm);
		} finally {
			os.close();
		}
	}

}
