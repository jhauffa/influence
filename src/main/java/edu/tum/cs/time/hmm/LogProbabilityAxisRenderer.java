package edu.tum.cs.time.hmm;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Set;

import de.erichseifert.gral.plots.axes.Axis;
import de.erichseifert.gral.plots.axes.LogarithmicRenderer2D;
import de.erichseifert.gral.plots.axes.Tick;
import de.erichseifert.gral.util.PointND;

public class LogProbabilityAxisRenderer extends LogarithmicRenderer2D {

	private static final long serialVersionUID = 1L;

	private final double minY;

	public LogProbabilityAxisRenderer(double minY) {
		super();
		this.minY = minY;
	}

	@Override
	protected void createTicks(List<Tick> ticks, Axis axis, double min, double max, Set<Double> tickPositions,
			boolean isAutoSpacing) {
		List<Tick> defaultTicks = new ArrayList<Tick>();
		super.createTicks(defaultTicks, axis, min, max, tickPositions, isAutoSpacing);
		int idx = 0;
		double p = minY;
		PointND<Double> prevPos = null;
		for (Tick tick : defaultTicks) {
			if ((tick.type == Tick.TickType.MAJOR) && (tick.position != null) && (tick.normal != null)) {
				if ((prevPos != null) && prevPos.equals(tick.position))
					continue;
				String label = "";
				if ((idx++ % 9) == 0) {
					label = formatProbability(p, Double.MIN_VALUE, true);
					p *= 10;
				} else
					tick.normal.set(PointND.X, tick.normal.get(PointND.X) * 0.6);
				ticks.add(new Tick(tick.type, tick.position, tick.normal, tick.drawable, tick.shape, label));
				prevPos = tick.position;
			}
		}
	}

	public static String formatProbability(double p, double minValue, boolean scientificNotation) {
		DecimalFormat fmt = (DecimalFormat) NumberFormat.getInstance(Locale.US);
		String pattern = "0.0##";
		if (scientificNotation && (p < 0.001))
			pattern += "E0";
		fmt.applyPattern(pattern);
		if (p < minValue)
			return "&lt; " + fmt.format(minValue);
		return fmt.format(p);
	}

}
