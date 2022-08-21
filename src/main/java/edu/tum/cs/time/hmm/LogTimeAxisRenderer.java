package edu.tum.cs.time.hmm;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import de.erichseifert.gral.plots.axes.Axis;
import de.erichseifert.gral.plots.axes.LogarithmicRenderer2D;
import de.erichseifert.gral.plots.axes.Tick;
import de.erichseifert.gral.util.PointND;

public class LogTimeAxisRenderer extends LogarithmicRenderer2D {

	private static final long serialVersionUID = 1L;

	public LogTimeAxisRenderer() {
		super();
		setTickLabelDistance(getLabelDistance() - 0.4);
		Map<Double, String> ticks = new HashMap<Double, String>();
		ticks.put(0.001, "1ms");
		ticks.put(1.0, "1s");
		ticks.put(60.0, "1m");
		ticks.put(60 * 60.0, "1h");
		ticks.put(24 * 60 * 60.0, "1d");
		ticks.put(30 * 24 * 60 * 60.0, "1M");
		ticks.put(365 * 24 * 60 * 60.0, "1y");
		ticks.put(10 * 365 * 24 * 60 * 60.0, "10y");
		setCustomTicks(ticks);
	}

	@Override
	protected void createTicksCustom(List<Tick> ticks, Axis axis, double min, double max, Set<Double> tickPositions) {
		// do nothing, custom ticks will be inserted later
	}

	@Override
	protected void createTicks(List<Tick> ticks, Axis axis, double min, double max, Set<Double> tickPositions,
			boolean isAutoSpacing) {
		List<Tick> defaultTicks = new ArrayList<Tick>();
		super.createTicks(defaultTicks, axis, min, max, tickPositions, isAutoSpacing);
		int idx = 0;
		for (Tick tick : defaultTicks) {
			if ((tick.type == Tick.TickType.MAJOR) && (tick.position != null) && (tick.normal != null)) {
				if ((idx++ % 9) != 0)
					tick.normal.set(PointND.Y, tick.normal.get(PointND.Y) * 0.6);
				ticks.add(new Tick(tick.type, tick.position, tick.normal, tick.drawable, tick.shape, ""));
			}
		}
		super.createTicksCustom(ticks, axis, min, max, tickPositions);
	}

}
