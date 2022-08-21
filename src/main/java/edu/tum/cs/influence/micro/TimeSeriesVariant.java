package edu.tum.cs.influence.micro;

import java.io.Serializable;
import java.util.Queue;

public class TimeSeriesVariant implements Serializable {

	private static final long serialVersionUID = 5102699619361504848L;

	private final int timePeriodLength; 	// days

	public TimeSeriesVariant(int timePeriodLength) {
		this.timePeriodLength = timePeriodLength;
	}

	public int getTimePeriodLength() {
		return timePeriodLength;
	}

	@Override
	public int hashCode() {
		return Integer.valueOf(timePeriodLength).hashCode();
	}

	@Override
	public boolean equals(Object o) {
		if (o == this)
			return true;
		if (!(o instanceof TimeSeriesVariant))
			return false;
		TimeSeriesVariant other = (TimeSeriesVariant) o;
		return (timePeriodLength == other.timePeriodLength);
	}

	@Override
	public String toString() {
		return Integer.toString(timePeriodLength);
	}

	public static String getCsvHeader() {
		return "timePeriodLength";
	}

	public static TimeSeriesVariant readCsv(Queue<String> parts) {
		return new TimeSeriesVariant(Integer.parseInt(parts.poll()));
	}

}
