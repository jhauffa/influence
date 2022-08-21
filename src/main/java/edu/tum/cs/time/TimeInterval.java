package edu.tum.cs.time;

import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;

public class TimeInterval {

	private static final SimpleDateFormat df = new SimpleDateFormat("yyyy.MM.dd");

	private final Date start;
	private final int length;

	public TimeInterval(Date start, int length) {
		if (length >= 0) {
			this.start = start;
			this.length = length;
		} else {
			this.start = shift(start, length);
			this.length = -length;
		}
	}

	private static Date shift(Date center, int offset) {
		Calendar cal = new GregorianCalendar();
		cal.setTime(center);
		cal.add(Calendar.DAY_OF_YEAR, offset);
		return cal.getTime();
	}

	public Date getStartDate() {
		return start;
	}

	public Date getEndDate() {
		return shift(start, length);
	}

	public long key() {
		assert length < 256;
		return (start.getTime() << 8) | length;
	}

	@Override
	public String toString() {
		return df.format(start) + ";" + length;
	}

	@Override
	public int hashCode() {
		return 31 * (31 + start.hashCode()) + length;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof TimeInterval))
			return false;
		TimeInterval other = (TimeInterval) obj;
		return ((start.equals(other.start)) && (length == other.length));
	}

}
