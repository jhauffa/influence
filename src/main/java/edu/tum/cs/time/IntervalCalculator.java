package edu.tum.cs.time;

import java.text.ParseException;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;

import edu.tum.cs.util.ExperimentConfiguration;

public class IntervalCalculator {

	private static final ExperimentConfiguration cfg;
	private static final int numIntervals, weeksPerInterval;
	static {
		cfg = new ExperimentConfiguration(IntervalCalculator.class);
		numIntervals = cfg.getLocalIntProperty("numIntervals");
		weeksPerInterval = cfg.getLocalIntProperty("weeksPerInterval");
	}

	private final Calendar cal = new GregorianCalendar();
	private final Date endDate;

	public IntervalCalculator(Date endDate) {
		this.endDate = endDate;
	}

	public Date getStartDate() throws ParseException {
		int globalLength = numIntervals * weeksPerInterval;
		cal.setTime(endDate);
		cal.add(Calendar.WEEK_OF_MONTH, -globalLength);
		return cal.getTime();
	}

	public Date[] getCenterDates() throws ParseException {
		cal.setTime(endDate);
		Date[] dates = new Date[numIntervals - 1];
		for (int i = 0; i < dates.length; i++) {
			cal.add(Calendar.WEEK_OF_MONTH, -weeksPerInterval);
			dates[i] = cal.getTime();
		}
		return dates;
	}

}
