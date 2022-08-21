package edu.tum.cs.db.stats;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Map;

public class EventJsonConverter {

	private static final boolean boundaryOnly = true;

	private static final void addToHisto(Map<Date, Integer> histo, Date date, int numEvents) {
		Integer n = histo.get(date);
		if (n == null)
			n = 0;
		histo.put(date, n + numEvents);
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 1) {
			System.err.println("usage: " + EventJsonConverter.class.getSimpleName() + " CSV");
			return;
		}

		// read event data from CSV file
		Map<Date, Integer> eventHisto = new HashMap<Date, Integer>();
		SimpleDateFormat df = new SimpleDateFormat("dd.MM.yyyy");
		Calendar cal = new GregorianCalendar();
		BufferedReader reader = new BufferedReader(new FileReader(args[0]));
		try {
			String line;
			while ((line = reader.readLine()) != null) {
				String[] parts = line.split(";");
				if (parts.length == 3) {
					Date startDate = df.parse(parts[0]);
					Date endDate = df.parse(parts[1]);
					int numTopics = Integer.parseInt(parts[2]);

					if (boundaryOnly) {
						addToHisto(eventHisto, startDate, numTopics);
						addToHisto(eventHisto, endDate, numTopics);
					} else {
						cal.setTime(startDate);
						while (!cal.getTime().after(endDate)) {
							addToHisto(eventHisto, cal.getTime(), numTopics);
							cal.add(Calendar.DAY_OF_YEAR, 1);
						}
					}
				}
			}
		} finally {
			reader.close();
		}

		// print JSON
		System.out.println("{");
		int idx = 0;
		for (Map.Entry<Date, Integer> e : eventHisto.entrySet()) {
			System.out.print("\t\"" + (e.getKey().getTime() / 1000L) + "\": " + e.getValue());
			if (++idx < eventHisto.size())
				System.out.println(",");
			else
				System.out.println();
		}
		System.out.println("}");
	}

}
