package edu.tum.cs.db.stats;

import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.Map;
import java.util.TreeMap;

import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.SocialMediaDaoFactory;
import edu.tum.cs.db.TwitterDao;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.math.dist.HistogramImpl.IdentityHistogram;

public class MessageVolumeHisto {

	private static final SimpleDateFormat sdf = new SimpleDateFormat("yyyy.MM.dd");
	private static final Date firstTweetDate = new Date(1142954693);	// slightly before the first available tweet

	private static void printHisto(String caption, IdentityHistogram histo) {
		System.out.println(caption);
		System.out.println(histo);
		System.out.println();
	}

	public static void main(String[] args) throws Exception {
		boolean coreOnly = true, detectBackdating = false;
		Date startDate = null, endDate = null;
		if (args.length > 0) {
			startDate = sdf.parse(args[0]);
			if (args.length > 1) {
				endDate = sdf.parse(args[1]);
				if (args.length > 2) {
					coreOnly = Boolean.parseBoolean(args[2]);
					if (args.length > 3)
						detectBackdating = Boolean.parseBoolean(args[3]);
				}
			}
		}

		Map<Integer, IdentityHistogram> indivYearHisto = new TreeMap<Integer, IdentityHistogram>();
		IdentityHistogram aggYearHisto = new IdentityHistogram();
		IdentityHistogram weekHisto = new IdentityHistogram();
		IdentityHistogram dayHisto = new IdentityHistogram();

		int numMessages = 0;
		int numUsers = 0;
		int numSkipped = 0;

		Calendar cal = GregorianCalendar.getInstance();
		SocialMediaDao dao = SocialMediaDaoFactory.createDao();
		boolean isTwitter = dao instanceof TwitterDao;
		for (long userId : dao.getUserIds(coreOnly)) {
			for (SocialMediaMessage message : dao.getMessages(userId, startDate, endDate)) {
				Date d = message.getDate();
				if (detectBackdating) {
					if (isTwitter) {
						cal.setTime(TwitterDao.fixMessageDate(message.getId(), d, firstTweetDate));
					} else {
						cal.setTime(d);
						if ((cal.get(Calendar.SECOND) == 0) &&
						    ((cal.get(Calendar.MINUTE) == 0) || (cal.get(Calendar.MINUTE) == 30))) {
							numSkipped++;
							continue;
						}
					}
				} else
					cal.setTime(d);

				int year = cal.get(Calendar.YEAR);
				IdentityHistogram yearHisto = indivYearHisto.get(year);
				if (yearHisto == null) {
					yearHisto = new IdentityHistogram();
					indivYearHisto.put(year, yearHisto);
				}
				int month = cal.get(Calendar.MONTH);
				yearHisto.addValue(month);
				aggYearHisto.addValue(month);
				weekHisto.addValue(cal.get(Calendar.DAY_OF_WEEK));
				dayHisto.addValue(cal.get(Calendar.HOUR_OF_DAY));

				if ((++numMessages % 1000000) == 0)
					System.err.println("processed " + numMessages + " messages");
			}

			if ((++numUsers % 100) == 0)
				System.err.println("processed " + numUsers + " users");
		}

		for (Map.Entry<Integer, IdentityHistogram> e : indivYearHisto.entrySet())
			printHisto("year " + e.getKey(), e.getValue());
		printHisto("yearly", aggYearHisto);
		printHisto("weekly", weekHisto);
		printHisto("daily", dayHisto);
		System.out.println("total: " + numMessages + " messages");
		if (numSkipped > 0)
			System.out.println("ignored " + numSkipped + " potentially backdated messages");
	}

}
