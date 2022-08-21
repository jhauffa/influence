package edu.tum.cs.db.stats;

import java.util.Date;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.SocialMediaDaoFactory;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.nlp.tokenizer.TwitterTokenizer;
import edu.tum.cs.time.IntervalCalculator;
import edu.tum.cs.util.ExperimentConfiguration;

public class TokenStatistics {

	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(TokenStatistics.class);

	public static void main(String[] args) throws Exception {
		boolean coreOnly = true;
		if (args.length > 0)
			coreOnly = Boolean.parseBoolean(args[0]);

		Date endDate = cfg.getDateProperty(ExperimentConfiguration.PROP_END_DATE);
		IntervalCalculator ic = new IntervalCalculator(endDate);
		Date startDate = ic.getStartDate();

		TwitterTokenizer tokenizer = new TwitterTokenizer();

		int numMessages = 0;
		int numMessagesInObsPeriod = 0;
		int numTokens = 0;
		int numTokensInObsPeriod = 0;
		int numUsers = 0;
		DescriptiveStatistics tokenStats = new DescriptiveStatistics();
		DescriptiveStatistics tokenStatsInObsPeriod = new DescriptiveStatistics();

		SocialMediaDao dao = SocialMediaDaoFactory.createDao();
		for (long userId : dao.getUserIds(coreOnly)) {
			for (SocialMediaMessage message : dao.getMessages(userId, null, null)) {
				List<String> words = tokenizer.tokenize(message.getText(), true);
				int n = words.size();
				numTokens += n;
				tokenStats.addValue(n);

				if (message.getDate().after(startDate) && message.getDate().before(endDate)) {
					numTokensInObsPeriod += n;
					tokenStatsInObsPeriod.addValue(n);
					numMessagesInObsPeriod++;
				}

				if ((++numMessages % 1000000) == 0)
					System.err.println("processed " + numMessages + " messages");
			}

			if ((++numUsers % 100) == 0)
				System.err.println("processed " + numUsers + " users");
		}
		System.out.println("finished processing " + numUsers + " users");
		System.out.println();

		System.out.println("all messages:");
		System.out.println("n = " + numMessages);
		System.out.println("#tokens = " + numTokens);
		System.out.println(tokenStats);

		System.out.println("messages in observation period (" + startDate + " - " + endDate + "):");
		System.out.println("n = " + numMessagesInObsPeriod);
		System.out.println("#tokens = " + numTokensInObsPeriod);
		System.out.println(tokenStatsInObsPeriod);
	}

}
