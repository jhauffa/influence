package edu.tum.cs.db.stats;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import edu.tum.cs.db.TwitterDao;
import edu.tum.cs.db.entities.SocialMediaMessage;

public class TweetTaxonomy {

	private static class TweetFeatures {
		private static final Pattern indicators = Pattern.compile("(?:(?:RT|via) )?@([a-zA-Z0-9_]{2,15})(?:\\s|$)");

		public static int numRetweetIndicatorRT = 0;
		public static int numRetweetIndicatorVia = 0;
		public static int numRetweetsOfSelf = 0;
		public static int numRepliesToSelf = 0;

		private final boolean hasRetweetUserId;
		private final boolean hasReplyUserId;
		private final boolean startsWithAtMention;
		private final int numAtMentions;
		private final boolean startsWithRetweetIndicator;
		private final int numRetweetIndicators;

		public TweetFeatures(boolean hasRetweetUserId, boolean hasReplyUserId, boolean startsWithAtMention,
				int numAtMentions, boolean startsWithRetweetIndicator, int numRetweetIndicators) {
			this.hasRetweetUserId = hasRetweetUserId;
			this.hasReplyUserId = hasReplyUserId;
			if (hasRetweetUserId && hasReplyUserId)
				throw new RuntimeException("tweet has both reply and retweet ID");
			this.startsWithAtMention = startsWithAtMention;
			this.numAtMentions = numAtMentions;
			this.startsWithRetweetIndicator = startsWithRetweetIndicator;
			this.numRetweetIndicators = numRetweetIndicators;
		}

		@Override
		public int hashCode() {
			int result = 31 + (hasReplyUserId ? 1231 : 1237);
			result = 31 * result + (hasRetweetUserId ? 1231 : 1237);
			result = 31 * result + numAtMentions;
			result = 31 * result + numRetweetIndicators;
			result = 31 * result + (startsWithAtMention ? 1231 : 1237);
			return 31 * result + (startsWithRetweetIndicator ? 1231 : 1237);
		}

		@Override
		public boolean equals(Object obj) {
			if (!(obj instanceof TweetFeatures))
				return false;
			TweetFeatures other = (TweetFeatures) obj;
			if ((hasReplyUserId == other.hasReplyUserId) && (hasRetweetUserId == other.hasRetweetUserId) &&
				(numAtMentions == other.numAtMentions) && (numRetweetIndicators == other.numRetweetIndicators) &&
				(startsWithAtMention == other.startsWithAtMention) &&
				(startsWithRetweetIndicator == other.startsWithRetweetIndicator))
				return true;
			return false;
		}

		@Override
		public String toString() {
			StringBuffer sb = new StringBuffer();
			if (hasReplyUserId)
				sb.append("has reply ID, ");
			else if (hasRetweetUserId)
				sb.append("has retweet ID, ");
			sb.append(numAtMentions).append(" @-mentions");
			if (startsWithAtMention)
				sb.append(" (starts with @-mention)");
			sb.append(", ").append(numRetweetIndicators).append(" retweet indicators");
			if (startsWithRetweetIndicator)
				sb.append(" (starts with retweet indicator)");
			return sb.toString();
		}

		public boolean isAddressive() {
			return (hasReplyUserId || (numAtMentions > 0));
		}

		public boolean isOriginal() {
			return !(hasRetweetUserId || (numRetweetIndicators > 0));
		}

		public static TweetFeatures classifyTweet(SocialMediaMessage tweet, long userId, String userName) {
			boolean hasRetweetUserId = false;
			if (tweet.getRepostOfUserId() > 0L) {
				hasRetweetUserId = true;
				if (tweet.getRepostOfUserId() == userId)
					numRetweetsOfSelf++;
			}
			boolean hasReplyUserId = false;
			if (tweet.getInReplyToUserId() > 0L) {
				hasReplyUserId = true;
				if (tweet.getInReplyToUserId() == userId)
					numRepliesToSelf++;
			}

			boolean startsWithAtMention = false;
			int numAtMentions = 0;
			boolean startsWithRetweetIndicator = false;
			int numRetweetIndicators = 0;
			Matcher m = indicators.matcher(tweet.getText());
			while (m.find()) {
				boolean startsWith = (m.start() == 0);
				if (m.group().charAt(0) == '@') {
					if (startsWith)
						startsWithAtMention = true;
					numAtMentions++;
					if (m.group(1).equals(userName))
						numRepliesToSelf++;
				} else {
					if (startsWith)
						startsWithRetweetIndicator = true;
					numRetweetIndicators++;
					if (m.group(1).equals(userName))
						numRetweetsOfSelf++;

					if (m.group().startsWith("RT"))
						numRetweetIndicatorRT++;
					else
						numRetweetIndicatorVia++;
				}
			}

			return new TweetFeatures(hasRetweetUserId, hasReplyUserId, startsWithAtMention, numAtMentions,
					startsWithRetweetIndicator, numRetweetIndicators);
		}
	}

	private static class AggregatedTweetFeatures {
		private final boolean hasRetweetUserId;
		private final boolean hasReplyUserId;
		private final boolean startsWithAtMention;
		private final boolean containsAtMentions;
		private final boolean multipleAtMentions;
		private final boolean startsWithRetweetIndicator;
		private final boolean containsRetweetIndicators;
		private final boolean multipleRetweetIndicators;

		public AggregatedTweetFeatures(TweetFeatures f) {
			hasRetweetUserId = f.hasRetweetUserId;
			hasReplyUserId = f.hasReplyUserId;
			if (hasRetweetUserId || hasReplyUserId) {
				startsWithAtMention = containsAtMentions = multipleAtMentions = false;
				startsWithRetweetIndicator = containsRetweetIndicators = multipleRetweetIndicators = false;
			} else {
				startsWithAtMention = f.startsWithAtMention;
				containsAtMentions = (f.numAtMentions > 0);
				multipleAtMentions = (f.numAtMentions > 1);
				startsWithRetweetIndicator = f.startsWithRetweetIndicator;
				containsRetweetIndicators = (f.numRetweetIndicators > 0);
				multipleRetweetIndicators = (f.numRetweetIndicators > 1);
			}
		}

		@Override
		public int hashCode() {
			int result = 31 + (hasReplyUserId ? 1231 : 1237);
			result = 31 * result + (hasRetweetUserId ? 1231 : 1237);
			result = 31 * result + (startsWithAtMention ? 1231 : 1237);
			result = 31 * result + (containsAtMentions ? 1231 : 1237);
			result = 31 * result + (multipleAtMentions ? 1231 : 1237);
			result = 31 * result + (startsWithRetweetIndicator ? 1231 : 1237);
			result = 31 * result + (containsRetweetIndicators ? 1231 : 1237);
			return 31 * result + (multipleRetweetIndicators ? 1231 : 1237);
		}

		@Override
		public boolean equals(Object obj) {
			if (!(obj instanceof AggregatedTweetFeatures))
				return false;
			AggregatedTweetFeatures other = (AggregatedTweetFeatures) obj;
			if ((hasReplyUserId == other.hasReplyUserId) && (hasRetweetUserId == other.hasRetweetUserId) &&
				(startsWithAtMention == other.startsWithAtMention) &&
				(containsAtMentions == other.containsAtMentions) && (multipleAtMentions == other.multipleAtMentions) &&
				(startsWithRetweetIndicator == other.startsWithRetweetIndicator) &&
				(containsRetweetIndicators == other.containsRetweetIndicators) &&
				(multipleRetweetIndicators == other.multipleRetweetIndicators))
				return true;
			return false;
		}

		@Override
		public String toString() {
			StringBuffer sb = new StringBuffer();
			if (hasReplyUserId)
				sb.append("reply via UI");
			else if (hasRetweetUserId)
				sb.append("retweet via UI");
			else {
				sb.append("Tweet");
				if (multipleAtMentions)
					sb.append(", multiple @-mentions");
				else if (containsAtMentions)
					sb.append(", single @-mention");
				if (startsWithAtMention)
					sb.append(" (at beginning)");

				if (multipleRetweetIndicators)
					sb.append(", multiple RT indicators");
				else if (containsRetweetIndicators)
					sb.append(", single RT indicator");
				if (startsWithRetweetIndicator)
					sb.append(" (at beginning)");
			}
			return sb.toString();
		}
	}

	private final Map<TweetFeatures, Integer> histo = new HashMap<TweetFeatures, Integer>();
	private final SummaryStatistics perUserAddressivity = new SummaryStatistics();
	private final SummaryStatistics perUserOriginality = new SummaryStatistics();
	private int numTweets = 0;

	private void processTweets(TwitterDao dao) throws Exception {
		Set<Long> userIds = dao.getUserIds(false);
		Map<Long, String> userNames = dao.getScreenNames(userIds);

		int idx = 0;
		for (Long userId : userIds) {
			String userName = userNames.get(userId);
			int numAddressiveTweets = 0;
			int numOriginalTweets = 0;
			int numTotalUserTweets = 0;

			List<SocialMediaMessage> tweets = dao.getMessages(userId, null, null);
			for (SocialMediaMessage tweet : tweets) {
				TweetFeatures f = TweetFeatures.classifyTweet(tweet, userId, userName);
				Integer cnt = histo.get(f);
				if (cnt == null)
					cnt = 0;
				histo.put(f, cnt + 1);

				if (f.isAddressive())
					numAddressiveTweets++;
				if (f.isOriginal())
					numOriginalTweets++;
				numTotalUserTweets++;
				numTweets++;
			}

			perUserAddressivity.addValue((double) numAddressiveTweets / numTotalUserTweets);
			perUserOriginality.addValue((double) numOriginalTweets / numTotalUserTweets);

			if ((++idx % 100) == 0)
				System.err.println("processed " + idx + " users");
		}
	}

	private void printReport() {
		System.out.println("overall number of tweets: " + numTweets);
		System.out.println("retweet indicators: 'RT' = " + TweetFeatures.numRetweetIndicatorRT +
				", 'via' = " + TweetFeatures.numRetweetIndicatorVia);
		System.out.println("found " + TweetFeatures.numRetweetsOfSelf + " retweets of own posts, " +
				TweetFeatures.numRepliesToSelf + " replies to own posts");
		System.out.println();

		System.out.println("per-user addressivity: mean = " + perUserAddressivity.getMean() +
				", std.dev. = " + perUserAddressivity.getStandardDeviation());
		System.out.println("per-user originality: mean = " + perUserOriginality.getMean() +
				", std.dev. = " + perUserOriginality.getStandardDeviation());
		System.out.println();

		for (Map.Entry<TweetFeatures, Integer> e : histo.entrySet())
			System.out.println(e.getValue() + ";" + e.getKey());
		System.out.println();

		Map<AggregatedTweetFeatures, Integer> aggHisto = new HashMap<AggregatedTweetFeatures, Integer>();
		for (Map.Entry<TweetFeatures, Integer> e : histo.entrySet()) {
			AggregatedTweetFeatures f = new AggregatedTweetFeatures(e.getKey());
			Integer cnt = aggHisto.get(f);
			if (cnt == null)
				cnt = 0;
			aggHisto.put(f, cnt + e.getValue());
		}
		for (Map.Entry<AggregatedTweetFeatures, Integer> e : aggHisto.entrySet())
			System.out.println(e.getValue() + ";" + e.getKey());
		System.out.println();
	}

	public static void main(String[] args) throws Exception {
		TwitterDao dao = new TwitterDao();
		TweetTaxonomy tax = new TweetTaxonomy();
		tax.processTweets(dao);
		tax.printReport();
	}

}
