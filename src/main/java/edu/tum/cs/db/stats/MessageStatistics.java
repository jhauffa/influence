package edu.tum.cs.db.stats;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.SocialMediaDaoFactory;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.db.loader.MessageLoader;
import edu.tum.cs.db.loader.MessageLoaderFactory;

public class MessageStatistics {

	private static final Logger logger = Logger.getLogger(MessageStatistics.class.getName());

	private static class ClassifyMessages implements MessageLoader.MessageTypeCallback {
		private final Set<Long> userIdSet;

		public long numMessages, numAddressive, numNonAddressive, numSharedInside, numSharedAll;
		public SummaryStatistics perMessageRecipients = new SummaryStatistics();

		// messages may be counted both as addressive and non-addressive
		public HashSet<Long> seenMessages = new HashSet<Long>();

		public long numUserMessages, numUserAddressive, numUserNonAddressive, numUserSharedInside, numUserSharedAll;

		public ClassifyMessages(Set<Long> userIdSet) {
			this.userIdSet = userIdSet;
		}

		public void nextUser() {
			seenMessages.clear();
			numUserMessages = numUserAddressive = numUserNonAddressive = numUserSharedInside = numUserSharedAll = 0;
		}

		private void countMessage(SocialMediaMessage message) {
			if (seenMessages.add(message.getId())) {
				numMessages++;
				numUserMessages++;
			}
		}

		@Override
		public boolean foundNonAddressiveMessage(SocialMediaMessage message, long senderId,
				Collection<Long> exposedIds) {
			numNonAddressive++;
			numUserNonAddressive++;
			countMessage(message);
			return false;
		}

		@Override
		public boolean foundAddressiveMessage(SocialMediaMessage message, long senderId,
				Collection<Long> recipientIds) {
			int numCoreRecipients = 0;
			for (Long id : recipientIds) {
				if (userIdSet.contains(id))
					numCoreRecipients++;
			}
			if (numCoreRecipients > 0) {
				perMessageRecipients.addValue(numCoreRecipients);
				numAddressive++;
				numUserAddressive++;
				countMessage(message);
			}
			return false;
		}

		@Override
		public boolean foundSharedMessage(SocialMediaMessage message, long senderId, long originalSenderId) {
			if (userIdSet.contains(originalSenderId)) {
				numSharedInside++;
				numUserSharedInside++;
			}
			numSharedAll++;
			numUserSharedAll++;
			numNonAddressive++;
			numUserNonAddressive++;
			countMessage(message);
			return false;
		}
	}

	public static void main(String[] args) throws Exception {
		SocialMediaDao dao = SocialMediaDaoFactory.createDao();
		Set<Long> userIdSet = dao.getUserIds(true);
		MessageLoader loader = MessageLoaderFactory.createMessageLoader(dao, userIdSet);
		ClassifyMessages classifier = new ClassifyMessages(userIdSet);
		loader.setCallback(classifier);

		SummaryStatistics perUserAddressive = new SummaryStatistics();
		SummaryStatistics perUserNonAddressive = new SummaryStatistics();
		SummaryStatistics perUserSharedInside = new SummaryStatistics();
		SummaryStatistics perUserSharedAll = new SummaryStatistics();

		int i = 0;
		for (Long userId : userIdSet) {
			loader.loadMessages(userId, null, null);
			if (classifier.numUserMessages > 0) {
				perUserAddressive.addValue((double) classifier.numUserAddressive / classifier.numUserMessages);
				perUserNonAddressive.addValue((double) classifier.numUserNonAddressive / classifier.numUserMessages);
				perUserSharedInside.addValue((double) classifier.numUserSharedInside / classifier.numUserMessages);
				perUserSharedAll.addValue((double) classifier.numUserSharedAll / classifier.numUserMessages);
			}

			classifier.nextUser();
			if ((++i % 1000) == 0)
				logger.info("Done with " + i + " of " + userIdSet.size());
		}

		System.out.println("# messages = " + classifier.numMessages);
		System.out.println("% addressive = " + ((double) classifier.numAddressive / classifier.numMessages));
		System.out.println("% non-addressive = " + ((double) classifier.numNonAddressive / classifier.numMessages));
		System.out.println("% shared (inside) = " + ((double) classifier.numSharedInside / classifier.numMessages));
		System.out.println("% shared (all) = " + ((double) classifier.numSharedAll / classifier.numMessages));
		System.out.println("# recipients per message: mean = " + classifier.perMessageRecipients.getMean() +
				", std.dev. = " + classifier.perMessageRecipients.getStandardDeviation());
		System.out.println("% addressive per user: mean = " + perUserAddressive.getMean() +
				", std.dev. = " + perUserAddressive.getStandardDeviation());
		System.out.println("% non-addressive per user: mean = " + perUserNonAddressive.getMean() +
				", std.dev. = " + perUserNonAddressive.getStandardDeviation());
		System.out.println("% shared (inside) per user: mean = " + perUserSharedInside.getMean() +
				", std.dev. = " + perUserSharedInside.getStandardDeviation());
		System.out.println("% shared (all) per user: mean = " + perUserSharedAll.getMean() +
				", std.dev. = " + perUserSharedAll.getStandardDeviation());
		logger.info("Done");
	}

}
