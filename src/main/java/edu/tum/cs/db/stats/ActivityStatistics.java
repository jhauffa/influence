package edu.tum.cs.db.stats;

import java.util.Collection;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.SocialMediaDaoFactory;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.db.entities.SocialMediaUser;
import edu.tum.cs.db.loader.MessageLoader;
import edu.tum.cs.db.loader.MessageLoaderFactory;
import edu.tum.cs.time.IntervalCalculator;
import edu.tum.cs.util.ExperimentConfiguration;

public class ActivityStatistics {

	private static final Logger logger = Logger.getLogger(ActivityStatistics.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(ActivityStatistics.class);

	private static class CompactSocialEdge {
		public long fromUserId;
		public long toUserId;

		public CompactSocialEdge(long fromUserId, long toUserId) {
			this.fromUserId = fromUserId;
			this.toUserId = toUserId;
		}

		@Override
		public int hashCode() {
			return Long.valueOf(fromUserId).hashCode() + 31 * Long.valueOf(toUserId).hashCode();
		}

		@Override
		public boolean equals(Object o) {
			return ((o instanceof CompactSocialEdge) &&
					(((CompactSocialEdge) o).fromUserId == fromUserId) &&
					(((CompactSocialEdge) o).toUserId == toUserId));
		}
	}

	private static class ClassifyMessages implements MessageLoader.MessageTypeCallback {
		private final Set<Long> userIdSet;
		private final Date startDate, endDate;
		private final Set<CompactSocialEdge> explicitEdges;
		private final boolean includeSharing;

		private boolean isActive, isActiveObs;
		private Set<CompactSocialEdge> activeEdges, activeEdgesObs, activeEdgesNotExp, activeEdgesNotExpObs;

		public ClassifyMessages(Set<Long> userIdSet, Date startDate, Date endDate,
				Set<CompactSocialEdge> explicitEdges, boolean includeSharing) {
			this.userIdSet = userIdSet;
			this.startDate = startDate;
			this.endDate = endDate;
			this.explicitEdges = explicitEdges;
			this.includeSharing = includeSharing;
			nextUser();
		}

		public void nextUser() {
			isActive = isActiveObs = false;
			activeEdges = new HashSet<CompactSocialEdge>();
			activeEdgesObs = new HashSet<CompactSocialEdge>();
			activeEdgesNotExp = new HashSet<CompactSocialEdge>();
			activeEdgesNotExpObs = new HashSet<CompactSocialEdge>();
		}

		public int getNumActiveEdges() {
			return activeEdges.size();
		}

		public int getNumActiveEdgesObs() {
			return activeEdgesObs.size();
		}

		public int getNumActiveEdgesNotExp() {
			return activeEdgesNotExp.size();
		}

		public int getNumActiveEdgesNotExpObs() {
			return activeEdgesNotExpObs.size();
		}

		@Override
		public boolean foundNonAddressiveMessage(SocialMediaMessage message, long senderId,
				Collection<Long> exposedIds) {
			isActive = true;
			if (message.getDate().after(startDate) && message.getDate().before(endDate))
				isActiveObs = true;
			return false;
		}

		private void registerEdge(long senderId, long recipientId, boolean isInObs) {
			if (userIdSet.contains(recipientId)) {
				CompactSocialEdge edge = new CompactSocialEdge(senderId, recipientId);
				activeEdges.add(edge);
				if (isInObs)
					activeEdgesObs.add(edge);
				if (!explicitEdges.contains(edge)) {
					activeEdgesNotExp.add(edge);
					if (isInObs)
						activeEdgesNotExpObs.add(edge);
				}
			}
		}

		@Override
		public boolean foundAddressiveMessage(SocialMediaMessage message, long senderId,
				Collection<Long> recipientIds) {
			boolean isInObs = (message.getDate().after(startDate) && message.getDate().before(endDate));
			for (long recipientId : recipientIds)
				registerEdge(senderId, recipientId, isInObs);
			return false;
		}

		@Override
		public boolean foundSharedMessage(SocialMediaMessage message, long senderId, long originalSenderId) {
			if (includeSharing) {
				boolean isInObs = (message.getDate().after(startDate) && message.getDate().before(endDate));
				registerEdge(senderId, originalSenderId, isInObs);
			}
			return false;
		}
	}

	public static void main(String[] args) throws Exception {
		boolean includeSharing = false;
		if (args.length > 0)
			includeSharing = Boolean.parseBoolean(args[0]);

		SocialMediaDao dao = SocialMediaDaoFactory.createDao();
		Set<Long> userIdSet = dao.getUserIds(true);

		logger.info("loading explicit edges...");
		Set<CompactSocialEdge> explicitEdges = new HashSet<CompactSocialEdge>();
		int i = 0;
		for (long userId : userIdSet) {
			SocialMediaUser user;
			try {
				user = dao.getUser(userId);
			} catch (Exception ex) {
				logger.log(Level.SEVERE, "error loading user " + userId, ex);
				throw ex;
			}

			if (user.getOutgoingEdges() != null) {
				user.getOutgoingEdges().retainAll(userIdSet);
				for (long otherUserId : user.getOutgoingEdges())
					explicitEdges.add(new CompactSocialEdge(userId, otherUserId));
			}

			if (user.getIncomingEdges() != null) {
				user.getIncomingEdges().retainAll(userIdSet);
				for (long otherUserId : user.getIncomingEdges())
					explicitEdges.add(new CompactSocialEdge(otherUserId, userId));
			}

			if ((++i % 1000) == 0)
				logger.info("processed " + i + " of " + userIdSet.size() + " users");
		}
		System.out.println("# explicit nodes = " + userIdSet.size() + ", # explicit edges = " + explicitEdges.size());

		logger.info("processing messages...");
		MessageLoader loader = MessageLoaderFactory.createMessageLoader(dao, userIdSet);
		Date endDate = cfg.getDateProperty(ExperimentConfiguration.PROP_END_DATE);
		IntervalCalculator ic = new IntervalCalculator(endDate);
		Date startDate = ic.getStartDate();
		ClassifyMessages classifier = new ClassifyMessages(userIdSet, startDate, endDate, explicitEdges,
				includeSharing);
		loader.setCallback(classifier);

		int numActiveNodes = 0, numActiveNodesObs = 0;
		int numActiveEdges = 0, numActiveEdgesObs = 0;
		int numActiveEdgesNotExp = 0, numActiveEdgesNotExpObs = 0;

		i = 0;
		for (Long userId : userIdSet) {
			loader.loadMessages(userId, null, null);
			if (classifier.isActive) {
				numActiveNodes++;
				if (classifier.isActiveObs)
					numActiveNodesObs++;
			}
			numActiveEdges += classifier.getNumActiveEdges();
			numActiveEdgesObs += classifier.getNumActiveEdgesObs();
			numActiveEdgesNotExp += classifier.getNumActiveEdgesNotExp();
			numActiveEdgesNotExpObs += classifier.getNumActiveEdgesNotExpObs();

			classifier.nextUser();
			if ((++i % 1000) == 0)
				logger.info("processed " + i + " of " + userIdSet.size() + " users");
		}

		System.out.println("# active nodes = " + numActiveNodes + ", in observation period = " + numActiveNodesObs);
		System.out.println("# active edges = " + numActiveEdges + ", in observation period = " + numActiveEdgesObs);
		System.out.println("# active edges (not explicit) = " + numActiveEdgesNotExp + ", in observation period = " +
				numActiveEdgesNotExpObs);
		logger.info("done");
	}

}
