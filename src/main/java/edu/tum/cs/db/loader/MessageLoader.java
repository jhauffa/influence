package edu.tum.cs.db.loader;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.nlp.corpus.TokenizedMessage;

public abstract class MessageLoader {

	public interface MessageTypeCallback {
		public boolean foundNonAddressiveMessage(SocialMediaMessage message, long senderId,
				Collection<Long> exposedIds);
		public boolean foundAddressiveMessage(SocialMediaMessage message, long senderId, Collection<Long> recipientIds);
		public boolean foundSharedMessage(SocialMediaMessage message, long senderId, long orginalSenderId);
	}

	protected static final Logger logger = Logger.getLogger(MessageLoader.class.getName());

	protected final SocialMediaDao dao;
	protected final Set<Long> userIds;

	protected MessageTypeCallback callback;
	protected Map<Long, List<SocialMediaMessage>> userId2Messages;

	protected MessageLoader(SocialMediaDao dao, Set<Long> userIds) {
		this.dao = dao;
		this.userIds = userIds;
	}

	public void setCallback(MessageTypeCallback callback) {
		this.callback = callback;
	}

	public abstract List<TokenizedMessage> loadMessages(long userId, Date startDate, Date endDate,
			List<SocialMediaMessage> origMessages) throws Exception;

	public List<TokenizedMessage> loadMessages(long userId, Date startDate, Date endDate) throws Exception {
		return loadMessages(userId, startDate, endDate, null);
	}

	protected List<SocialMediaMessage> loadRawMessages(long userId, Date startDate, Date endDate) throws Exception {
		return dao.getMessages(userId, startDate, endDate);
	}

	public void cacheMessages(Date startDate, Date endDate) {
		userId2Messages = new HashMap<Long, List<SocialMediaMessage>>(userIds.size());
		logger.info("Loading messages");
		try {
			int numUsers = 0;
			for (Long userId : userIds) {
				userId2Messages.put(userId, loadRawMessages(userId, startDate, endDate));
				if (++numUsers % 1000 == 0)
					logger.info("Processed " + numUsers + " users");
			}
			logger.info("Processed " + numUsers + " users, done");
		} catch (Exception ex) {
			logger.log(Level.WARNING, "Caching messages failed, falling back to direct DB access", ex);
			userId2Messages = null;
		}
	}

	/**
	 * @return the messages within > startDate and <= endDate
	 */
	public List<SocialMediaMessage> getRawMessages(long userId, Date startDate, Date endDate) throws Exception {
		if (userId2Messages != null) {
			List<SocialMediaMessage> messagesInTimeSlice = new ArrayList<SocialMediaMessage>();
			for (SocialMediaMessage message : userId2Messages.get(userId)) {
				if (message.getDate().after(startDate)) {
					if (message.getDate().after(endDate))
						break;	// messages are temporally ordered
					messagesInTimeSlice.add(message);
				}
			}
			return messagesInTimeSlice;
		}
		return loadRawMessages(userId, startDate, endDate);
	}

}
