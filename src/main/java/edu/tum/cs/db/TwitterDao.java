package edu.tum.cs.db;

import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

import edu.tum.cs.db.entities.MessageTimestamp;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.db.entities.SocialMediaUrl;
import edu.tum.cs.db.entities.SocialMediaUser;
import edu.tum.cs.util.ExperimentConfiguration;

public class TwitterDao extends DatabaseAccessor implements SocialMediaDao {

	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(TwitterDao.class);

	public TwitterDao() throws SQLException {
	}

	@Override
	public Set<Long> getUserIds(boolean coreOnly) throws IOException {
		// All users on the ID list are fully crawled. The Twitter dataset does not contain any messages from non-core
		// to core users, so non-core users do not need to be considered.
		return ExperimentConfiguration.loadUserIds(cfg.getLocalProperty("usersList"),
				cfg.getLocalIntProperty("numUsers"));
	}

	@Override
	public SocialMediaUser getUser(long id) throws SQLException, IOException, ClassNotFoundException {
		SocialMediaUser user = null;
		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement("SELECT TOTAL_FOLLOWER_COUNT, TOTAL_FRIEND_COUNT, " +
					"TOTAL_TWEET_COUNT, CREATED_AT, UNCOMPRESS(`FOLLOWERS_SER`), UNCOMPRESS(`FRIENDS_SER`) " +
					"FROM USER WHERE ID = ?");
			ps.setLong(1, id);

			ResultSet rs = ps.executeQuery();
			if (rs.next()) {
				user = new SocialMediaUser();
				user.setId(id);
				user.setNumIncomingEdges(rs.getInt(1));
				user.setNumOutgoingEdges(rs.getInt(2));
				user.setNumMessages(rs.getInt(3));
				user.setCreatedAt(new Date(rs.getTimestamp(4).getTime()));

				Collection<Long> followers = deserializeObject(rs.getBytes(5), Collections.<Long>emptyList());
				user.getIncomingEdges().addAll(followers);
				Collection<Long> friends = deserializeObject(rs.getBytes(6), Collections.<Long>emptyList());
				user.getOutgoingEdges().addAll(friends);
			}
		} finally {
			c.close();
		}
		return user;
	}

	@Override
	public BiMap<Long, String> getScreenNames(Set<Long> userIdSet) throws SQLException {
		BiMap<Long, String> screenNames = HashBiMap.<Long, String>create(userIdSet.size());
		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement("SELECT SCREEN_NAME FROM USER WHERE ID = ?");
			for (Long id : userIdSet) {
				ps.setLong(1, id);
				ResultSet rs = ps.executeQuery();
				if (rs.next()) {
					String name = rs.getString(1);
					try {
						screenNames.put(id, name);
					} catch (IllegalArgumentException ex) {
						Long oldId = screenNames.inverse().get(name);
						if (oldId != null) {
							logger.warning("screen name '" + name + "' already present, old id = " + oldId +
									", new id = " + id);
						} else
							logger.log(Level.WARNING, "error storing screen name", ex);
					}
				}
			}
		} finally {
			c.close();
		}
		return screenNames;
	}

	public List<SocialMediaMessage> getMessages(long userId, Date startDate, Date endDate, boolean preferOriginal)
			throws SQLException {
		String sql = "SELECT T1.ID, T1.STATUS_TEXT, T1.CREATED_AT, " +
			"T1.USER_ID, T1.IN_REPLY_TO_USER_ID, T1.RETWEET_OF_USER_ID";
		if (preferOriginal)
			sql += ", T2.STATUS_TEXT";
		sql += " FROM TWEET AS T1";
		if (preferOriginal)
			sql += " LEFT JOIN TWEET AS T2 ON T2.ID = T1.RETWEET_OF_STATUS_ID";
		sql += " WHERE T1.USER_ID = ?";
		if (startDate != null)
			sql += " AND T1.CREATED_AT >= ?";
		if (endDate != null)
			sql += " AND T1.CREATED_AT < ?";
		sql += " ORDER BY T1.CREATED_AT ASC";

		List<SocialMediaMessage> tweets = new ArrayList<SocialMediaMessage>();
		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement(sql);
			ps.setFetchSize(maxFetchSize);
			int idx = 1;
			ps.setLong(idx++, userId);
			if (startDate != null)
				ps.setTimestamp(idx++, new Timestamp(startDate.getTime()));
			if (endDate != null)
				ps.setTimestamp(idx++, new Timestamp(endDate.getTime()));

			ResultSet rs = ps.executeQuery();
			while (rs.next()) {
				SocialMediaMessage tweet = new SocialMediaMessage();
				tweet.setId(rs.getLong(1));
				if (preferOriginal)
					tweet.setText(rs.getString(7));
				if (tweet.getText() == null)
					tweet.setText(rs.getString(2));
				tweet.setDate(new Date(rs.getTimestamp(3).getTime()));
				tweet.setSenderId(rs.getLong(4));
				tweet.setInReplyToUserId(rs.getLong(5));
				tweet.setRepostOfUserId(rs.getLong(6));
				tweets.add(tweet);
			}
		} finally {
			c.close();
		}
		return tweets;
	}

	@Override
	public List<SocialMediaMessage> getMessages(long userId, Date startDate, Date endDate) throws SQLException {
		return getMessages(userId, startDate, endDate, false);
	}

	private static final long minSnowflakeId = 419430400000L;
	private static final long snowflakeOffset = 1288834974657L;
	private static final Date firstTweetDate = new Date(1142974214000L);
	private static final int reportInterval = 1000000;

	public static Date fixMessageDate(long messageId, Date messageDate, Date knownPriorDate) {
		if (messageId > minSnowflakeId)
			return new Date(snowflakeOffset + (messageId >> 22));
		if (messageDate.before(knownPriorDate))
			logger.warning(messageId + " has suspicious date " + messageDate);
		return messageDate;
	}

	private static final int maxValuesPerQuery = 10000;

	private static String generateArgumentList(int n) {
		StringBuilder sb = new StringBuilder();
		if (n > 0) {
			sb.append("?");
			for (int i = 1; i < n; i++)
				sb.append(",?");
		}
		return sb.toString();
	}

	@Override
	public void buildMessageChain(MessageTimestamp.MessageType type,
			Map<MessageTimestamp, Collection<MessageTimestamp>> chain, Set<MessageTimestamp> roots,
			Set<MessageTimestamp> incompleteRoots, int minSeqLength) throws SQLException {
		// Using *_USER_ID for filtering, because these columns have indices. However, due to a bug in Twitter's API,
		// there are tweets with IN_REPLY_TO_USER_ID > 0 and IN_REPLY_TO_STATUS_ID = -1 that need to be filtered out.
		String sqlLink;
		if (type == null) {
			sqlLink = "SELECT ID, CREATED_AT, IN_REPLY_TO_STATUS_ID, RETWEET_OF_STATUS_ID FROM TWEET WHERE " +
					"((IN_REPLY_TO_USER_ID > 0) AND (IN_REPLY_TO_STATUS_ID > 0)) OR (RETWEET_OF_USER_ID > 0)";
		} else if (type == MessageTimestamp.MessageType.REPLY) {
			// For unknown reasons, using the index in this particular query causes a massive slowdown.
			// TODO: It is conceivable that a linear table scan is actually faster, since > 50% of records match the
			//	selection criteria, and the difference in speed can be explained by MySQL falling back to a linear scan
			//	for the above query.
			sqlLink = "SELECT ID, CREATED_AT, IN_REPLY_TO_STATUS_ID FROM TWEET WHERE IN_REPLY_TO_STATUS_ID > 0";
		} else {
			sqlLink = "SELECT ID, CREATED_AT, RETWEET_OF_STATUS_ID FROM TWEET WHERE RETWEET_OF_USER_ID > 0";
		}
		String sqlRootPre = "SELECT ID, CREATED_AT FROM TWEET WHERE ID IN (";
		String sqlRootPost = ")";

		Map<MessageTimestamp, MessageTimestamp> entities = new HashMap<MessageTimestamp, MessageTimestamp>();
		Connection c = getConnection();
		try {
			// collect backlinks
			logger.info("collecting backlinks...");
			int n = 0;
			PreparedStatement psLink = c.prepareStatement(sqlLink);
			psLink.setFetchSize(maxFetchSize);
			ResultSet rs = psLink.executeQuery();
			try {
				while (rs.next()) {
					MessageTimestamp protoChild = new MessageTimestamp(rs.getLong(1));
					Date childDate = fixMessageDate(protoChild.getId(), new Date(rs.getTimestamp(2).getTime()),
							firstTweetDate);
					if (childDate != null) {
						long parentId;
						MessageTimestamp.MessageType curType = type;
						if (curType == null) {
							parentId = rs.getLong(3);
							curType = MessageTimestamp.MessageType.REPLY;
							if (parentId <= 0) {
								parentId = rs.getLong(4);
								curType = MessageTimestamp.MessageType.SHARE;
							}
						} else {
							parentId = rs.getLong(3);
							curType = type;
						}

						MessageTimestamp child = entities.get(protoChild);
						if (child == null) {
							child = protoChild;
							entities.put(child, child);
						}
						child.setDate(childDate);
						child.setType(curType);
						roots.remove(child);

						MessageTimestamp protoParent = new MessageTimestamp(parentId);
						MessageTimestamp parent = entities.get(protoParent);
						if (parent == null) {
							parent = protoParent;
							entities.put(parent, parent);
							parent.setType(MessageTimestamp.MessageType.OTHER);
							roots.add(parent);
						}

						Collection<MessageTimestamp> children = chain.get(parent);
						if (children == null) {
							// appears to be more memory efficient than a LinkedList
							children = new ArrayList<MessageTimestamp>(1);
							chain.put(parent, children);
						}
						children.add(child);
					}

					if ((++n % reportInterval) == 0)
						logger.info("processed " + n + " records");
				}
			} finally {
				rs.close();
			}

			// use minimum sequence length to prune candidate chains
			logger.info("pruning candidate chains...");
			if (minSeqLength > 1) {	// at this point, each root has at least one child
				Iterator<MessageTimestamp> it = roots.iterator();
				while (it.hasNext()) {
					MessageTimestamp root = it.next();
					List<MessageTimestamp> elem = new ArrayList<MessageTimestamp>();
					elem.add(root);
					int idx = 0;
					while (idx < elem.size()) {
						MessageTimestamp parent = elem.get(idx++);
						Collection<MessageTimestamp> children = chain.get(parent);
						if (children != null)
							elem.addAll(children);
					}
					if ((elem.size() - 1) < minSeqLength) {
						for (MessageTimestamp e : elem) {
							chain.remove(e);
							entities.remove(e);
						}
						it.remove();
					}
				}
			}

			// identify roots and retrieve their timestamps
			logger.info("checking potential roots...");
			n = 0;
			long rootIds[] = new long[roots.size()];
			int idx = 0;
			for (MessageTimestamp r : roots)
				rootIds[idx++] = r.getId();
			Arrays.sort(rootIds);
			Set<MessageTimestamp> rootsToKeep = new HashSet<MessageTimestamp>();
			PreparedStatement psRoot;
			int numRemainingRoots = roots.size();
			if (numRemainingRoots > maxValuesPerQuery)
				psRoot = c.prepareStatement(sqlRootPre + generateArgumentList(maxValuesPerQuery) +
						sqlRootPost);
			else
				psRoot = c.prepareStatement(sqlRootPre + generateArgumentList(numRemainingRoots) + sqlRootPost);
			int curIdx = 1;
			for (long id : rootIds) {
				psRoot.setLong(curIdx++, id);
				numRemainingRoots--;
				if ((curIdx > maxValuesPerQuery) || (numRemainingRoots == 0)) {
					rs = psRoot.executeQuery();
					try {
						while (rs.next()) {
							MessageTimestamp candRoot = entities.get(new MessageTimestamp(rs.getLong(1)));
							Date d = fixMessageDate(candRoot.getId(), new Date(rs.getTimestamp(2).getTime()),
									firstTweetDate);
							if (d != null) {
								candRoot.setDate(d);
								rootsToKeep.add(candRoot);
								roots.remove(candRoot);
							}

							if ((++n % reportInterval) == 0)
								logger.info("processed " + n + " records");
						}
					} finally {
						rs.close();
					}
					if (numRemainingRoots < maxValuesPerQuery) {
						psRoot.close();
						psRoot = c.prepareStatement(sqlRootPre + generateArgumentList(numRemainingRoots) + sqlRootPost);
					}
					curIdx = 1;
				}
			}
			entities = null;
			for (MessageTimestamp r : roots)
				incompleteRoots.addAll(chain.get(r));
			roots.clear();
			roots.addAll(rootsToKeep);
		} finally {
			c.close();
		}
		logger.info("done (" + roots.size() + " complete, " + incompleteRoots.size() + " incomplete)");
	}

	public List<SocialMediaUrl> getResolvedUrls(Collection<SocialMediaMessage> messages) throws SQLException {
		List<SocialMediaUrl> urls = new ArrayList<SocialMediaUrl>();
		if (!messages.isEmpty()) {
			Connection c = getConnection();
			try {
				// We're not expecting more than a few hundred messages per call, so splitting across multiple queries
				// is not necessary.
				StringBuilder sb = new StringBuilder("SELECT TWEET_ID, ORIGINAL_URL, RESOLVED_URL FROM WEBSITE " +
						"WHERE TWEET_ID IN (?");
				for (int i = 1; i < messages.size(); i++)
					sb.append(",?");
				sb.append(")");

				PreparedStatement ps = c.prepareStatement(sb.toString());
				int idx = 1;
				Map<Long, SocialMediaMessage> idToMessage = new HashMap<Long, SocialMediaMessage>();
				for (SocialMediaMessage msg : messages) {
					ps.setLong(idx++, msg.getId());
					idToMessage.put(msg.getId(), msg);
				}

				ResultSet rs = ps.executeQuery();
				while (rs.next()) {
					SocialMediaUrl url = new SocialMediaUrl();
					url.setMessage(idToMessage.get(rs.getLong(1)));
					url.setOriginalUrl(rs.getString(2));
					url.setResolvedUrl(rs.getString(3));
					urls.add(url);
				}
			} finally {
				c.close();
			}
		}
		return urls;
	}

	public List<String> getWebsiteContent(long msgId, boolean useFilteredWebsites) throws SQLException {
		List<String> content = new ArrayList<String>();
		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement("SELECT CONTENT FROM WEBSITE" +
					(useFilteredWebsites ? "_FILTERED" : "") + " WHERE TWEET_ID = ?");
			ps.setLong(1, msgId);
			ResultSet rs = ps.executeQuery();
			while (rs.next())
				content.add(rs.getString(1));
		} finally {
			c.close();
		}
		return content;
	}

}
