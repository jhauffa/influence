package edu.tum.cs.db;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

import edu.tum.cs.db.entities.MessageTimestamp;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.db.entities.SocialMediaUser;

public class GenericSocialMediaDao extends DatabaseAccessor implements SocialMediaDao {

	private final String tablePrefix;
	private final boolean hasExplicitRelationships;

	public GenericSocialMediaDao(String tablePrefix) throws SQLException {
		this.tablePrefix = tablePrefix;
		hasExplicitRelationships = (getNumRelationships() > 0);
	}

	private int getNumRelationships() throws SQLException {
		int numRelationships;
		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement("SELECT COUNT(*) FROM " + tablePrefix + "_RELATIONSHIPS");
			ResultSet rs = ps.executeQuery();
			rs.next();
			numRelationships = rs.getInt(1);
		} finally {
			c.close();
		}
		return numRelationships;
	}

	@Override
	public Set<Long> getUserIds(boolean coreOnly) throws SQLException {
		Set<Long> userIds = new HashSet<Long>();

		String sql = "SELECT `id` FROM " + tablePrefix + "_USER";
		if (coreOnly)
			sql += " WHERE isFullyCrawled = TRUE";

		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement(sql);
			ResultSet rs = ps.executeQuery();
			while (rs.next())
				userIds.add(rs.getLong(1));
		} finally {
			c.close();
		}

		return userIds;
	}

	@Override
	public SocialMediaUser getUser(long userId) throws SQLException {
		SocialMediaUser user = new SocialMediaUser();
		user.setId(userId);
		// leaving date of account creation and number of sent messages unset, only used for building interaction chains

		if (hasExplicitRelationships) {
			// Facebook is the only medium currently covered by this DAO that has explicit relationships, so we assume
			// that all relationships are symmetric.
			Connection c = getConnection();
			try {
				PreparedStatement ps = c.prepareStatement("SELECT numRelationships FROM " + tablePrefix +
						"_USER WHERE `id` = ?");
				ps.setLong(1, userId);

				ResultSet rs = ps.executeQuery();
				if (!rs.next())
					return null;
				int numRelationships = rs.getInt(1);
				user.setNumIncomingEdges(numRelationships);
				user.setNumOutgoingEdges(numRelationships);

				ps = c.prepareStatement("SELECT `to` FROM " + tablePrefix + "_RELATIONSHIPS WHERE `from` = ?");
				ps.setLong(1, userId);

				rs = ps.executeQuery();
				while (rs.next()) {
					long otherUserId = rs.getLong(1);
					user.getIncomingEdges().add(otherUserId);
					user.getOutgoingEdges().add(otherUserId);
				}
			} finally {
				c.close();
			}
		}

		return user;
	}

	@Override
	public List<SocialMediaMessage> getMessages(long userId, Date startDate, Date endDate) throws SQLException {
		List<SocialMediaMessage> messages = new ArrayList<SocialMediaMessage>();

		String sql = "SELECT m.`id`, m.content, m.date, m.sender, p.sender, o.sender FROM " +
				tablePrefix + "_MESSAGE AS m " +
				"LEFT JOIN " + tablePrefix + "_MESSAGE AS p ON p.`id` = m.parent " +
				"LEFT JOIN " + tablePrefix + "_MESSAGE AS o ON o.`id` = m.origin WHERE m.sender = ?";
		if (startDate != null)
			sql += " AND m.date >= ?";
		if (endDate != null)
			sql += " AND m.date < ?";
		sql += " ORDER BY m.date ASC";

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
				SocialMediaMessage message = new SocialMediaMessage();
				message.setId(rs.getLong(1));
				message.setText(rs.getString(2));
				if (rs.wasNull())
					message.setText("");
				message.setDate(new Date(rs.getDate(3).getTime()));
				message.setSenderId(rs.getLong(4));
				message.setInReplyToUserId(rs.getLong(5));
				message.setRepostOfUserId(rs.getLong(6));
				messages.add(message);
			}
			rs.close();

			// fetch recipients
			ps = c.prepareStatement("SELECT recipient FROM " + tablePrefix + "_RECIPIENTS WHERE message = ?");
			ps.setFetchSize(maxFetchSize);
			for (SocialMediaMessage message : messages) {
				ps.setLong(1, message.getId());
				rs = ps.executeQuery();
				while (rs.next())
					message.getRecipients().add(rs.getLong(1));
				rs.close();
			}
		} finally {
			c.close();
		}

		return messages;
	}

	@Override
	public void buildMessageChain(MessageTimestamp.MessageType type,
			Map<MessageTimestamp, Collection<MessageTimestamp>> chain, Set<MessageTimestamp> roots,
			Set<MessageTimestamp> incompleteRoots, int minSeqLength) throws SQLException {
		Map<MessageTimestamp, MessageTimestamp> entities = new HashMap<MessageTimestamp, MessageTimestamp>();
		if ((type == null) || (type == MessageTimestamp.MessageType.REPLY))
			buildMessageChain(MessageTimestamp.MessageType.REPLY, chain, entities, roots, incompleteRoots);
		if ((type == null) || (type == MessageTimestamp.MessageType.SHARE))
			buildMessageChain(MessageTimestamp.MessageType.SHARE, chain, entities, roots, incompleteRoots);
	}

	private void buildMessageChain(MessageTimestamp.MessageType type,
			Map<MessageTimestamp, Collection<MessageTimestamp>> chain, Map<MessageTimestamp, MessageTimestamp> entities,
			Set<MessageTimestamp> roots, Set<MessageTimestamp> incompleteRoots) throws SQLException {
		String sql = "SELECT child.`id`, child.date, parent.`id`, parent.date FROM " +
			tablePrefix + "_MESSAGE AS child LEFT JOIN " + tablePrefix + "_MESSAGE parent ON ";
		if (type == MessageTimestamp.MessageType.REPLY)
			sql += "(parent.`id` = child.parent) WHERE (child.isReply = TRUE)";
		else
			sql += "(parent.`id` = child.origin) WHERE (child.isShared = TRUE)";

		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement(sql);
			ps.setFetchSize(maxFetchSize);
			ResultSet rs = ps.executeQuery();
			while (rs.next()) {
				MessageTimestamp proto = new MessageTimestamp(rs.getLong(1));
				MessageTimestamp child = entities.get(proto);
				if (child == null) {
					child = proto;
					entities.put(child, child);
				}
				child.setType(type);
				child.setDate(new Date(rs.getTimestamp(2).getTime()));
				roots.remove(child);

				proto = new MessageTimestamp(rs.getLong(3));
				if (proto.getId() > 0) {
					MessageTimestamp parent = entities.get(proto);
					if (parent == null) {
						parent = proto;
						entities.put(parent, parent);
						parent.setType(MessageTimestamp.MessageType.OTHER);
						parent.setDate(new Date(rs.getTimestamp(4).getTime()));
						roots.add(parent);
					}

					Collection<MessageTimestamp> children = chain.get(parent);
					if (children == null) {
						children = new LinkedList<MessageTimestamp>();
						chain.put(parent, children);
					}
					children.add(child);
				} else {
					incompleteRoots.add(child);
				}
			}
		} finally {
			c.close();
		}
	}

	@Override
	public BiMap<Long, String> getScreenNames(Set<Long> userIdSet) {
		return HashBiMap.<Long, String>create();
	}

}
