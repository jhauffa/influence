package edu.tum.cs.db;

import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import edu.tum.cs.graph.SocialEdge;
import edu.tum.cs.graph.SocialNode;
import edu.tum.cs.util.ExperimentConfiguration;

public class SocialGraphDao extends DatabaseAccessor {

	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(SocialGraphDao.class);
	private static final String targetTableSuffix = cfg.getLocalProperty("tableSuffix", "");
	private static final int BATCH_SIZE = 10000;

	private static final long globalNodeId = -1;

	public SocialGraphDao() throws SQLException {
		logger.info("table suffix is '" + targetTableSuffix + "'");
	}

	public void createTables() throws SQLException {
		Connection c = getConnection();
		try {
			Statement s = c.createStatement();
			s.executeUpdate("CREATE TABLE `SOCIAL_NODE" + targetTableSuffix + "` (" +
					"`ID` bigint(20) NOT NULL," +
					"`IS_CORE` bit(1) NOT NULL," +
					"`REPORTED_NUM_INCOMING` int(11) NOT NULL," +
					"`REPORTED_NUM_OUTGOING` int(11) NOT NULL," +
					"`NUM_NON_ADDR_MAP` blob," +
					"`NUM_GOT_SHARED_MAP` blob," +
					"`THETA_MAP` blob," +
					"`EXPOSURE_THETA_MAP` blob," +
					"`SENDER_THETA_MAP` blob," +
					"`RECIPIENT_THETA_MAP` blob," +
					"`PAGE_RANK` double DEFAULT NULL," +
					"`BETWEENNESS_CENTRALITY` double DEFAULT NULL," +
					"PRIMARY KEY (`ID`))");
			s.executeUpdate("CREATE TABLE `SOCIAL_EDGE" + targetTableSuffix + "` (" +
					"`FROM_USER_ID` bigint(20) NOT NULL," +
					"`TO_USER_ID` bigint(20) NOT NULL," +
					"`IS_EXPLICIT` bit(1) NOT NULL," +
					"`NUM_ADDR_MAP` blob," +
					"`NUM_SHARED_MAP` blob," +
					"`THETA_MAP` blob," +
					"PRIMARY KEY (`FROM_USER_ID`,`TO_USER_ID`)," +
					"KEY `TO_USER_ID_INDEX` (`TO_USER_ID`))");
		} finally {
			c.close();
		}
	}

	public Set<Long> getUserIds() throws SQLException {
		Set<Long> ids = new HashSet<Long>();
		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement("SELECT ID FROM SOCIAL_NODE" + targetTableSuffix);
			ResultSet rs = ps.executeQuery();
			while (rs.next())
				ids.add(rs.getLong(1));
		} finally {
			c.close();
		}
		ids.remove(globalNodeId);
		return ids;
	}

	/** @return ID of the node representing the global topic distributions */
	public long getGlobalNodeId() {
		return globalNodeId;
	}

	public boolean hasExplicitEdges() throws SQLException {
		int numExplicitEdges;
		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement("SELECT COUNT(*) FROM SOCIAL_EDGE" + targetTableSuffix +
					" WHERE IS_EXPLICIT = TRUE");
			ResultSet rs = ps.executeQuery();
			rs.next();
			numExplicitEdges = rs.getInt(1);
		} finally {
			c.close();
		}
		return (numExplicitEdges > 0);
	}

	public boolean isExplicitSymmetric() throws SQLException {
		boolean isSymmetric = true;
		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement("SELECT COUNT(*) FROM SOCIAL_EDGE" + targetTableSuffix +
					" AS FWD LEFT JOIN SOCIAL_EDGE" + targetTableSuffix +
					" AS REV ON (FWD.FROM_USER_ID = REV.TO_USER_ID) AND (FWD.TO_USER_ID = REV.FROM_USER_ID) " +
					" WHERE (REV.FROM_USER_ID IS NULL) AND (FWD.IS_EXPLICIT = TRUE)");
			ResultSet rs = ps.executeQuery();
			rs.next();
			if (rs.getInt(1) > 0)
				isSymmetric = false;
		} finally {
			c.close();
		}
		return isSymmetric;
	}

	public void saveSocialEdges(Collection<SocialEdge> edges) throws SQLException, IOException {
		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement("INSERT INTO SOCIAL_EDGE" + targetTableSuffix +
					" (FROM_USER_ID, TO_USER_ID, IS_EXPLICIT, NUM_ADDR_MAP, NUM_SHARED_MAP, THETA_MAP) " +
					"VALUES (?, ?, ?, ?, ?, ?)");
			int i = 0;
			for (SocialEdge edge : edges) {
				ps.setLong(1, edge.fromUserId);
				ps.setLong(2, edge.toUserId);
				ps.setBoolean(3, edge.isExplicit);

				ps.setBytes(4, serializeObject(edge.numAddrMap));
				ps.setBytes(5, serializeObject(edge.numSharedMap));
				ps.setBytes(6, serializeObject(edge.thetaMap));

				ps.addBatch();
				if ((++i % BATCH_SIZE) == 0)
					ps.executeBatch();
			}
			ps.executeBatch();
		} finally {
			c.close();
		}
	}

	public void saveSocialNodes(Collection<SocialNode> nodes) throws SQLException, IOException {
		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement("INSERT INTO SOCIAL_NODE" + targetTableSuffix + " (ID, " +
					"IS_CORE, REPORTED_NUM_INCOMING, REPORTED_NUM_OUTGOING, PAGE_RANK, BETWEENNESS_CENTRALITY, " +
					"NUM_NON_ADDR_MAP, NUM_GOT_SHARED_MAP, THETA_MAP, EXPOSURE_THETA_MAP, SENDER_THETA_MAP, " +
					"RECIPIENT_THETA_MAP) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
			int i = 0;
			for (SocialNode node : nodes) {
				ps.setLong(1, node.id);
				ps.setBoolean(2, node.isCore);
				ps.setInt(3, node.reportedNumIncoming);
				ps.setInt(4, node.reportedNumOutgoing);
				ps.setDouble(5, node.pageRank);
				ps.setDouble(6, node.betweennessCentrality);

				ps.setBytes(7, serializeObject(node.numNonAddrMap));
				ps.setBytes(8, serializeObject(node.numGotSharedMap));

				ps.setBytes(9, serializeObject(node.thetaMap));
				ps.setBytes(10, serializeObject(node.exposureThetaMap));
				ps.setBytes(11, serializeObject(node.senderThetaMap));
				ps.setBytes(12, serializeObject(node.recipientThetaMap));

				ps.addBatch();
				if ((++i % BATCH_SIZE) == 0)
					ps.executeBatch();
			}
			ps.executeBatch();
		} finally {
			c.close();
		}
	}

	public Map<Long, SocialNode> getSocialNodes(Collection<Long> userIds)
			throws SQLException, IOException, ClassNotFoundException {
		Map<Long, SocialNode> nodes = new HashMap<>();
		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement("SELECT IS_CORE, REPORTED_NUM_INCOMING, REPORTED_NUM_OUTGOING, " +
					"PAGE_RANK, BETWEENNESS_CENTRALITY, NUM_NON_ADDR_MAP, NUM_GOT_SHARED_MAP, THETA_MAP, " +
					"EXPOSURE_THETA_MAP, SENDER_THETA_MAP, RECIPIENT_THETA_MAP FROM SOCIAL_NODE" + targetTableSuffix +
					" WHERE ID = ?");
			for (Long userId : userIds) {
				ps.setLong(1, userId);

				ResultSet rs = ps.executeQuery();
				try {
					if (rs.next()) {
						boolean isCore = rs.getBoolean(1);
						SocialNode node = new SocialNode(userId, isCore);

						node.reportedNumIncoming = rs.getInt(2);
						node.reportedNumOutgoing = rs.getInt(3);
						node.pageRank = rs.getDouble(4);
						node.betweennessCentrality = rs.getDouble(5);

						node.numNonAddrMap = deserializeObject(rs.getBytes(6), node.numNonAddrMap);
						node.numGotSharedMap = deserializeObject(rs.getBytes(7), node.numGotSharedMap);

						node.thetaMap = deserializeObject(rs.getBytes(8), node.thetaMap);
						node.exposureThetaMap = deserializeObject(rs.getBytes(9), node.exposureThetaMap);
						node.senderThetaMap = deserializeObject(rs.getBytes(10), node.senderThetaMap);
						node.recipientThetaMap = deserializeObject(rs.getBytes(11), node.recipientThetaMap);

						nodes.put(userId, node);
					}
				} finally {
					rs.close();
				}
			}
		} finally {
			c.close();
		}
		return nodes;
	}

	public SocialNode getSocialNode(long userId) throws SQLException, IOException, ClassNotFoundException {
		Map<Long, SocialNode> result = getSocialNodes(Arrays.asList(userId));
		return result.get(userId);
	}

	public Map<Long, SocialEdge> getOutgoingSocialEdges(long userId)
			throws SQLException, IOException, ClassNotFoundException {
		Map<Long, SocialEdge> edges = new HashMap<Long, SocialEdge>();
		Connection c = getConnection();
		try {
			PreparedStatement ps = c.prepareStatement("SELECT FROM_USER_ID, TO_USER_ID, IS_EXPLICIT, NUM_ADDR_MAP, " +
					"NUM_SHARED_MAP, THETA_MAP FROM SOCIAL_EDGE" + targetTableSuffix + " WHERE FROM_USER_ID = ?");
			ps.setLong(1, userId);

			ResultSet rs = ps.executeQuery();
			while (rs.next()) {
				long fromUserId = rs.getLong(1);
				long toUserId = rs.getLong(2);

				SocialEdge edge = new SocialEdge(fromUserId, toUserId);
				edge.isExplicit = rs.getBoolean(3);
				edge.numAddrMap = deserializeObject(rs.getBytes(4), edge.numAddrMap);
				edge.numSharedMap = deserializeObject(rs.getBytes(5), edge.numSharedMap);
				edge.thetaMap = deserializeObject(rs.getBytes(6), edge.thetaMap);

				edges.put(edge.toUserId, edge);
			}
		} finally {
			c.close();
		}
		return edges;
	}

}
