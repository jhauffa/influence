package edu.tum.cs.graph;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import com.carrotsearch.hppc.LongObjectScatterMap;
import com.carrotsearch.hppc.cursors.LongCursor;
import com.carrotsearch.hppc.cursors.LongObjectCursor;
import com.carrotsearch.hppc.LongHashSet;

import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.SocialMediaDaoFactory;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.db.entities.SocialMediaUser;
import edu.tum.cs.db.loader.MessageLoader;
import edu.tum.cs.db.loader.MessageLoaderFactory;
import edu.tum.cs.time.IntervalCalculator;
import edu.tum.cs.util.ExperimentConfiguration;

public class DumpNetlist {

	private static final Logger logger = Logger.getLogger(DumpNetlist.class.getName());
	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(DumpNetlist.class);

	private static enum GraphType { EXPLICIT_RAW, EXPLICIT, COMMUNICATION };
	private static enum CommunicationGraphType { ADDRESSIVE_SHARE, ADDRESSIVE, SHARE };

	private static class ClassifyEdges implements MessageLoader.MessageTypeCallback {
		public final Map<CommunicationGraphType, Set<Long>> incidentNodes =
				new HashMap<CommunicationGraphType, Set<Long>>();
		private final Set<Long> userIds;

		public ClassifyEdges(Set<Long> userIds) {
			this.userIds = userIds;
			reset();
		}

		public void reset() {
			for (CommunicationGraphType type : CommunicationGraphType.values())
				incidentNodes.put(type, new HashSet<Long>());
		}

		@Override
		public boolean foundNonAddressiveMessage(SocialMediaMessage message, long senderId,
				Collection<Long> exposedIds) {
			return false;
		}

		@Override
		public boolean foundAddressiveMessage(SocialMediaMessage message, long senderId,
				Collection<Long> recipientIds) {
			for (long recipientId : recipientIds) {
				if (userIds.contains(recipientId)) {
					incidentNodes.get(CommunicationGraphType.ADDRESSIVE).add(recipientId);
					incidentNodes.get(CommunicationGraphType.ADDRESSIVE_SHARE).add(recipientId);
				}
			}
			return false;
		}

		@Override
		public boolean foundSharedMessage(SocialMediaMessage message, long senderId, long originalSenderId) {
			if (userIds.contains(originalSenderId)) {
				incidentNodes.get(CommunicationGraphType.SHARE).add(originalSenderId);
				incidentNodes.get(CommunicationGraphType.ADDRESSIVE_SHARE).add(originalSenderId);
			}
			return false;
		}
	}

	private static void writeCommunicationGraphs(SocialMediaDao dao, Set<Long> userIds, Date startDate, Date endDate)
			throws Exception {
		Map<CommunicationGraphType, Writer> out = new HashMap<CommunicationGraphType, Writer>();
		try {
			// A communication graph cannot contain isolated nodes, so the netlist format is fine. Use "ncol" extension
			// so that python-igraph automatically selects the correct import routine.
			for (CommunicationGraphType type : CommunicationGraphType.values())
				out.put(type, new BufferedWriter(new FileWriter("graph-" + type.name() + ".ncol")));

			MessageLoader loader = MessageLoaderFactory.createMessageLoader(dao, userIds);
			ClassifyEdges classifier = new ClassifyEdges(userIds);
			loader.setCallback(classifier);
			int idx = 0;
			for (long userId : userIds) {
				loader.loadMessages(userId, startDate, endDate);
				for (Map.Entry<CommunicationGraphType, Set<Long>> e : classifier.incidentNodes.entrySet()) {
					Writer w = out.get(e.getKey());
					for (Long toUserId : e.getValue())
						w.write(userId + "\t" + toUserId + "\n");
				}
				classifier.reset();

				if ((++idx % 1000) == 0)
					logger.info("processed " + idx + " of " + userIds.size() + " users");
			}
		} finally {
			for (Writer w : out.values())
				w.close();
		}
	}

	private static void writeExplicitSocialGraphRaw(SocialMediaDao dao, Set<Long> userIds) throws Exception {
		Writer out = new BufferedWriter(new FileWriter("graph.raw"));
		try {
			int idx = 0;
			for (long userId : userIds) {
				SocialMediaUser user = dao.getUser(userId);
				if (user.getOutgoingEdges() != null) {
					user.getOutgoingEdges().retainAll(userIds);
					for (long otherUserId : user.getOutgoingEdges())
						out.write(userId + "\t" + otherUserId + "\n");
				}
				if (user.getIncomingEdges() != null) {
					user.getIncomingEdges().retainAll(userIds);
					for (long otherUserId : user.getIncomingEdges())
						out.write(otherUserId + "\t" + userId + "\n");
				}

				if ((++idx % 1000) == 0)
					logger.info("processed " + idx + " of " + userIds.size() + " users");
			}
		} finally {
			out.close();
		}
	}

	private static void writeExplicitSocialGraph(SocialMediaDao dao, Set<Long> userIds) throws Exception {
		// Cannot avoid storing the whole graph in memory, because inconsistencies that arise from the crawling process
		// have to be resolved. For example, user a may be recorded as follower, but not friend of b, while b is
		// bidirectionally connected to a.
		LongObjectScatterMap<LongHashSet> adjLists = new LongObjectScatterMap<LongHashSet>();
		for (Long userId : userIds)
			adjLists.put(userId, new LongHashSet());
		for (Long userId : userIds) {
			SocialMediaUser user = dao.getUser(userId);
			user.getOutgoingEdges().retainAll(userIds);
			for (Long otherUserId : user.getOutgoingEdges())
				adjLists.get(userId).add(otherUserId);
			user.getIncomingEdges().retainAll(userIds);
			for (Long otherUserId : user.getIncomingEdges())
				adjLists.get(otherUserId).add(userId);
		}

		// LGL is the most primitive storage format that explicitly supports isolated nodes.
		Writer out = new BufferedWriter(new FileWriter("graph.lgl"));
		try {
			for (LongObjectCursor<LongHashSet> csr1 : adjLists) {
				out.write("# " + csr1.key + "\n");
				for (LongCursor csr2 : csr1.value)
					out.write(csr2.value + "\n");
			}
		} finally {
			out.close();
		}
	}

	public static void main(String[] args) throws Exception {
		if (args.length == 0) {
			System.err.println("usage: " + DumpNetlist.class.getSimpleName() + " graph [core?] [start day] [end day]");
			System.err.println("graph types:");
			for (GraphType t : GraphType.values())
				System.err.print(t + " ");
			System.err.println();
			return;
		}

		GraphType graphType = GraphType.valueOf(args[0]);
		boolean coreOnly = false;
		if (args.length > 1)
			coreOnly = Boolean.parseBoolean(args[1]);
		logger.info("writing " + graphType + " graph" + (coreOnly ? " (core users only)" : ""));
		SocialMediaDao dao = SocialMediaDaoFactory.createDao();
		Set<Long> userIds = dao.getUserIds(coreOnly);
		if (graphType == GraphType.COMMUNICATION) {
			Date cfgDate = cfg.getDateProperty(ExperimentConfiguration.PROP_END_DATE);
			Date startDate, endDate;
			if (args.length > 3) {
				SimpleDateFormat df = new SimpleDateFormat("yyyy.MM.dd");
				startDate = df.parse(args[2]);
				endDate = df.parse(args[3]);
			} else {
				endDate = cfgDate;
				IntervalCalculator ic = new IntervalCalculator(endDate);
				startDate = ic.getStartDate();
			}
			logger.info("for days " + startDate + "-" + endDate);

			writeCommunicationGraphs(dao, userIds, startDate, endDate);
		} else if (graphType == GraphType.EXPLICIT_RAW) {
			writeExplicitSocialGraphRaw(dao, userIds);
			logger.warning("Generated a raw netlist: edges are in no particular order, and the list may contain " +
					"duplicates. Use 'sort' and 'uniq' to build a proper netlist.");
			logger.info("Use 'sort list1 list2 | uniq -d' to compute the intersection.");
		} else {
			writeExplicitSocialGraph(dao, userIds);
		}
		logger.info("done");
	}

}
