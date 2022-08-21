package edu.tum.cs.time.hmm;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import edu.tum.cs.db.GenericSocialMediaDao;
import edu.tum.cs.db.SocialMediaDao;
import edu.tum.cs.db.TwitterDao;
import edu.tum.cs.db.entities.MessageTimestamp;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.db.entities.SocialMediaUser;

public class BuildInteractionChains {

	private static enum InteractionType { PER_USER, REPLY, SHARE, BOTH };

	private static final double scaleDiv = 1000.0;	// use seconds as base unit for output
	private static final long oneYear = 365 * 24 * 60 * 60;
	private static final long fiveYears = 5 * oneYear;
	private static final int defaultMinSeqLength = 2;

	private static final DateFormat fmt = SimpleDateFormat.getDateInstance(SimpleDateFormat.SHORT, Locale.US);
	private static Date minDate = null, maxDate = null;
	private static int minSeqLength = defaultMinSeqLength;
	private static long numShareZeroes = 0, numInteractionsZeroes = 0;
	private static long numShareNoZeroes = 0, numInteractionsNoZeroes = 0;

	private static void printSequence(PrintWriter w, List<MessageTimestamp> seq, boolean isComplete,
			double warnMinDelta, boolean excludeZeroes) throws IOException {
		List<ObservationReal> deltas = new ArrayList<ObservationReal>(Math.max(0, seq.size() - 1));
		MessageTimestamp prevTs = null;
		int seqNumShare = 0, seqNumInteractions = 0;
		for (MessageTimestamp ts : seq) {
			Date curDate = ts.getDate();
			if (((minDate != null) && curDate.before(minDate)) || ((maxDate != null) && curDate.after(maxDate))) {
				System.err.println("warning: date " + fmt.format(curDate) + " (" + curDate.getTime() +  ") of event " +
						ts.getId() + " out of range");
			}
			if (prevTs != null) {
				double delta = (curDate.getTime() - prevTs.getDate().getTime()) / scaleDiv;
				if (!excludeZeroes || (delta > 0.0)) {
					if (delta >= warnMinDelta) {
						System.err.println("warning: parent = " + prevTs.getId() + ", child = " + ts.getId() +
								", delta = " + delta);
					}
					deltas.add(new ObservationReal(delta));

					if (ts.getType() == MessageTimestamp.MessageType.SHARE)
						seqNumShare++;
					seqNumInteractions++;
				}
			}
			prevTs = ts;
		}

		if (seqNumInteractions >= minSeqLength) {
			TimeSequence deltaSeq = new TimeSequence(deltas, !isComplete);
			w.println(deltaSeq);
			if (excludeZeroes) {
				numShareNoZeroes += seqNumShare;
				numInteractionsNoZeroes += seqNumInteractions;
			} else {
				numShareZeroes += seqNumShare;
				numInteractionsZeroes += seqNumInteractions;
			}
		}
	}

	private static void printTwitterUserSequences(PrintWriter wz, PrintWriter wnz, TwitterDao dao) throws Exception {
		Collection<Long> userIds = dao.getUserIds(true);
		for (Long userId : userIds) {
			SocialMediaUser user = dao.getUser(userId);
			// The database contains a number of tweets with obviously bogus dates (before the first official
			// tweet in 2006) - reject all tweets that claim to have been sent before the sender's account was
			// created.
			List<SocialMediaMessage> tweets = dao.getMessages(userId, user.getCreatedAt(), null);
			if (!tweets.isEmpty()) {
				List<MessageTimestamp> seq = new ArrayList<MessageTimestamp>(tweets.size() + 1);

				// Does the database contain all Tweets of the user at the time of crawling? Not 100% accurate,
				// as new tweets could have been sent between querying the total number of tweets and retrieving
				// the tweets.
				boolean isComplete = false;
				if (tweets.size() >= user.getNumMessages()) {
					seq.add(new MessageTimestamp(-1, user.getCreatedAt(), MessageTimestamp.MessageType.OTHER));
					isComplete = true;
				}

				for (SocialMediaMessage tweet : tweets)
					seq.add(new MessageTimestamp(tweet.getId(), tweet.getDate(), MessageTimestamp.MessageType.OTHER));
				printSequence(wz, seq, isComplete, fiveYears, false);
				printSequence(wnz, seq, isComplete, Double.MAX_VALUE, true);
			}
		}
	}

	private static boolean isLikelyBackdatedPost(Date d) {
		Calendar cal = new GregorianCalendar();
		cal.setTime(d);
		return ((cal.get(Calendar.SECOND) == 0) &&
				((cal.get(Calendar.MINUTE) == 0) || (cal.get(Calendar.MINUTE) == 30)));
	}

	private static void printSocialUserSequences(PrintWriter wz, PrintWriter wnz, GenericSocialMediaDao dao,
			boolean isFacebook) throws Exception {
		Collection<Long> userIds = dao.getUserIds(true);
		for (Long userId : userIds) {
			List<SocialMediaMessage> messages = dao.getMessages(userId, null, null);
			if (!messages.isEmpty()) {
				// Special case for Facebook: Users can backdate wall posts. Such posts will always appear to have been
				// posted exactly at a full or half hour. Treat the first message that does not meet these criteria as
				// the actual first message. We only chop off the beginning, and do not care about the initial HMM
				// probabilities anyway, so the impact of a false positive should be low.
				if (isFacebook) {
					int firstTrueMessageIdx = 0;
					for (SocialMediaMessage message : messages) {
						if (!isLikelyBackdatedPost(message.getDate()))
							break;
						firstTrueMessageIdx++;
					}
					messages = messages.subList(firstTrueMessageIdx, messages.size());
				}

				// The account creation date is unknown, so we treat all sequences as complete and ignore the initial
				// probabilities.
				List<MessageTimestamp> seq = new ArrayList<MessageTimestamp>(messages.size());
				for (SocialMediaMessage message : messages) {
					seq.add(new MessageTimestamp(message.getId(), message.getDate(),
							MessageTimestamp.MessageType.OTHER));
				}
				printSequence(wz, seq, true, fiveYears, false);
				printSequence(wnz, seq, true, Double.MAX_VALUE, true);
			}
		}
	}

	private static void followChain(PrintWriter wz, PrintWriter wnz,
			Map<MessageTimestamp, Collection<MessageTimestamp>> chain, Set<MessageTimestamp> roots, boolean isComplete,
			boolean isFacebook) throws IOException {
		for (MessageTimestamp root : roots) {
			List<MessageTimestamp> seq = new ArrayList<MessageTimestamp>();
			boolean isCurChainComplete = isComplete;
			Queue<MessageTimestamp> messages = new LinkedList<MessageTimestamp>();
			messages.add(root);
			while (!messages.isEmpty()) {
				MessageTimestamp parent = messages.poll();
				// On Facebook, the first post of a reply or sharing chain may be backdated and should be skipped.
				if (isFacebook && (parent == root) && isLikelyBackdatedPost(parent.getDate()))
					isCurChainComplete = false;
				else
					seq.add(parent);

				Collection<MessageTimestamp> childMessages = chain.get(parent);
				if (childMessages != null) {
					for (MessageTimestamp child : childMessages) {
						// Reject all child messages that claim to have been sent before their parents.
						if (!child.getDate().before(parent.getDate()))
							messages.add(child);
						else
							System.err.println("warning: " + child.getId() + " before parent " + parent.getId());
					}
				}
			}

			Collections.sort(seq, new Comparator<MessageTimestamp>() {
				@Override
				public int compare(MessageTimestamp o1, MessageTimestamp o2) {
					return Long.compare(o1.getDate().getTime(), o2.getDate().getTime());
				}
			});
			printSequence(wz, seq, isCurChainComplete, oneYear, false);
			printSequence(wnz, seq, isCurChainComplete, Double.MAX_VALUE, true);
		}
	}

	private static void printSocialEventSequences(PrintWriter wz, PrintWriter wnz, InteractionType type,
			SocialMediaDao dao, boolean isFacebook) throws Exception {
		Map<MessageTimestamp, Collection<MessageTimestamp>> chain =
				new HashMap<MessageTimestamp, Collection<MessageTimestamp>>();
		Set<MessageTimestamp> roots = new HashSet<MessageTimestamp>();
		Set<MessageTimestamp> incompleteRoots = new HashSet<MessageTimestamp>();
		MessageTimestamp.MessageType msgType = null;
		if (type == InteractionType.REPLY)
			msgType = MessageTimestamp.MessageType.REPLY;
		else if (type == InteractionType.SHARE)
			msgType = MessageTimestamp.MessageType.SHARE;
		dao.buildMessageChain(msgType, chain, roots, incompleteRoots, minSeqLength);
		followChain(wz, wnz, chain, roots, true, isFacebook);
		followChain(wz, wnz, chain, incompleteRoots, false, isFacebook);
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 3) {
			System.err.println("usage: " + BuildInteractionChains.class.getSimpleName() + " output source type" +
					" [min. length] [min. date] [max. date]");
			System.err.print("values for type = ");
			for (InteractionType v : InteractionType.values())
				System.err.print(v.toString() + " ");
			System.err.println("date format mm/dd/yyyy");
			System.err.println();
			return;
		}

		String source = args[1];
		InteractionType type = InteractionType.valueOf(args[2]);
		if (args.length > 3) {
			minSeqLength = Integer.parseInt(args[3]);
			if (args.length > 4) {
				minDate = fmt.parse(args[4]);
				if (args.length > 5)
					maxDate = fmt.parse(args[5]);
			}
		}

		System.err.println("extracting " + source + "-" + type + " sequences of length >= " + minSeqLength);
		long startTime = System.currentTimeMillis();
		PrintWriter wz = TimeSequence.prepareOutput(new File(args[0] + "-zeroes.txt.gz"));
		try {
			PrintWriter wnz = TimeSequence.prepareOutput(new File(args[0] + "-no-zeroes.txt.gz"));
			try {
				if (source.equalsIgnoreCase("twitter")) {
					TwitterDao dao = new TwitterDao();
					if (type == InteractionType.PER_USER)
						printTwitterUserSequences(wz, wnz, dao);
					else
						printSocialEventSequences(wz, wnz, type, dao, false);
				} else {
					GenericSocialMediaDao dao = new GenericSocialMediaDao(source);
					boolean isFacebook = (source.equalsIgnoreCase("fb") || source.equalsIgnoreCase("facebook"));
					if (type == InteractionType.PER_USER)
						printSocialUserSequences(wz, wnz, dao, isFacebook);
					else
						printSocialEventSequences(wz, wnz, type, dao, isFacebook);
				}
			} finally {
				wnz.close();
			}
		} finally {
			wz.close();
		}

		System.err.println("zeroes: #interactions = " + numInteractionsZeroes + ", #share = " + numShareZeroes + " (" +
				(((double) numShareZeroes / numInteractionsZeroes) * 100.0) + "%)");
		System.err.println("no zeroes: #interactions = " + numInteractionsNoZeroes + ", #share = " + numShareNoZeroes +
				" (" + (((double) numShareNoZeroes / numInteractionsNoZeroes) * 100.0) + "%)");
		System.err.println(System.currentTimeMillis() - startTime);
		System.err.println("done");
	}

}
