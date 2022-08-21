package edu.tum.cs.db.loader;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.Set;

import edu.tum.cs.db.TwitterDao;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.db.entities.SocialMediaUrl;
import edu.tum.cs.nlp.corpus.TokenizedMessage;
import edu.tum.cs.nlp.tokenizer.TwitterTokenizer;
import edu.tum.cs.nlp.tokenizer.TwitterTokenizerUrlAware;
import edu.tum.cs.util.ExperimentConfiguration;

public class TwitterMessageLoader extends MessageLoader {

	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(TwitterMessageLoader.class);

	private static final boolean tokenizeUrls = cfg.getBooleanProperty(ExperimentConfiguration.PROP_TOKENIZE_URLS);
	private static final boolean resolveUrls = cfg.getLocalBooleanProperty("resolveURLs");
	private static final boolean augmentTweets = cfg.getLocalBooleanProperty("augmentTweets");

	private final TwitterDao dao;
	private final TwitterTokenizer tokenizer;
	private final TwitterTokenizer webTokenizer;
	private Map<Long, List<SocialMediaUrl>> userId2Urls;

	public TwitterMessageLoader(TwitterDao dao, Set<Long> userIdSet) throws Exception {
		super(dao, userIdSet);
		this.dao = dao;
		if (resolveUrls || tokenizeUrls)
			tokenizer = new TwitterTokenizerUrlAware(tokenizeUrls);
		else
			tokenizer = new TwitterTokenizer();
		webTokenizer = augmentTweets ? new TwitterTokenizer() : null;
	}

	@Override
	protected List<SocialMediaMessage> loadRawMessages(long userId, Date startDate, Date endDate) throws Exception {
		return dao.getMessages(userId, startDate, endDate, (resolveUrls || tokenizeUrls));
	}

	@Override
	public void cacheMessages(Date startDate, Date endDate) {
		super.cacheMessages(startDate, endDate);
		if (!resolveUrls || (userId2Messages == null))
			return;

		userId2Urls = new HashMap<Long, List<SocialMediaUrl>>(userIds.size());
		logger.info("Loading URLs");
		try {
			int numUsers = 0;
			for (Long userId : userIds) {
				userId2Urls.put(userId, dao.getResolvedUrls(userId2Messages.get(userId)));
				if (++numUsers % 1000 == 0)
					logger.info("Processed " + numUsers + " users");
			}
			logger.info("Processed " + numUsers + " users, done");
		} catch (Exception ex) {
			logger.log(Level.WARNING, "Caching URLs failed, falling back to direct DB access", ex);
			userId2Urls = null;
		}
	}

	@Override
	public List<TokenizedMessage> loadMessages(long userId, Date startDate, Date endDate,
			List<SocialMediaMessage> origMessages) throws Exception {
		List<SocialMediaMessage> tweets = getRawMessages(userId, startDate, endDate);
		if (resolveUrls) {
			List<SocialMediaUrl> resolvedUrls;
			if (userId2Urls != null) {
				resolvedUrls = new ArrayList<SocialMediaUrl>();
				for (SocialMediaUrl url : userId2Urls.get(userId)) {
					if (url.getMessageDate().after(startDate) && !url.getMessageDate().after(endDate))
						resolvedUrls.add(url);
				}
			} else
				resolvedUrls = dao.getResolvedUrls(tweets);
			((TwitterTokenizerUrlAware) tokenizer).setResolvedUrls(resolvedUrls);
		}

		Map<Long, List<SocialMediaMessage>> recipientToTweets = new HashMap<Long, List<SocialMediaMessage>>();
		for (SocialMediaMessage tweet : tweets) {
			// reply to different user via UI -> only visible to recipient (and followers of both sender and recipient,
			// but we pretend this doesn't happen) -> addressive communication
			long otherUserId = tweet.getInReplyToUserId();
			if ((otherUserId > 0L) && (otherUserId != userId)) {
				if ((callback == null) || callback.foundAddressiveMessage(tweet, userId, Arrays.asList(otherUserId)))
					addToRecipientList(recipientToTweets, otherUserId, tweet);
			} else {
				// retweet via UI -> visible to all followers -> information sharing, discard
				otherUserId = tweet.getRepostOfUserId();
				if ((otherUserId > 0L) && (otherUserId != userId)) {
					if (callback != null)
						callback.foundSharedMessage(tweet, userId, otherUserId);
					continue;
				}

				// all other tweets -> non-addressive communication (includes tweets with @-mentions, which are visible
				// to all followers)
				if ((callback == null) || callback.foundNonAddressiveMessage(tweet, userId,
						Collections.<Long>emptySet()))
					addToRecipientList(recipientToTweets, userId, tweet);
			}
		}

		List<TokenizedMessage> documents = new ArrayList<TokenizedMessage>();
		for (Entry<Long, List<SocialMediaMessage>> e : recipientToTweets.entrySet()) {
			List<Long> recipientIds = new ArrayList<Long>(1);
			recipientIds.add(e.getKey());
			for (SocialMediaMessage tweet : e.getValue()) {
				long originalUserId = tweet.getRepostOfUserId();
				boolean keepUrls = (resolveUrls || tokenizeUrls) ||
						((originalUserId > 0L) && userIds.contains(originalUserId));

				List<String> words = tokenizer.tokenize(tweet.getText(), keepUrls);
				if (augmentTweets) {
					List<String> webContent = dao.getWebsiteContent(tweet.getId(), true);
					for (String content : webContent) {
						if ((content != null) && !content.isEmpty())
							words.addAll(webTokenizer.tokenize(content, false));
					}
				}

				TokenizedMessage msg = new TokenizedMessage(words, userId, recipientIds, tweet.getDate());
				if (msg.size() > 0) {
					documents.add(msg);
					if (origMessages != null)
						origMessages.add(tweet);
				}
			}
		}
		return documents;
	}

	private void addToRecipientList(Map<Long, List<SocialMediaMessage>> recipientToTweets, Long recipient,
			SocialMediaMessage tweet) {
		if (!userIds.contains(recipient))
			return;

		List<SocialMediaMessage> tweetsForRecipient = recipientToTweets.get(recipient);
		if (tweetsForRecipient == null) {
			tweetsForRecipient = new ArrayList<SocialMediaMessage>();
			recipientToTweets.put(recipient, tweetsForRecipient);
		}
		tweetsForRecipient.add(tweet);
	}

}
