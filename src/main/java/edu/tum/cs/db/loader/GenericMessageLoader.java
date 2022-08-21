package edu.tum.cs.db.loader;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import edu.tum.cs.db.GenericSocialMediaDao;
import edu.tum.cs.db.entities.SocialMediaMessage;
import edu.tum.cs.nlp.corpus.TokenizedMessage;
import edu.tum.cs.nlp.tokenizer.TwitterTokenizer;
import edu.tum.cs.util.ExperimentConfiguration;

public class GenericMessageLoader extends MessageLoader {

	private static final ExperimentConfiguration cfg = new ExperimentConfiguration(GenericMessageLoader.class);
	private static final int minRecipientGroupSize = cfg.getLocalIntProperty("minRecipientGroupSize");

	private final TwitterTokenizer tokenizer;
	private final Set<Long> coreUserIds;

	public GenericMessageLoader(GenericSocialMediaDao dao, Set<Long> userIds) throws Exception {
		super(dao, userIds);
		tokenizer = new TwitterTokenizer();
		coreUserIds = dao.getUserIds(true);
	}

	@Override
	public List<TokenizedMessage> loadMessages(long userId, Date startDate, Date endDate,
			List<SocialMediaMessage> origMessages) throws Exception {
		List<SocialMediaMessage> messages = getRawMessages(userId, startDate, endDate);

		Map<Set<Long>, List<SocialMediaMessage>> recipientsToMessages =
				new HashMap<Set<Long>, List<SocialMediaMessage>>();
		for (SocialMediaMessage message : messages) {
			// information sharing -> discard
			long otherUserId = message.getRepostOfUserId();
			if ((otherUserId > 0L) && (otherUserId != userId)) {
				if (callback != null)
					callback.foundSharedMessage(message, userId, otherUserId);
				continue;
			}

			// treat messages without explicit recipients or a sufficiently large group of recipients as non-addressive
			Set<Long> recipients = message.getRecipients();
			recipients.remove(userId);
			if (recipients.isEmpty() || (recipients.size() > minRecipientGroupSize)) {
				if ((callback == null) || callback.foundNonAddressiveMessage(message, userId, recipients))
					addMessage(recipientsToMessages, message, userId, new HashSet<Long>(Arrays.asList(userId)));
			} else if ((callback == null) || callback.foundAddressiveMessage(message, userId, recipients))
				addMessage(recipientsToMessages, message, userId, recipients);
		}

		List<TokenizedMessage> documents = new ArrayList<TokenizedMessage>();
		for (Entry<Set<Long>, List<SocialMediaMessage>> e : recipientsToMessages.entrySet()) {
			for (SocialMediaMessage message : e.getValue()) {
				List<String> words = tokenizer.tokenize(message.getText(), true);
				TokenizedMessage msg = new TokenizedMessage(words, userId, e.getKey(), message.getDate());
				if (msg.size() > 0) {
					documents.add(msg);
					if (origMessages != null)
						origMessages.add(message);
				}
			}
		}
		return documents;
	}

	private void addMessage(Map<Set<Long>, List<SocialMediaMessage>> recipientsToMessages, SocialMediaMessage message,
			long sender, Set<Long> recipients) {
		recipients.retainAll(userIds);
		if (recipients.isEmpty() || (!coreUserIds.contains(sender) && Collections.disjoint(coreUserIds, recipients)))
			return;

		List<SocialMediaMessage> messages = recipientsToMessages.get(recipients);
		if (messages == null) {
			messages = new ArrayList<SocialMediaMessage>();
			recipientsToMessages.put(recipients, messages);
		}
		messages.add(message);
	}

}
