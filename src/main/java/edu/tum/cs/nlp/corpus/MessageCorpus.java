package edu.tum.cs.nlp.corpus;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class MessageCorpus<N extends Comparable<N>, R extends Message<N>> extends Corpus<ProcessedMessage> {

	private static final long serialVersionUID = -173361159826330253L;

	public MessageCorpus(Index<String> bagOfWords, Index<N> bagOfPersons, Collection<R> rawMessages) {
		super();
		assignDocuments(processMessages(rawMessages, bagOfWords, bagOfPersons));
	}

	protected MessageCorpus(List<ProcessedMessage> messages) {
		super(messages);
	}

	protected Comparator<R> buildMessageComparator() {
		// sort by sender in preparation for slicing
		return new Comparator<R>() {
			@Override
			public int compare(R o1, R o2) {
				return o1.getSender().compareTo(o2.getSender());
			}
		};
	}

	protected ProcessedMessage processMessage(Message<N> rawMessage, int docIdx, ArrayList<R> sortedRawMessages,
			Index<String> bagOfWords, Index<N> bagOfPersons) {
		int[] wordIds = new int[rawMessage.size()];
		int wordIdx = 0;
		for (String word : rawMessage)
			wordIds[wordIdx++] = bagOfWords.add(word);
		int senderId = bagOfPersons.add(rawMessage.getSender());
		int[] recipientIds = new int[rawMessage.getRecipients().size()];
		int recipientIdx = 0;
		for (N recipient : rawMessage.getRecipients())
			recipientIds[recipientIdx++] = bagOfPersons.add(recipient);
		return new ProcessedMessage(docIdx, -1, -1, wordIds, senderId, recipientIds);
	}

	protected List<ProcessedMessage> processMessages(Collection<R> rawMessages, Index<String> bagOfWords,
			Index<N> bagOfPersons) {
		ArrayList<R> sortedRawMessages = new ArrayList<R>(rawMessages);
		Collections.sort(sortedRawMessages, buildMessageComparator());

		List<ProcessedMessage> messages = new ArrayList<ProcessedMessage>(sortedRawMessages.size());
		int docIdx = 0;
		for (Message<N> rawMessage : sortedRawMessages) {
			messages.add(processMessage(rawMessage, docIdx, sortedRawMessages, bagOfWords, bagOfPersons));
			docIdx++;
		}
		return messages;
	}

	/**
 	 * If multiple messages may contribute to a single document-topic distribution theta, as is the case for ART
 	 * (multiple messages per author-recipient pair), that set of messages must not be split across different slices.
 	 * Ideally, we would compute the transitive closure of messages with the same author and at least one shared
 	 * recipient, and treat it as an indivisible unit. Since this is computationally expensive and requires extra
 	 * storage space, we simply keep all messages of an author together.
	 */
	@Override
	public List<Corpus<ProcessedMessage>> slice(int numSlices, boolean keepOriginal) {
		List<List<ProcessedMessage>> sliceDocuments = new ArrayList<List<ProcessedMessage>>(numSlices);
		for (int i = 0; i < numSlices; i++)
			sliceDocuments.add(new ArrayList<ProcessedMessage>());

		int startIndex = 0, curIndex = 0;
		int prevAuthor = -1;
		while (curIndex < documents.size()) {
			int curAuthor = documents.get(curIndex).getSenderId();
			if (prevAuthor < 0)
				prevAuthor = curAuthor;
			if ((curAuthor != prevAuthor) || (curIndex == (documents.size() - 1))) {
				int minIdx = 0;
				int minSize = Integer.MAX_VALUE;
				for (int i = 0; i < numSlices; i++) {
					if (sliceDocuments.get(i).size() < minSize) {
						minSize = sliceDocuments.get(i).size();
						minIdx = i;
					}
				}
				if (curIndex == (documents.size() - 1))
					curIndex++;
				sliceDocuments.get(minIdx).addAll(documents.subList(startIndex, curIndex));
				startIndex = curIndex;
				prevAuthor = curAuthor;
			}
			curIndex++;
		}

		List<Corpus<ProcessedMessage>> slices = new ArrayList<Corpus<ProcessedMessage>>(numSlices);
		for (List<ProcessedMessage> docs : sliceDocuments)
			slices.add(new MessageCorpus<N, R>(docs));
		if (!keepOriginal)
			documents.clear();
		return slices;
	}

}
