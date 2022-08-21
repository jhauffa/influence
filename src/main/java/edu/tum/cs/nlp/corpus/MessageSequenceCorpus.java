package edu.tum.cs.nlp.corpus;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;

public class MessageSequenceCorpus<N extends Comparable<N>, R extends Message<N>> extends MessageCorpus<N, R> {

	private static final long serialVersionUID = 8140352517332987797L;

	public MessageSequenceCorpus(Index<String> bagOfWords, Index<N> bagOfPersons, Collection<R> rawMessages) {
		super(bagOfWords, bagOfPersons, rawMessages);
	}

	public MessageSequenceCorpus(MessageSequenceCorpus<N, R> other, boolean discardSingleMessages) {
		this(filterMessages(other.documents, discardSingleMessages));
	}

	protected MessageSequenceCorpus(List<ProcessedMessage> messages) {
		super(messages);
	}

	private static List<ProcessedMessage> filterMessages(List<ProcessedMessage> messages,
			boolean discardSingleMessages) {
		List<ProcessedMessage> filteredMessages = new ArrayList<ProcessedMessage>();
		int pos = 0;
		Map<Integer, Integer> posMap = new HashMap<Integer, Integer>();
		posMap.put(-1, -1);
		for (ProcessedMessage msg : messages) {
			if (!discardSingleMessages || (msg.getPredecessor() != -1) || (msg.getSuccessor() != -1)) {
				filteredMessages.add(new ProcessedMessage(pos, msg.getPredecessor(), msg.getSuccessor(),
						msg.getWordIds(), msg.getSenderId(), msg.getRecipientIds()));
				posMap.put(msg.getIndex(), pos);
				pos++;
			}
		}
		for (ProcessedMessage msg : filteredMessages) {
			msg.setPredecessor(posMap.get(msg.getPredecessor()));
			msg.setSuccessor(posMap.get(msg.getSuccessor()));
		}
		return filteredMessages;
	}

	@Override
	protected Comparator<R> buildMessageComparator() {
		// sort by sender and date in preparation for assignment of predecessor and slicing
		return new Comparator<R>() {
			@Override
			public int compare(R o1, R o2) {
				int result = o1.getSender().compareTo(o2.getSender());
				if (result == 0)
					result = o1.getDate().compareTo(o2.getDate());
				return result;
			}
		};
	}

	@Override
	protected ProcessedMessage processMessage(Message<N> rawMessage, int docIdx, ArrayList<R> sortedRawMessages,
			Index<String> bagOfWords, Index<N> bagOfPersons) {
		ProcessedMessage msg = super.processMessage(rawMessage, docIdx, sortedRawMessages, bagOfWords, bagOfPersons);

		// search for closest previous messages with same sender and first recipient in either direction, use as
		// predecessor and successor
		if (rawMessage.getRecipients().size() == 1) {
			N singleRecipient = rawMessage.getRecipients().iterator().next();

			int candPredecessor = docIdx - 1;
			while (candPredecessor >= 0) {
				Message<N> candMessage = sortedRawMessages.get(candPredecessor);
				if (!candMessage.getSender().equals(rawMessage.getSender()))
					break;
				if ((candMessage.getRecipients().size() == 1) &&
					candMessage.getRecipients().contains(singleRecipient)) {
					msg.setPredecessor(candPredecessor);
					break;
				}
				candPredecessor--;
			}

			int candSuccessor = docIdx + 1;
			while (candSuccessor < sortedRawMessages.size()) {
				Message<N> candMessage = sortedRawMessages.get(candSuccessor);
				if (!candMessage.getSender().equals(rawMessage.getSender()))
					break;
				if ((candMessage.getRecipients().size() == 1) &&
					candMessage.getRecipients().contains(singleRecipient)) {
					msg.setSuccessor(candSuccessor);
					break;
				}
				candSuccessor++;
			}
		}

		return msg;
	}

	public static class MessageSequence implements Iterator<ProcessedMessage> {
		private List<ProcessedMessage> documents;
		private int pos;

		private MessageSequence(List<ProcessedMessage> documents, int startPos) {
			this.documents = documents;
			pos = startPos;
		}

		@Override
		public boolean hasNext() {
			return (pos >= 0);
		}

		@Override
		public ProcessedMessage next() {
			ProcessedMessage msg = documents.get(pos);
			pos = msg.getSuccessor();
			return msg;
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
	}

	private class SequenceIterator implements Iterator<MessageSequence> {
		private int pos = -1;

		public SequenceIterator() {
			fetchNext();
		}

		private void fetchNext() {
			do {
				pos++;
				if (pos >= documents.size())
					pos = -1;
			} while ((pos >= 0) && (documents.get(pos).getPredecessor() >= 0));
		}

		@Override
		public boolean hasNext() {
			return (pos >= 0);
		}

		@Override
		public MessageSequence next() {
			int curPos = pos;
			fetchNext();
			return new MessageSequence(documents, curPos);
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
	}

	public Iterator<MessageSequence> sequenceIterator() {
		return new SequenceIterator();
	}

	public void shuffle() {
		// Need to temporarily store all message sequences to be shuffled, as the iterator operates directly on the
		// indices we are modifying.
		Queue<List<ProcessedMessage>> sequences = new LinkedList<List<ProcessedMessage>>();
		Iterator<MessageSequence> it = new SequenceIterator();
		while (it.hasNext()) {
			MessageSequence seq = it.next();
			ProcessedMessage msg = seq.next();
			if (seq.hasNext()) {	// skip sequences of length 1
				List<ProcessedMessage> messages = new ArrayList<ProcessedMessage>();
				messages.add(msg);
				while (seq.hasNext())
					messages.add(seq.next());
				sequences.add(messages);
			}
		}

		while (!sequences.isEmpty()) {
			List<ProcessedMessage> messages = sequences.poll();
			Collections.shuffle(messages);
			int prevIdx = -1, curIdx = messages.get(0).getIndex(), nextIdx;
			for (int i = 0; i < messages.size(); i++) {
				if ((i + 1) < messages.size())
					nextIdx = messages.get(i + 1).getIndex();
				else
					nextIdx = -1;
				messages.get(i).setPredecessor(prevIdx);
				messages.get(i).setSuccessor(nextIdx);
				prevIdx = curIdx;
				curIdx = nextIdx;
			}
		}
	}

}
