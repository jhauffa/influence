package edu.tum.cs.nlp.corpus;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

public class ProcessedMessage extends ProcessedDocument {

	private static final long serialVersionUID = 3353762393501672309L;

	private final int senderId;
	private final int[] recipientIds;
	private int prevIndex, nextIndex;

	public ProcessedMessage(int index, int prevIndex, int nextIndex, int[] wordIds, int senderId, int[] recipientIds) {
		super(index, wordIds);
		this.prevIndex = prevIndex;
		this.nextIndex = nextIndex;
		this.senderId = senderId;
		this.recipientIds = recipientIds;
	}

	private ProcessedMessage(ProcessedDocument doc, int prevIndex, int nextIndex, int senderId, int[] recipientIds) {
		this(doc.getIndex(), prevIndex, nextIndex, doc.getWordIds(), senderId, recipientIds);
	}

	// TODO: when implementing a social topic model with multiple recipients, this method (and getSuccessor) should take
	//	the recipient index as a parameter
	public int getPredecessor() {
		return prevIndex;
	}

	public void setPredecessor(int prevIndex) {
		this.prevIndex = prevIndex;
	}

	public int getSuccessor() {
		return nextIndex;
	}

	public void setSuccessor(int nextIndex) {
		this.nextIndex = nextIndex;
	}

	public int getSenderId() {
		return senderId;
	}

	public int[] getRecipientIds() {
		return recipientIds;
	}

	@Override
	public void write(ObjectOutputStream os) throws IOException {
		os.writeInt(prevIndex);
		os.writeInt(nextIndex);
		os.writeInt(senderId);
		os.writeInt(recipientIds.length);
		for (int id : recipientIds)
			os.writeInt(id);
		super.write(os);
	}

	protected static class ReaderImpl extends ProcessedDocument.ReaderImpl {
		@Override
		public ProcessedMessage read(ObjectInputStream is) throws IOException {
			int prevIndex = is.readInt();
			int nextIndex = is.readInt();
			int senderId = is.readInt();
			int[] recipientIds = new int[is.readInt()];
			for (int i = 0; i < recipientIds.length; i++)
				recipientIds[i] = is.readInt();
			return new ProcessedMessage(super.read(is), prevIndex, nextIndex, senderId, recipientIds);
		}
	}

	public static class Reader extends ReaderImpl implements ProcessedDocumentReader<ProcessedMessage> {
	}

}
