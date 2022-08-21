package edu.tum.cs.nlp.corpus;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

public class ProcessedDocument implements Serializable {

	private static final long serialVersionUID = -4453057434136658437L;

	private final int index;
	private final int[] wordIds;

	public ProcessedDocument(int index, int[] wordIds) {
		this.index = index;
		this.wordIds = wordIds;
	}

	/** @return the number of words in the document. */
	public int size() {
		return wordIds.length;
	}

	public int getIndex() {
		return index;
	}

	public int[] getWordIds() {
		return wordIds;
	}

	public void write(ObjectOutputStream os) throws IOException {
		os.writeInt(index);
		os.writeInt(wordIds.length);
		for (int id : wordIds)
			os.writeInt(id);
	}

	protected static class ReaderImpl {
		public ProcessedDocument read(ObjectInputStream is) throws IOException {
			int index = is.readInt();
			int[] wordIds = new int[is.readInt()];
			for (int i = 0; i < wordIds.length; i++)
				wordIds[i] = is.readInt();
			return new ProcessedDocument(index, wordIds);
		}
	}

	public static class Reader extends ReaderImpl implements ProcessedDocumentReader<ProcessedDocument> {
	}

}
