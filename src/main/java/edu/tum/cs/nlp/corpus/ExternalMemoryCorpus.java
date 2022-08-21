package edu.tum.cs.nlp.corpus;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.concurrent.LinkedBlockingQueue;

public class ExternalMemoryCorpus<P extends ProcessedDocument> extends Corpus<P> {

	private static final long serialVersionUID = 912261828084994992L;

	private class BufferFiller extends Thread {
		private final LinkedBlockingQueue<P> buf;
		private int numRemaining;
		public Exception ex = null;

		public BufferFiller(LinkedBlockingQueue<P> buf, int numMessages) {
			this.buf = buf;
			numRemaining = numMessages;
		}

		@Override
		public void run() {
			try {
				ObjectInputStream ois = new ObjectInputStream(new FileInputStream(storageFile));
				try {
					while ((numRemaining > 0) && !interrupted()) {
						P m = reader.read(ois);
						numRemaining--;
						buf.put(m);
					}
				} finally {
					ois.close();
				}
			} catch (Exception ex) {
				this.ex = ex;
				try {
					buf.put(null);
				} catch (InterruptedException innerEx) {
					// nothing we can do
				}
			}
		}
	}

	private class BufferIterator implements Iterator<P> {
		private static final int bufSize = 1000;
		private final LinkedBlockingQueue<P> buf = new LinkedBlockingQueue<P>(bufSize);
		private final BufferFiller filler;
		private int pos = 0;

		public BufferIterator() {
			filler = new BufferFiller(buf, numDocuments);
			filler.start();
		}

		@Override
		public boolean hasNext() {
			return (pos < numDocuments);
		}

		@Override
		public P next() {
			if (pos++ >= numDocuments)
				throw new NoSuchElementException();
			P m = null;
			try {
				m = buf.take();
				if (m == null)
					throw new RuntimeException("error reading from document storage", filler.ex);
			} catch (InterruptedException ex) {
				// ignore
			}
			return m;
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
	}

	private final ProcessedDocumentReader<P> reader;
	private final int numDocuments;
	private final File storageFile;

	public ExternalMemoryCorpus(Corpus<P> corpus, ProcessedDocumentReader<P> reader) {
		super((List<P>) null);
		this.reader = reader;
		numDocuments = corpus.size();

		try {
			storageFile = File.createTempFile("doc-", ".bin");
			storageFile.deleteOnExit();
			ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(storageFile));
			try {
				for (P m : corpus)
					m.write(oos);
			} finally {
				oos.close();
			}
		} catch (IOException ex) {
			throw new RuntimeException("error creating document storage", ex);
		}
	}

	@Override
	public int size() {
		return numDocuments;
	}

	@Override
	public Iterator<P> iterator() {
		return new BufferIterator();
	}

	@Override
	public List<Corpus<P>> slice(int numSlices, boolean keepOriginal) {
		throw new UnsupportedOperationException("external memory corpus cannot be subdivided");
	}

}
