package edu.tum.cs.nlp.corpus;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/** helper class for passing corpus slices through methods that operate on a single corpus */
public class SlicedCorpus<P extends ProcessedDocument> extends Corpus<P> {

	private static final long serialVersionUID = -2996546982297203264L;

	private final int numDocuments;
	private final Corpus<P> original;
	private final ProcessedDocumentReader<P> reader;
	private final boolean useExternalMemoryCorpus;
	private List<Corpus<P>> slices;

	public SlicedCorpus(Corpus<P> corpus, ProcessedDocumentReader<P> reader, int numSlices,
			boolean useExternalMemoryCorpus, boolean keepOriginal) {
		super(null);
		numDocuments = corpus.size();
		original = corpus;
		this.reader = reader;
		this.useExternalMemoryCorpus = useExternalMemoryCorpus;
		slices = buildSlices(corpus, numSlices, keepOriginal);
	}

	@Override
	public int size() {
		return numDocuments;
	}

	@Override
	public Iterator<P> iterator() {
		return new Iterator<P>() {
			private int idxDoc = 0;
			private final Iterator<Corpus<P>> itSlice = slices.iterator();
			private Corpus<P> curSlice = null;
			private Iterator<P> itDoc;

			@Override
			public boolean hasNext() {
				return (idxDoc < numDocuments);
			}

			@Override
			public P next() {
				while ((curSlice == null) || !itDoc.hasNext()) {
					curSlice = itSlice.next();
					itDoc = curSlice.iterator();
				}
				idxDoc++;
				return itDoc.next();
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	private List<Corpus<P>> buildSlices(Corpus<P> corpus, int numSlices, boolean keepOriginal) {
		List<Corpus<P>> slices = new ArrayList<Corpus<P>>(numSlices);
		for (Corpus<P> slice : corpus.slice(numSlices, keepOriginal)) {
			if (useExternalMemoryCorpus)
				slices.add(new ExternalMemoryCorpus<P>(slice, reader));
			else
				slices.add(slice);
		}
		return slices;
	}

	@Override
	public List<Corpus<P>> slice(int numSlices, boolean keepOriginal) {
		if ((slices == null) || (numSlices != slices.size())) {
			if (original.size() == 0)
				throw new UnsupportedOperationException("original corpus is empty");
			slices = buildSlices(original, numSlices, keepOriginal);
		}
		return slices;
	}

}
