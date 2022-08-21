package edu.tum.cs.nlp.corpus;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public abstract class Corpus<P extends ProcessedDocument> implements Iterable<P>, Serializable {

	private static final long serialVersionUID = 2639473522439224146L;

	protected final List<P> documents;

	public Corpus(List<P> documents) {
		this.documents = documents;
	}

	protected Corpus() {
		this.documents = new ArrayList<P>();
	}

	protected void assignDocuments(List<P> documents) {
		this.documents.clear();
		this.documents.addAll(documents);
	}

	/** @return the number of documents in the corpus */
	public int size() {
		return documents.size();
	}

	@Override
	public Iterator<P> iterator() {
		return documents.iterator();
	}

	/**
	 * Required for parallel LDA. Returns a list of numSlices iterables over pairwise disjoint subsets of documents.
	 * Ranges may be of arbitrary size, even empty, but for optimal performance the ranges should be of equal size (in
	 * terms of word count).
	 */
	public abstract List<Corpus<P>> slice(int numSlices, boolean keepOriginal);

	public int countTokens() {
		int numTokens = 0;
		for (P doc : this)
			numTokens += doc.size();
		return numTokens;
	}

}
