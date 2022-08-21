package edu.tum.cs.nlp.corpus;

public interface Document extends Iterable<String> {

	/** @return the number of words in the document */
	public int size();

}
