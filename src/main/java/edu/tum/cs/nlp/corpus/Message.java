package edu.tum.cs.nlp.corpus;

import java.util.Collection;
import java.util.Date;

public interface Message<T extends Comparable<T>> extends Document {

	public T getSender();
	public Collection<T> getRecipients();
	public Date getDate();

}
