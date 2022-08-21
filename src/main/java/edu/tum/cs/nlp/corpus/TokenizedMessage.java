package edu.tum.cs.nlp.corpus;

import java.io.Serializable;
import java.util.Collection;
import java.util.Date;
import java.util.Iterator;

public class TokenizedMessage implements Message<Long>, Serializable {

	private static final long serialVersionUID = 7166647663625861441L;

	private final Collection<String> words;
	private final Long sender;
	private final Collection<Long> recipients;
	private final Date date;

	public TokenizedMessage(Collection<String> words, Long sender, Collection<Long> recipients, Date date) {
		this.words = words;
		this.sender = sender;
		this.recipients = recipients;
		this.date = date;
	}

	@Override
	public Iterator<String> iterator() {
		return words.iterator();
	}

	@Override
	public int size() {
		return words.size();
	}

	@Override
	public Long getSender() {
		return sender;
	}

	@Override
	public Collection<Long> getRecipients() {
		return recipients;
	}

	@Override
	public Date getDate() {
		return date;
	}

}
