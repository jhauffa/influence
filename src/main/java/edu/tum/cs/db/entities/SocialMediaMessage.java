package edu.tum.cs.db.entities;

import java.io.Serializable;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;

public class SocialMediaMessage implements Serializable {

	private static final long serialVersionUID = 4505468243010927065L;

	private long id;
	private long senderId;
	private String text;
	private Date date;
	private long inReplyToUserId;
	private long repostOfUserId;

	private final Set<Long> recipients = new HashSet<Long>();

	public long getId() {
		return id;
	}

	public void setId(long id) {
		this.id = id;
	}

	public long getSenderId() {
		return senderId;
	}

	public void setSenderId(long senderId) {
		this.senderId = senderId;
	}

	public String getText() {
		return text;
	}

	public void setText(String text) {
		this.text = text;
	}

	public Date getDate() {
		return date;
	}

	public void setDate(Date date) {
		this.date = date;
	}

	public long getInReplyToUserId() {
		return inReplyToUserId;
	}

	public void setInReplyToUserId(long inReplyToUserId) {
		this.inReplyToUserId = inReplyToUserId;
	}

	public long getRepostOfUserId() {
		return repostOfUserId;
	}

	public void setRepostOfUserId(long repostOfUserId) {
		this.repostOfUserId = repostOfUserId;
	}

	public Set<Long> getRecipients() {
		return recipients;
	}

	@Override
	public String toString() {
		return getText();
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof SocialMediaMessage))
			return false;
		return getId() == ((SocialMediaMessage) o).getId();
	}

	@Override
	public int hashCode() {
		return Long.valueOf(id).hashCode();
	}

}
