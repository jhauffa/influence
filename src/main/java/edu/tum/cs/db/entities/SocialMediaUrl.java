package edu.tum.cs.db.entities;

import java.io.Serializable;
import java.util.Date;

public class SocialMediaUrl implements Serializable {

	private static final long serialVersionUID = -2151984463798606493L;

	private String originalUrl, resolvedUrl;
	private long messageId;
	private Date messageDate;

	public String getOriginalUrl() {
		return originalUrl;
	}

	public void setOriginalUrl(String originalUrl) {
		this.originalUrl = originalUrl;
	}

	public String getResolvedUrl() {
		return resolvedUrl;
	}

	public void setResolvedUrl(String resolvedUrl) {
		this.resolvedUrl = resolvedUrl;
	}

	public long getMessageId() {
		return messageId;
	}

	public Date getMessageDate() {
		return messageDate;
	}

	public void setMessage(SocialMediaMessage message) {
		messageId = message.getId();
		messageDate = message.getDate();
	}

	@Override
	public String toString() {
		return resolvedUrl;
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof SocialMediaUrl))
			return false;
		return getResolvedUrl().equals(((SocialMediaUrl) o).getResolvedUrl());
	}

	@Override
	public int hashCode() {
		return resolvedUrl.hashCode();
	}

}
