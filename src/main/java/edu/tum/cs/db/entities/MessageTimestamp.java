package edu.tum.cs.db.entities;

import java.io.Serializable;
import java.util.Date;

public class MessageTimestamp implements Serializable {

	private static final long serialVersionUID = 3343255655235420240L;

	public static enum MessageType { REPLY, SHARE, OTHER };

	private long id;
	private Date date;
	private MessageType type;

	public MessageTimestamp() {
	}

	public MessageTimestamp(long id) {
		this(id, null, null);
	}

	public MessageTimestamp(long id, Date date, MessageType type) {
		this.id = id;
		this.date = date;
		this.type = type;
	}

	public long getId() {
		return id;
	}

	public void setId(long statusId) {
		this.id = statusId;
	}

	public Date getDate() {
		return date;
	}

	public void setDate(Date date) {
		this.date = date;
	}

	public MessageType getType() {
		return type;
	}

	public void setType(MessageType type) {
		this.type = type;
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof MessageTimestamp))
			return false;
		return id == ((MessageTimestamp) o).id;
	}

	@Override
	public int hashCode() {
		return Long.valueOf(id).hashCode();
	}

}
