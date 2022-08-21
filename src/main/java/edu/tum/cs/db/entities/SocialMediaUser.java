package edu.tum.cs.db.entities;

import java.io.Serializable;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;

public class SocialMediaUser implements Serializable {

	private static final long serialVersionUID = -3030811933512498507L;

	private long id;
	private Date createdAt;

	private int numIncomingEdges;
	private int numOutgoingEdges;
	private int numMessages;

	private final Set<Long> incomingEdges = new HashSet<Long>();
	private final Set<Long> outgoingEdges = new HashSet<Long>();

	public long getId() {
		return id;
	}

	public void setId(long userId) {
		this.id = userId;
	}

	public Date getCreatedAt() {
		return createdAt;
	}

	public void setCreatedAt(Date createdAt) {
		this.createdAt = createdAt;
	}

	public int getNumIncomingEdges() {
		return numIncomingEdges;
	}

	public void setNumIncomingEdges(int numIncomingEdges) {
		this.numIncomingEdges = numIncomingEdges;
	}

	public int getNumOutgoingEdges() {
		return numOutgoingEdges;
	}

	public void setNumOutgoingEdges(int numOutgoingEdges) {
		this.numOutgoingEdges = numOutgoingEdges;
	}

	public int getNumMessages() {
		return numMessages;
	}

	public void setNumMessages(int numMessages) {
		this.numMessages = numMessages;
	}

	public Set<Long> getIncomingEdges() {
		return incomingEdges;
	}

	public Set<Long> getOutgoingEdges() {
		return outgoingEdges;
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof SocialMediaUser))
			return false;
		return getId() == ((SocialMediaUser) o).getId();
	}

	@Override
	public int hashCode() {
		return Long.valueOf(id).hashCode();
	}

}
