package edu.tum.cs.graph;

import java.io.Serializable;

import edu.tum.cs.time.TimeInterval;
import edu.tum.cs.util.io.SerializableLongIntHashMap;
import edu.tum.cs.util.io.SerializableLongObjectHashMap;

/**
 * Represents an edge of the social network graph
 */
public class SocialEdge implements Serializable {

	private static final long serialVersionUID = -7866875481454206476L;

	public long fromUserId;
	public long toUserId;

	public SerializableLongObjectHashMap<double[]> thetaMap = new SerializableLongObjectHashMap<>();
	public SerializableLongIntHashMap numAddrMap = new SerializableLongIntHashMap();
	public SerializableLongIntHashMap numSharedMap = new SerializableLongIntHashMap();

	/** true if there is an explicit directed social edge from fromUser to toUser */
	public boolean isExplicit = false;

	public SocialEdge(long fromUserId, long toUserId) {
		this.fromUserId = fromUserId;
		this.toUserId = toUserId;
	}

	@Override
	public int hashCode() {
		return Long.valueOf(fromUserId).hashCode() + 31 * Long.valueOf(toUserId).hashCode();
	}

	@Override
	public boolean equals(Object o) {
		return ((o instanceof SocialEdge) &&
				((SocialEdge) o).fromUserId == fromUserId && ((SocialEdge) o).toUserId == toUserId);
	}

	public double[] getTheta(TimeInterval iv) {
		return thetaMap.get(iv.key());
	}

	public void setTheta(TimeInterval iv, double[] theta) {
		thetaMap.put(iv.key(), theta);
	}

	public boolean hasTheta() {
		return !thetaMap.isEmpty();
	}

	public int getNumAddr(TimeInterval iv) {
		return numAddrMap.get(iv.key());
	}

	public void incNumAddr(TimeInterval iv) {
		numAddrMap.addTo(iv.key(), 1);
	}

	public int getNumShared(TimeInterval iv) {
		return numSharedMap.get(iv.key());
	}

	public void incNumShared(TimeInterval iv) {
		numSharedMap.addTo(iv.key(), 1);
	}

}
