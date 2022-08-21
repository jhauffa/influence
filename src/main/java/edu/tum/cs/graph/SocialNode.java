package edu.tum.cs.graph;

import java.io.Serializable;
import java.util.HashSet;
import java.util.Set;

import edu.tum.cs.time.TimeInterval;
import edu.tum.cs.util.io.SerializableLongIntHashMap;
import edu.tum.cs.util.io.SerializableLongObjectHashMap;

public class SocialNode implements Serializable {

	private static final long serialVersionUID = 1551215629731732514L;

	public final long id;
	public final boolean isCore;

	public SerializableLongObjectHashMap<double[]> thetaMap = new SerializableLongObjectHashMap<>();
	public SerializableLongObjectHashMap<double[]> exposureThetaMap = new SerializableLongObjectHashMap<>();
	public SerializableLongObjectHashMap<double[]> senderThetaMap = new SerializableLongObjectHashMap<>();
	public SerializableLongObjectHashMap<double[]> recipientThetaMap = new SerializableLongObjectHashMap<>();

	public SerializableLongIntHashMap numNonAddrMap = new SerializableLongIntHashMap();
	public SerializableLongIntHashMap numGotSharedMap = new SerializableLongIntHashMap();

	// only used temporarily by GenerateSocialEdgesInDB
	public transient SerializableLongObjectHashMap<Set<Long>> exposersMap = new SerializableLongObjectHashMap<>();

	public int reportedNumIncoming;
	public int reportedNumOutgoing;

	public double pageRank;
	public double betweennessCentrality;

	public SocialNode(long id, boolean isCore) {
		this.id = id;
		this.isCore = isCore;
	}

	@Override
	public int hashCode() {
		return Long.valueOf(id).hashCode();
	}

	@Override
	public boolean equals(Object o) {
		return (o instanceof SocialNode && ((SocialNode) o).id == id);
	}

	/** Topic distribution of the sent non-addressive messages */
	public double[] getTheta(TimeInterval iv) {
		return thetaMap.get(iv.key());
	}

	public void setTheta(TimeInterval iv, double[] theta) {
		thetaMap.put(iv.key(), theta);
	}

	/** Topic distribution of the received non-addressive messages */
	public double[] getExposureTheta(TimeInterval iv) {
		return exposureThetaMap.get(iv.key());
	}

	public void setExposureTheta(TimeInterval iv, double[] theta) {
		exposureThetaMap.put(iv.key(), theta);
	}

	/** Topic distribution of the sent addressive messages */
	public double[] getSenderTheta(TimeInterval iv) {
		return senderThetaMap.get(iv.key());
	}

	public void setSenderTheta(TimeInterval iv, double[] theta) {
		senderThetaMap.put(iv.key(), theta);
	}

	/** Topic distribution of the received addressive messages */
	public double[] getRecipientTheta(TimeInterval iv) {
		return recipientThetaMap.get(iv.key());
	}

	public void setRecipientTheta(TimeInterval iv, double[] theta) {
		recipientThetaMap.put(iv.key(), theta);
	}

	public int getNumNonAddr(TimeInterval iv) {
		return numNonAddrMap.get(iv.key());
	}

	public void incNumNonAddr(TimeInterval iv) {
		numNonAddrMap.addTo(iv.key(), 1);
	}

	public int getNumGotShared(TimeInterval iv) {
		return numGotSharedMap.get(iv.key());
	}

	public void setNumGotShared(TimeInterval iv, int gotShared) {
		numGotSharedMap.put(iv.key(), gotShared);
	}

	public Set<Long> getExposers(TimeInterval iv) {
		return exposersMap.get(iv.key());
	}

	public void addExposer(TimeInterval iv, long exposerId) {
		long key = iv.key();
		Set<Long> exposers = exposersMap.get(key);
		if (exposers == null) {
			exposers = new HashSet<Long>();
			exposersMap.put(key, exposers);
		}
		exposers.add(exposerId);
	}

}
