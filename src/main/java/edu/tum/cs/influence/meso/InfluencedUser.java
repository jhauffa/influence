package edu.tum.cs.influence.meso;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import edu.tum.cs.graph.SocialEdge;
import edu.tum.cs.graph.SocialNode;
import edu.tum.cs.time.TimeInterval;

public abstract class InfluencedUser implements Cloneable {

	protected static final Logger logger = Logger.getLogger(InfluencedUser.class.getName());

	public static final double EPSILON = 0.001;

	private final SocialNode fromUserNode;
	private final long toUserId;

	protected final SocialGraphCache socialGraph;
	protected final Date centerDate;
	protected final int timePeriod;
	protected final boolean reverseTime;
	protected final int longestTimePeriod;
	protected final int numTopics;
	protected final double[] priorObserve;
	protected final double[] priorPredict;

	protected double[] thetaWritingUserNode;
	protected double[] thetaWritingUserEdge;
	protected double[] thetaWritingUserExposure;
	protected double[] thetaWritingUserIncoming;
	protected double[] thetaIncomingEdges;
	protected double[] thetaOutgoingEdges;
	protected double[] thetaOtherUsersNode;
	protected double[] thetaOtherUsersEdge;
	protected double[] thetaGlobalNode;
	protected double[] thetaGlobalEdge;

	protected double[] thetaToPredict;

	public InfluencedUser(SocialNode fromUserNode, long toUserId, SocialGraphCache socialGraph, Date centerDate,
			int timePeriod, int longestTimePeriod, boolean reverseTime, double[] priorObserve,
			double[] priorPredict) {
		this.fromUserNode = fromUserNode;
		this.toUserId = toUserId;
		this.socialGraph = socialGraph;
		this.centerDate = centerDate;
		this.timePeriod = timePeriod;
		this.longestTimePeriod = longestTimePeriod;
		this.reverseTime = reverseTime;
		this.numTopics = priorObserve.length;
		this.priorObserve = priorObserve;
		this.priorPredict = priorPredict;
	}

	@Override
	public InfluencedUser clone() {
		try {
			// a shallow copy is sufficient
			return (InfluencedUser) super.clone();
		} catch (CloneNotSupportedException ex) {
			throw new AssertionError(ex);
		}
	}

	public abstract void setThetas(SocialNode fromUserNode, SocialNode toUserNode, SocialEdge outgoingEdge,
			SocialEdge incomingEdge);
	public abstract double[][] getThetaMatrix();
	public abstract double[] getBaselineTheta();
	public abstract int getNumVariables();

	public void computeNeighbourhoodThetas(ModelVariant variant) {
		TimeInterval iv, ivExt;
		double[] priorDist;
		if (!reverseTime) {
			iv = new TimeInterval(centerDate, -timePeriod);
			ivExt = new TimeInterval(centerDate, -longestTimePeriod);
			priorDist = priorObserve;
		} else {
			iv = new TimeInterval(centerDate, timePeriod);
			ivExt = new TimeInterval(centerDate, longestTimePeriod);
			priorDist = priorPredict;
		}

		Map<Long, SocialEdge> outgoingEdges = socialGraph.getOutgoingEdges(fromUserNode.id);
		Set<Long> explicitOutDistTwo = null;
		if (variant.getIndicator().needsExplicitDistanceTwo())
			explicitOutDistTwo = socialGraph.getExplicitOutgoingDistanceTwo(outgoingEdges);
		Set<Long> commDistTwo = null;
		if (variant.getIndicator().needsCommunicationDistanceTwo())
			commDistTwo = socialGraph.getCommunicationDistanceTwo(outgoingEdges, iv);

		List<Double> weights = new ArrayList<Double>();
		double weightSum[] = new double[4];
		double[][] neighbourhoodTheta = new double[4][numTopics];
		for (SocialNode candidateNode : socialGraph.getNodes()) {
			if ((candidateNode.id == fromUserNode.id) || (candidateNode.id == toUserId))
				continue;

			// candidate --> origin
			SocialEdge incomingEdge = socialGraph.getOutgoingEdges(candidateNode.id).get(fromUserNode.id);
			// origin --> candidate
			SocialEdge outgoingEdge = outgoingEdges.get(candidateNode.id);

			int indicatorValue = variant.getIndicator().eval(fromUserNode, candidateNode, incomingEdge, outgoingEdge,
					socialGraph, iv, explicitOutDistTwo, commDistTwo);
			if (indicatorValue == 0)
				continue;

			double weightResult = variant.getWeight().eval(fromUserNode, candidateNode, incomingEdge, outgoingEdge,
					socialGraph, iv);
			weights.add(weightResult);
			if (weightResult == 0.0)
				continue;

			if (incomingEdge != null) {
				double[] theta = incomingEdge.getTheta(iv);
				if ((theta == null) && (timePeriod != longestTimePeriod))
					theta = incomingEdge.getTheta(ivExt);
				if (theta != null) {
					for (int i = 0; i < theta.length; i++)
						neighbourhoodTheta[0][i] += theta[i] * weightResult;
					weightSum[0] += weightResult;
				}
			}

			if (outgoingEdge != null) {
				double[] theta = outgoingEdge.getTheta(iv);
				if ((theta == null) && (timePeriod != longestTimePeriod))
					theta = outgoingEdge.getTheta(ivExt);
				if (theta != null) {
					for (int i = 0; i < theta.length; i++)
						neighbourhoodTheta[1][i] += theta[i] * weightResult;
					weightSum[1] += weightResult;
				}
			}

			double[] thetaNonAddr = candidateNode.getTheta(iv);
			if ((thetaNonAddr == null) && (timePeriod != longestTimePeriod))
				thetaNonAddr = candidateNode.getTheta(ivExt);
			if (thetaNonAddr != null) {
				for (int i = 0; i < thetaNonAddr.length; i++)
					neighbourhoodTheta[2][i] += thetaNonAddr[i] * weightResult;
				weightSum[2] += weightResult;
			}

			double[] thetaSender = candidateNode.getSenderTheta(iv);
			if ((thetaSender == null) && (timePeriod != longestTimePeriod))
				thetaSender = candidateNode.getSenderTheta(ivExt);
			if (thetaSender != null) {
				for (int i = 0; i < thetaSender.length; i++)
					neighbourhoodTheta[3][i] += thetaSender[i] * weightResult;
				weightSum[3] += weightResult;
			}
		}

		for (int i = 0; i < neighbourhoodTheta.length; i++) {
			double cdf = 0.0;
			if (weightSum[i] != 0.0) {
				for (int j = 0; j < neighbourhoodTheta[i].length; j++) {
					neighbourhoodTheta[i][j] /= weightSum[i];
					cdf += neighbourhoodTheta[i][j];
				}
			}

			if (cdf == 0.0)
				neighbourhoodTheta[i] = priorDist;
			else if (Math.abs(1.0 - cdf) > EPSILON)
				throw new RuntimeException("CDF sum != 1 (" + cdf + ")");
		}
		WeightFunction.addStats(variant, weights);

		thetaIncomingEdges = neighbourhoodTheta[0];
		thetaOutgoingEdges = neighbourhoodTheta[1];
		thetaOtherUsersNode = neighbourhoodTheta[2];
		thetaOtherUsersEdge = neighbourhoodTheta[3];
	}

	public void setEmptyNeighborhood() {
		thetaIncomingEdges = priorObserve;
		thetaOutgoingEdges = priorObserve;
		thetaOtherUsersNode = priorObserve;
		thetaOtherUsersEdge = priorObserve;
	}

	public long getFromUserId() {
		return fromUserNode.id;
	}

	public long getToUserId() {
		return toUserId;
	}

	public double[] getThetaToPredict() {
		return thetaToPredict;
	}

}
