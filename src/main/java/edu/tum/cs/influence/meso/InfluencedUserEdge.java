package edu.tum.cs.influence.meso;

import java.util.Date;

import edu.tum.cs.graph.SocialEdge;
import edu.tum.cs.graph.SocialNode;
import edu.tum.cs.time.TimeInterval;

public class InfluencedUserEdge extends InfluencedUser {

	private double[] thetaEdge;
	private double[] thetaReversedEdge;
	private double[] thetaReadingUserNode;
	private double[] thetaReadingUserEdge;

	public InfluencedUserEdge(SocialNode fromUserNode, long toUserId, SocialGraphCache socialGraph, Date centerDate,
			int timePeriod, int longestTimePeriod, boolean reverseTime, double[] priorObserve, double[] priorPredict) {
		super(fromUserNode, toUserId, socialGraph, centerDate, timePeriod, longestTimePeriod, reverseTime, priorObserve,
				priorPredict);
	}

	@Override
	public void setThetas(SocialNode fromUserNode, SocialNode toUserNode, SocialEdge edge, SocialEdge reversedEdge) {
		TimeInterval observe, observeExt;
		TimeInterval predict, predictExt;
		double[] predictPriorDist, observePriorDist;
		if (!reverseTime) {
			observe = new TimeInterval(centerDate, -timePeriod);
			observeExt = new TimeInterval(centerDate, -longestTimePeriod);
			predict = new TimeInterval(centerDate, timePeriod);
			predictExt = new TimeInterval(centerDate, longestTimePeriod);
			observePriorDist = priorObserve;
			predictPriorDist = priorPredict;
		} else {
			observe = new TimeInterval(centerDate, timePeriod);
			observeExt = new TimeInterval(centerDate, longestTimePeriod);
			predict = new TimeInterval(centerDate, -timePeriod);
			predictExt = new TimeInterval(centerDate, -longestTimePeriod);
			predictPriorDist = priorObserve;
			observePriorDist = priorPredict;
		}

		// THETA EDGE
		thetaEdge = edge.getTheta(observe);
		if ((thetaEdge == null) && (timePeriod != longestTimePeriod))
			thetaEdge = edge.getTheta(observeExt);
		if (thetaEdge == null) {
			/* Even though we selected the edge under the assumption that it has at least one associated message within
			   the current time period, that message may be empty after tokenization/preprocessing -> no edge theta. */
			logger.warning("'theta edge' missing for edge " + edge.fromUserId + " -> " + edge.toUserId);
			thetaEdge = observePriorDist;
		}

		// THETA EDGE REVERSED
		thetaReversedEdge = null;
		if (reversedEdge != null) {
			thetaReversedEdge = reversedEdge.getTheta(observe);
			if ((thetaReversedEdge == null) && (timePeriod != longestTimePeriod))
				thetaReversedEdge = reversedEdge.getTheta(observeExt);
		}
		if (thetaReversedEdge == null)
			thetaReversedEdge = observePriorDist;

		// THETA WRITING USER NODE
		thetaWritingUserNode = fromUserNode.getTheta(observe);
		if ((thetaWritingUserNode == null) && (timePeriod != longestTimePeriod))
			thetaWritingUserNode = fromUserNode.getTheta(observeExt);
		if (thetaWritingUserNode == null)
			thetaWritingUserNode = observePriorDist;

		// THETA WRITING USER EDGE
		thetaWritingUserEdge = fromUserNode.getSenderTheta(observe);
		if ((thetaWritingUserEdge == null) && (timePeriod != longestTimePeriod))
			thetaWritingUserEdge = fromUserNode.getSenderTheta(observe);
		if (thetaWritingUserEdge == null)
			thetaWritingUserEdge = observePriorDist;

		// THETA WRITING USER EXPOSURE
		thetaWritingUserExposure = fromUserNode.getExposureTheta(observe);
		if ((thetaWritingUserExposure == null) && (timePeriod != longestTimePeriod))
			thetaWritingUserExposure = fromUserNode.getExposureTheta(observeExt);
		if (thetaWritingUserExposure == null)
			thetaWritingUserExposure = observePriorDist;

		// THETA WRITING USER INCOMING
		thetaWritingUserIncoming = fromUserNode.getRecipientTheta(observe);
		if ((thetaWritingUserIncoming == null) && (timePeriod != longestTimePeriod))
			thetaWritingUserIncoming = fromUserNode.getRecipientTheta(observeExt);
		if (thetaWritingUserIncoming == null)
			thetaWritingUserIncoming = observePriorDist;

		// THETA READING USER INCOMING
		thetaReadingUserNode = toUserNode.getTheta(observe);
		if ((thetaReadingUserNode == null) && (timePeriod != longestTimePeriod))
			thetaReadingUserNode = toUserNode.getTheta(observeExt);
		if (thetaReadingUserNode == null)
			thetaReadingUserNode = observePriorDist;

		// THETA READING USER EDGE
		thetaReadingUserEdge = toUserNode.getSenderTheta(observe);
		if ((thetaReadingUserEdge == null) && (timePeriod != longestTimePeriod))
			thetaReadingUserEdge = toUserNode.getSenderTheta(observeExt);
		if (thetaReadingUserEdge == null)
			thetaReadingUserEdge = observePriorDist;

		// THETA GLOBAL NODE / EDGE
		thetaGlobalNode = socialGraph.getGlobalNode().getTheta(observe);
		thetaGlobalEdge = socialGraph.getGlobalNode().getSenderTheta(observe);

		// THETA TO PREDICT
		thetaToPredict = edge.getTheta(predict);
		if ((thetaToPredict == null) && (timePeriod != longestTimePeriod))
			thetaToPredict = edge.getTheta(predictExt);
		if (thetaToPredict == null) {
			// same issue as above
			logger.warning("'theta to predict' missing for edge " + edge.fromUserId + " -> " + edge.toUserId);
			thetaToPredict = predictPriorDist;
		}
	}

	@Override
	public double[] getBaselineTheta() {
		return thetaEdge;
	}

	@Override
	public double[][] getThetaMatrix() {
		double[][] inst = new double[numTopics][getNumVariables()];
		for (int i = 0; i < inst.length; i++) {
			inst[i][0] = thetaEdge[i];
			inst[i][1] = thetaReversedEdge[i];
			inst[i][2] = thetaWritingUserNode[i];
			inst[i][3] = thetaWritingUserEdge[i];
			inst[i][4] = thetaWritingUserExposure[i];
			inst[i][5] = thetaWritingUserIncoming[i];
			inst[i][6] = thetaReadingUserNode[i];
			inst[i][7] = thetaReadingUserEdge[i];
			inst[i][8] = thetaIncomingEdges[i];
			inst[i][9] = thetaOutgoingEdges[i];
			inst[i][10] = thetaOtherUsersNode[i];
			inst[i][11] = thetaOtherUsersEdge[i];
			inst[i][12] = thetaGlobalNode[i];
			inst[i][13] = thetaGlobalEdge[i];
		}
		return inst;
	}

	@Override
	public int getNumVariables() {
		return 14;
	}

}
