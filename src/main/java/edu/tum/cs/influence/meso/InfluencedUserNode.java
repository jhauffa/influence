package edu.tum.cs.influence.meso;

import java.util.Date;

import edu.tum.cs.graph.SocialEdge;
import edu.tum.cs.graph.SocialNode;
import edu.tum.cs.time.TimeInterval;

public class InfluencedUserNode extends InfluencedUser {

	public InfluencedUserNode(SocialNode fromUserNode, long toUserId, SocialGraphCache socialGraph, Date centerDate,
			int timePeriod, int longestTimePeriod, boolean reverseTime, double[] priorObserve, double[] priorPredict) {
		super(fromUserNode, toUserId, socialGraph, centerDate, timePeriod, longestTimePeriod, reverseTime, priorObserve,
				priorPredict);
	}

	@Override
	public void setThetas(SocialNode fromUserNode, SocialNode toUserNode, SocialEdge edge, SocialEdge reversedEdge) {
		setThetas(fromUserNode);	// the other parameters are not appropriate for nodes
	}

	public void setThetas(SocialNode fromUserNode) {
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

		// THETA WRITING USER NODE
		thetaWritingUserNode = fromUserNode.getTheta(observe);
		if ((thetaWritingUserNode == null) && (timePeriod != longestTimePeriod))
			thetaWritingUserNode = fromUserNode.getTheta(observeExt);
		if (thetaWritingUserNode == null) {
			/* Even though we selected the node under the assumption that it has sent at least one non-addressive
			   message within the current time period, that message may be empty after tokenization/preprocessing -> no
			   node theta. */
			logger.warning("'theta writing user node' missing for node " + fromUserNode.id);
			thetaWritingUserNode = observePriorDist;
		}

		// THETA WRITING USER EDGE
		thetaWritingUserEdge = fromUserNode.getSenderTheta(observe);
		if ((thetaWritingUserEdge == null) && (timePeriod != longestTimePeriod))
			thetaWritingUserEdge = fromUserNode.getSenderTheta(observeExt);
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

		// THETA GLOBAL NODE / EDGE
		thetaGlobalNode = socialGraph.getGlobalNode().getTheta(observe);
		thetaGlobalEdge = socialGraph.getGlobalNode().getSenderTheta(observe);

		// THETA TO PREDICT
		thetaToPredict = fromUserNode.getTheta(predict);
		if ((thetaToPredict == null) && (timePeriod != longestTimePeriod))
			thetaToPredict = fromUserNode.getTheta(predictExt);
		if (thetaToPredict == null) {
			// same issue as above
			logger.warning("'theta to predict' missing for node " + fromUserNode.id);
			thetaToPredict = predictPriorDist;
		}
	}

	@Override
	public double[] getBaselineTheta() {
		return thetaWritingUserNode;
	}

	@Override
	public double[][] getThetaMatrix() {
		double[][] inst = new double[numTopics][getNumVariables()];
		for (int i = 0; i < inst.length; i++) {
			inst[i][0] = thetaWritingUserNode[i];
			inst[i][1] = thetaWritingUserEdge[i];
			inst[i][2] = thetaWritingUserExposure[i];
			inst[i][3] = thetaWritingUserIncoming[i];
			inst[i][4] = thetaIncomingEdges[i];
			inst[i][5] = thetaOutgoingEdges[i];
			inst[i][6] = thetaOtherUsersNode[i];
			inst[i][7] = thetaOtherUsersEdge[i];
			inst[i][8] = thetaGlobalNode[i];
			inst[i][9] = thetaGlobalEdge[i];
		}
		return inst;
	}

	@Override
	public int getNumVariables() {
		return 10;
	}

}
