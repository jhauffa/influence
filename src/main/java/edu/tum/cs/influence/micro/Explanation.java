package edu.tum.cs.influence.micro;

import java.util.Arrays;
import java.util.Comparator;

public class Explanation {

	public boolean haveData = true;
	public int timePeriodIdx;
	public double[] topicsInfluencer, topicsInfluencee;
	public double localMagnitude, globalMagnitude;

	private static Integer[] sortTopics(final double[] topics) {
		Integer[] map = new Integer[topics.length];
		for (int i = 0; i < map.length; i++)
			map[i] = i;
		Arrays.sort(map, new Comparator<Integer>() {
			@Override
			public int compare(Integer o1, Integer o2) {
				return Double.compare(topics[o2], topics[o1]);
			}
		});
		return map;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		if (haveData) {
			if (timePeriodIdx > 0)
				sb.append("at time step ").append(timePeriodIdx);
			sb.append("\nprev. topic distribution of influencer:\n");
			Integer[] map = sortTopics(topicsInfluencer);
			for (int idx : map) {
				if (topicsInfluencer[idx] == 0.0)
					break;
				sb.append(idx).append(';').append(topicsInfluencer[idx]).append("\n");
			}
			sb.append("\ncur. topic distribution of influencee:\n");
			map = sortTopics(topicsInfluencee);
			for (int idx : map) {
				if (topicsInfluencee[idx] == 0.0)
					break;
				sb.append(idx).append(';').append(topicsInfluencee[idx]).append("\n");
			}
			sb.append("\nlocal magnitude = ").append(localMagnitude);
			sb.append(", global magnitude = ").append(globalMagnitude).append("\n");
		} else {
			sb.append("no data");
		}
		return sb.toString();
	}

}
