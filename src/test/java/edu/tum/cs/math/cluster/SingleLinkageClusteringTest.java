package edu.tum.cs.math.cluster;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Test;

public class SingleLinkageClusteringTest {

	private static class DummyClustering extends SingleLinkageClustering<Integer> {
		private static final long serialVersionUID = 1L;

		public DummyClustering(List<Integer> data) {
			super(data, false);
		}

		@Override
		protected double distance(Integer e1, Integer e2) {
			return Math.abs(e1 - e2);
		}
	}

	public static void verifyClustering(SingleLinkageClustering<Integer> clust, Collection<Integer> data) {
		Set<Integer> dataSet = new HashSet<Integer>(data);

		Collection<Collection<Integer>> prevClusters = null;
		for (int i = 0; i <= clust.getMaxLevel(); i++) {
			Collection<Collection<Integer>> clusters = clust.getClusters(i);
			if (i == 0) {
				assertEquals(data.size(), clusters.size());
			} else {
				if (i == clust.getMaxLevel())
					assertEquals(1, clusters.size());

				for (Collection<Integer> c1 : prevClusters) {
					boolean hasMatch = false;
					for (Collection<Integer> c2 : clusters) {
						if (c2.containsAll(c1)) {
							hasMatch = true;
							break;
						}
					}
					assertTrue(hasMatch);
				}
			}

			Set<Integer> allClusterData = new HashSet<Integer>();
			for (Collection<Integer> cluster : clusters)
				allClusterData.addAll(cluster);
			assertEquals(dataSet, allClusterData);

			prevClusters = clusters;
		}
	}

	@Test
	public void testSingleLinkageClustering() {
		Integer[] data = new Integer[] { 1, 10, 7, 4, 5, 3, 6, 8, 2, 9 };
		DummyClustering clust = new DummyClustering(Arrays.asList(data));
		clust.process();
		verifyClustering(clust, Arrays.asList(data));
		assertEquals(1, clust.getMaxLevel());

		data = new Integer[] { 1, 46, 22, 7, 11, 4, 16, 29, 2, 37 };
		clust = new DummyClustering(Arrays.asList(data));
		clust.process();
		verifyClustering(clust, Arrays.asList(data));
		assertEquals(data.length - 1, clust.getMaxLevel());
	}

}
