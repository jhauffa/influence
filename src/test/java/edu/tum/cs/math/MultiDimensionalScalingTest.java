package edu.tum.cs.math;

import static org.junit.Assert.assertEquals;

import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Test;

public class MultiDimensionalScalingTest {

	// city distances from Wickelmaier, 2003
	private static final double[][] data = {
			{ 0, 93, 82, 133 },
			{ 93, 0, 52, 60 },
			{ 82, 52, 0, 111 },
			{ 133, 60, 111, 0 }
	};

	private static double euclideanDist(double x1, double y1, double x2, double y2) {
		return Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	}

	@Test
	public void testScaling() {
		MultiDimensionalScaling mds = new MultiDimensionalScaling(data);
		RealMatrix emb = mds.getEmbedding(2);
		assertEquals(data.length, emb.getRowDimension());
		assertEquals(2, emb.getColumnDimension());

		for (int i = 0; i < data.length; i++) {
			for (int j = i + 1; j < data.length; j++) {
				double d1 = euclideanDist(emb.getEntry(i, 0), emb.getEntry(i, 1),
						emb.getEntry(j, 0), emb.getEntry(j, 1));
				for (int m = 0; m < data.length; m++) {
					for (int n = m + 1; n < data.length; n++) {
						double d2 = euclideanDist(emb.getEntry(m, 0), emb.getEntry(m, 1),
								emb.getEntry(n, 0), emb.getEntry(n, 1));
						assertEquals((data[i][j] < data[m][n]), (d1 < d2));
					}
				}
			}
		}
	}

}
