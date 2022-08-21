package edu.tum.cs.graph.clique;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.collections15.Factory;
import org.junit.Test;

import com.carrotsearch.hppc.ShortHashSet;
import com.carrotsearch.hppc.ShortSet;

import edu.uci.ics.jung.algorithms.generators.random.ErdosRenyiGenerator;
import edu.uci.ics.jung.graph.UndirectedGraph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;

public class CliquePercolationTest {

	private static short computeOverlap(short[] c1, short[] c2) {
		if (c1.length > c2.length) {
			short[] tmp = c1;
			c1 = c2;
			c2 = tmp;
		}

		short n = 0;
		for (int i = 0; i < c1.length; i++)
			if (Arrays.binarySearch(c2, c1[i]) >= 0)
				n++;
		return n;
	}

	private static short[][] buildCliqueOverlapMatrix(List<short[]> cliques) {
		int n = cliques.size();
		short[][] cliqueOverlap = new short[n - 1][];
		for (int i = 0; i < (n - 1); i++) {
			cliqueOverlap[i] = new short[i + 1];
			short[] clique = cliques.get(i + 1);
			for (int j = 0; j < (i + 1); j++)
				cliqueOverlap[i][j] = computeOverlap(clique, cliques.get(j));
		}
		return cliqueOverlap;
	}

	private static Collection<short[]> getReferenceCommunities(List<short[]> cliques, int size) {
		Collection<short[]> communities = new ArrayList<short[]>();

		for (short[] clique : cliques)
			Arrays.sort(clique);
		short[][] cliqueOverlap = buildCliqueOverlapMatrix(cliques);

		int n = cliques.size();
		int[] label = new int[n];
		for (int i = 0; i < n; i++) {
			ShortSet community = new ShortHashSet();
			LinkedList<Integer> queue = new LinkedList<Integer>();
			queue.add(i);
			while (!queue.isEmpty()) {
				int j = queue.poll();
				if ((label[j] != 0) || (cliques.get(j).length < size))
					continue;

				label[j] = i + 1;
				for (short m : cliques.get(j))
					community.add(m);

				for (int k = 0; k < j; k++)
					if (cliqueOverlap[j - 1][k] >= (size - 1))
						queue.add(k);
				for (int k = j + 1; k < n; k++)
					if (cliqueOverlap[k - 1][j] >= (size - 1))
						queue.add(k);
			}
			if (!community.isEmpty())
				communities.add(community.toArray());
		}

		return communities;
	}

	private static void verifyCommunities(Collection<short[]> communities, Set<Set<Short>> refCommunities) {
		assertEquals(refCommunities.size(), communities.size());
		for (short[] community : communities) {
			Set<Short> communitySet = new HashSet<Short>();
			for (short v : community)
				communitySet.add(v);
			assertTrue(refCommunities.contains(communitySet));
		}
	}

	@Test
	public void testCliquePercolation() {
		short[][] cliques = {
			{ 0, 1, 2, 3 }, { 1, 2, 3, 4 }, { 5, 8, 9 }, { 5, 6, 8 }, { 6, 7, 8 }, { 0, 5 }
		};
		Short[][] communities3 = {
			{ 0, 1, 2, 3, 4 }, { 5, 6, 7, 8, 9 }
		};

		Set<Set<Short>> communities3Set = new HashSet<Set<Short>>();
		for (Short[] community : communities3)
			communities3Set.add(new HashSet<Short>(Arrays.asList(community)));

		verifyCommunities(getReferenceCommunities(Arrays.asList(cliques), 3), communities3Set);

		CliquePercolation cmf = new CliquePercolation(3, false);
		cmf.addAllCliques(Arrays.asList(cliques));
		Collection<short[]> communities = cmf.getCommunities();
		verifyCommunities(communities, communities3Set);
		cmf = new CliquePercolation(3, true);
		cmf.addAllCliques(Arrays.asList(cliques));
		communities = cmf.getCommunities();
		verifyCommunities(communities, communities3Set);

		CliquePercolationParallel cmp = new CliquePercolationParallel(3, false);
		cmp.addAllCliques(Arrays.asList(cliques));
		communities = cmp.getCommunities();
		verifyCommunities(communities, communities3Set);
		cmp = new CliquePercolationParallel(3, true);
		cmp.addAllCliques(Arrays.asList(cliques));
		communities = cmp.getCommunities();
		verifyCommunities(communities, communities3Set);
	}

	private static final int numRandomGraphs = 1000;
	private static final int numVertices = 25;
	private static final double edgeProb = 0.2;

	private static int edgeIdx;
	private static int vertexIdx;

	@Test
	public void testRandomGraph() {
		ErdosRenyiGenerator<Integer, Integer> gen = new ErdosRenyiGenerator<Integer, Integer>(
				new Factory<UndirectedGraph<Integer, Integer>>() {
					@Override
					public UndirectedGraph<Integer, Integer> create() {
						return new UndirectedSparseGraph<Integer, Integer>();
					}
				},
				new Factory<Integer>() {
					@Override
					public Integer create() {
						return new Integer(vertexIdx++);
					}
				},
				new Factory<Integer>() {
					@Override
					public Integer create() {
						return new Integer(edgeIdx++);
					}
				},
				numVertices, edgeProb);
		gen.setSeed(123123L);

		for (int i = 0; i < numRandomGraphs; i++) {
			edgeIdx = vertexIdx = 0;
			UndirectedGraph<Integer, Integer> graph = (UndirectedGraph<Integer, Integer>) gen.create();

			VertexIndex<Integer> index = new VertexIndex<Integer>(graph);
			MaximalClique<Integer, Integer> mc = new MaximalClique<Integer, Integer>(graph, index);
			ArrayList<short[]> cliques = mc.getAllMaximalSubgraphs();

			Collection<short[]> refCommunities = getReferenceCommunities(cliques, 3);
			Set<Set<Short>> refCommunitiesSet = new HashSet<Set<Short>>();
			for (short[] community : refCommunities) {
				Set<Short> communitySet = new HashSet<Short>();
				for (short v : community)
					communitySet.add(v);
				refCommunitiesSet.add(communitySet);
			}

			CliquePercolation cmf = new CliquePercolation(3, false);
			cmf.addAllCliques(cliques);
			verifyCommunities(cmf.getCommunities(), refCommunitiesSet);

			CliquePercolationParallel cmp = new CliquePercolationParallel(3, false);
			cmp.addAllCliques(cliques);
			verifyCommunities(cmp.getCommunities(), refCommunitiesSet);
		}
	}

}
