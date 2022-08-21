package edu.tum.cs.graph;

import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

import edu.uci.ics.jung.graph.Graph;

/*
 * Adapted from GraphStream <http://graphstream-project.org> under the terms of the CeCILL-C license.
 * Modified by Jan Hauffa on 25.06.2019.
 *
 * Computation of betweenness centrality is slow, so we imported the implementation of GraphStream, which is ~1.6 times
 * faster than JUNG.
 */
/*
 * Copyright 2006 - 2016
 *     Stefan Balev     <stefan.balev@graphstream-project.org>
 *     Julien Baudry    <julien.baudry@graphstream-project.org>
 *     Antoine Dutot    <antoine.dutot@graphstream-project.org>
 *     Yoann Pigné      <yoann.pigne@graphstream-project.org>
 *     Guilhelm Savin   <guilhelm.savin@graphstream-project.org>
 * 
 * This program is free software distributed under the terms of two licenses, the
 * CeCILL-C license that fits European law, and the GNU Lesser General Public
 * License. You can  use, modify and/ or redistribute the software under the terms
 * of the CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
 * URL <http://www.cecill.info> or under the terms of the GNU LGPL as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
public class BetweennessCentrality<V, E> {

	private class NodeData {
		public int distance = Integer.MAX_VALUE;
		public double sigma = 0.0, delta = 0.0, centrality = 0.0;
		public Set<V> predecessors = new HashSet<V>();
	}

	private final Graph<V, E> graph;
	private final Map<V, NodeData> data = new HashMap<V, NodeData>();

	private class BrandesNodeComparatorLargerFirst implements Comparator<V> {
		public int compare(V x, V y) {
			int yy = data.get(y).distance;
			int xx = data.get(x).distance;

			if (xx > yy)
				return -1;
			else if (xx < yy)
				return 1;
			return 0;
		}
	}

	public BetweennessCentrality(Graph<V, E> graph) {
		this.graph = graph;
	}

	public double getCentrality(V v) {
		return data.get(v).centrality;
	}

	public void compute() {
		for (V node : graph.getVertices())
			data.put(node, new NodeData());

		for (V s : graph.getVertices()) {
			PriorityQueue<V> q = explore(s);

			// accumulation phase of Brandes' algorithm
			while (!q.isEmpty()) {
				V w = q.poll();
				NodeData wData = data.get(w);

				for (V v : wData.predecessors) {
					NodeData vData = data.get(v);
					double c = ((vData.sigma / wData.sigma) * (1.0 + wData.delta));
					vData.delta += c;
				}
				if (w != s)
					wData.centrality += wData.delta;
			}

			// reset per-node data
			for (V v : graph.getVertices()) {
				NodeData vData = data.get(v);
				vData.predecessors.clear();
				vData.sigma = vData.delta = 0.0;
				vData.distance = Integer.MAX_VALUE;
			}
		}
	}

	private PriorityQueue<V> explore(V source) {
		LinkedList<V> q = new LinkedList<V>();
		PriorityQueue<V> s = new PriorityQueue<V>(graph.getVertexCount(), new BrandesNodeComparatorLargerFirst());

		q.add(source);
		NodeData sourceData = data.get(source);
		sourceData.sigma = 1.0;
		sourceData.distance = 0;

		while (!q.isEmpty()) {
			V v = q.removeFirst();
			NodeData vData = data.get(v);
			s.add(v);

			for (E l : graph.getOutEdges(v)) {
				V w = graph.getOpposite(v, l);
				NodeData wData = data.get(w);

				if (wData.distance == Integer.MAX_VALUE) {
					wData.distance = vData.distance + 1;
					q.add(w);
				}

				if (wData.distance == (vData.distance + 1)) {
					wData.sigma += vData.sigma;
					wData.predecessors.add(v);
				}
			}
		}
		return s;
	}

}
