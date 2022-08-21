package edu.tum.cs.graph;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import edu.uci.ics.jung.graph.UndirectedSparseGraph;

/*
 * Adapted from GraphStream <http://graphstream-project.org> under the terms of the CeCILL-C license.
 * Modified by Jan Hauffa on 25.06.2019.
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
public class BetweennessCentralityTest {

	@Test
	public void test1() {
		//     F---E     Cb(A) = 1
		//    /|    \    Cb(B) = 1
		//   / |     \   Cb(C) = 3
		//  /  |      \  Cb(D) = 3
		// A---C-------D Cb(E) = 1
		//  \  |     _/  Cb(F) = 3
		//   \ |  __/
		//    \|_/
		//     B

		UndirectedSparseGraph<String, String> graph = new UndirectedSparseGraph<String, String>();
		graph.addVertex("A");
		graph.addVertex("B");
		graph.addVertex("C");
		graph.addVertex("D");
		graph.addVertex("E");
		graph.addVertex("F");

		graph.addEdge("AB", "A", "B");
		graph.addEdge("AC", "A", "C");
		graph.addEdge("AF", "A", "F");
		graph.addEdge("BC", "B", "C");
		graph.addEdge("FC", "F", "C");
		graph.addEdge("CD", "C", "D");
		graph.addEdge("FE", "F", "E");
		graph.addEdge("ED", "E", "D");
		graph.addEdge("BD", "B", "D");

		BetweennessCentrality<String, String> bcb = new BetweennessCentrality<String, String>(graph);
		bcb.compute();
		assertEquals(1.0, bcb.getCentrality("A"), 0.0);
		assertEquals(1.0, bcb.getCentrality("B"), 0.0);
		assertEquals(3.0, bcb.getCentrality("C"), 0.0);
		assertEquals(3.0, bcb.getCentrality("D"), 0.0);
		assertEquals(1.0, bcb.getCentrality("E"), 0.0);
		assertEquals(3.0, bcb.getCentrality("F"), 0.0);
	}

}
