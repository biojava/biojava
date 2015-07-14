package org.biojava.nbio.structure.symmetry;

import static org.junit.Assert.*;

import java.util.Iterator;

import org.biojava.nbio.structure.symmetry.utils.DirectedAcyclicGraph;
import org.biojava.nbio.structure.symmetry.utils.DirectedGraph;
import org.biojava.nbio.structure.symmetry.utils.Graph;
import org.biojava.nbio.structure.symmetry.utils.GraphIteratorBFS;
import org.biojava.nbio.structure.symmetry.utils.GraphIteratorDFS;
import org.biojava.nbio.structure.symmetry.utils.SimpleGraph;
import org.junit.Test;

/**
 * Test the correctness of the {@link Graph} implementations
 * and Graph iterators.<p>
 * Currently testing for undirected, directed and DAG graph 
 * implementations and for DFS and BFS graph iterators.
 * 
 * @author Aleix Lafita
 *
 */
public class TestGraphIterators {

	@Test
	public void testDFS(){

		//Test directed graph
		Graph<Character> g1 = createDirectedGraph();
		Iterator<Character> iter = new GraphIteratorDFS<Character>(g1, 0);
		String result = "";
		while (iter.hasNext()){
			result += iter.next();
		}
		assertEquals("ACDEB", result);

		//Test DAG
		Graph<Character> g2 = createDAG();
		Iterator<Character> iter2 = new GraphIteratorDFS<Character>(g2, 0);
		String result2 = "";
		while (iter2.hasNext()){
			result2 += iter2.next();
		}
		assertEquals("ABCDE", result2);

		//Test undirected graph
		Graph<Character> g3 = createUndirectedGraph();
		Iterator<Character> iter3 = new GraphIteratorDFS<Character>(g3, 0);
		String result3 = "";
		while (iter3.hasNext()){
			result3 += iter3.next();
		}
		assertEquals("ACBDE", result3);
	}

	@Test
	public void testBFS(){
		
		//Test directed graph
		Graph<Character> g1 = createDirectedGraph();
		Iterator<Character> iter = new GraphIteratorBFS<Character>(g1, 0);
		String result = "";
		while (iter.hasNext()){
			result += iter.next();
		}
		assertEquals("ABCDE", result);

		//Test DAG
		Graph<Character> g2 = createDAG();
		Iterator<Character> iter2 = new GraphIteratorBFS<Character>(g2, 0);
		String result2 = "";
		while (iter2.hasNext()){
			result2 += iter2.next();
		}
		assertEquals("ABCDE", result2);

		//Test undirected graph
		Graph<Character> g3 = createUndirectedGraph();
		Iterator<Character> iter3 = new GraphIteratorBFS<Character>(g3, 0);
		String result3 = "";
		while (iter3.hasNext()){
			result3 += iter3.next();
		}
		assertEquals("ABCDE", result3);
	}

	/**
	 * Create a DAG consisting of 5 character nodes: A,B,C,D,E.
	 * They are connected in their alphabetical order.
	 */
	private Graph<Character> createDAG(){

		Graph<Character> graph = new DirectedAcyclicGraph<Character>();

		//Add vertices
		char c = 'A';
		for (int i=0; i<5; i++){
			graph.addVertex(c);
			c++;
		}
		//Add edges
		graph.addEdge('A', 'B');
		graph.addEdge('B', 'C');
		graph.addEdge('C', 'D');
		graph.addEdge('D', 'E');

		return graph;
	}

	/**
	 * Create a directed graph consisting of 5 character 
	 * nodes: A,B,C,D,E.
	 * They are connected randomly.
	 */
	private Graph<Character> createDirectedGraph(){

		Graph<Character> graph = new DirectedGraph<Character>();

		//Add vertices
		char c = 'A';
		for (int i=0; i<5; i++){
			graph.addVertex(c);
			c++;
		}
		//Add edges
		graph.addEdge('A', 'B');
		graph.addEdge('A', 'C');
		graph.addEdge('C', 'D');
		graph.addEdge('D', 'B');
		graph.addEdge('D', 'E');

		return graph;
	}

	/**
	 * Create an undirected graph consisting of 5 character 
	 * nodes: A,B,C,D,E.
	 * They are connected randomly.
	 */
	private Graph<Character> createUndirectedGraph(){

		Graph<Character> graph = new SimpleGraph<Character>();

		//Add vertices
		char c = 'A';
		for (int i=0; i<5; i++){
			graph.addVertex(c);
			c++;
		}
		//Add edges
		graph.addEdge('A', 'B');
		graph.addEdge('B', 'C');
		graph.addEdge('C', 'A');
		graph.addEdge('B', 'D');
		graph.addEdge('D', 'E');

		return graph;
	}

}
