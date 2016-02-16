package org.biojava.nbio.structure;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.Test;

/**
 * This tests the correctness of the {@link Structure} data structure, in terms
 * of correct implementation, given the expected behavior, of the methods.
 * 
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class TestStructureImpl {

	/**
	 * This tests the correctness of the Links between the different objects
	 * that create a Structure (Chains, Groups, Atoms). The correct behavior is
	 * that objects higher in the hierarchy should set the links of lower level
	 * objects.
	 */
	@Test
	public void testLinks() {

		// Create new Structure and assign one Chain to it
		Structure s1 = new StructureImpl();
		s1.addModel(new ArrayList<Chain>(1));

		Chain c1 = new ChainImpl();
		s1.addChain(c1);

		// Test that chains parent is s1
		assertEquals(s1, c1.getStructure());

		// Populate the Chain with one Group
		Group g1 = new HetatomImpl();
		c1.addGroup(g1);

		// Test that the group parent is c1
		assertEquals(c1, g1.getChain());

		// Add a single Atom to the Group
		Atom a1 = new AtomImpl();
		g1.addAtom(a1);

		// Test that the atom parent is g1
		assertEquals(c1, g1.getChain());

		// Now clone the Atom and assign it to a second structure
		Atom a2 = (Atom) a1.clone();

		Structure s2 = new StructureImpl();
		s2.addModel(new ArrayList<Chain>(1));
		Chain c2 = new ChainImpl();
		s2.addChain(c2);
		Group g2 = new HetatomImpl();
		c2.addGroup(g2);
		g2.addAtom(a2);

		// Test correct parent links in Atoms
		assertEquals(g1, a1.getGroup());
		assertEquals(g2, a2.getGroup());

		// Now clone the Group and assign it to new third structure
		Group g3 = (Group) g1.clone();
		Atom a3 = g3.getAtom(0);
		Structure s3 = new StructureImpl();
		s3.addModel(new ArrayList<Chain>(1));
		Chain c3 = new ChainImpl();
		s3.addChain(c3);
		c3.addGroup(g3);

		// Test correct parent links in Groups and Atoms
		assertEquals(c1, g1.getChain());
		assertEquals(g1, a1.getGroup());
		assertEquals(c3, g3.getChain());
		assertEquals(g3, a3.getGroup());

	}
}
