/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure.test.io;

import static org.junit.Assert.*;

import java.io.IOException;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

public class TestNcsOpsParsing {

	/**
	 * A structure (viral capsid) with struct_ncs_ops
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void test1auy() throws IOException, StructureException {
		
		AtomCache cache = new AtomCache();
		StructureIO.setAtomCache(cache);

		Structure s = StructureIO.getStructure("1auy");
		
		Matrix4d[] ops = s.getPDBHeader().getCrystallographicInfo().getNcsOperators();
		
		assertNotNull(ops);
		
		// the given operator must not be in our list, only the "generate" operators
		assertEquals(14, ops.length);
		
		for (Matrix4d op:ops) {
			assertEquals(0, op.m30, 0.000001);
			assertEquals(0, op.m31, 0.000001);
			assertEquals(0, op.m32, 0.000001);
			assertEquals(1, op.m33, 0.000001);
		}

	}
	
	/**
	 * A structure with struct_ncs_ops which is not a viral capsid
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void test1a37() throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		StructureIO.setAtomCache(cache);

		Structure s = StructureIO.getStructure("1a37");
		
		Matrix4d[] ops = s.getPDBHeader().getCrystallographicInfo().getNcsOperators();
		
		assertNotNull(ops);
		
		// the given operator must not be in our list, only the "generate" operators
		assertEquals(1, ops.length);
		
	}
	
	/**
	 * A structure without struct_ncs_ops
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void test1smt() throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		StructureIO.setAtomCache(cache);

		Structure s = StructureIO.getStructure("1smt");
		
		Matrix4d[] ops = s.getPDBHeader().getCrystallographicInfo().getNcsOperators();
		
		assertNull(ops);
		
		
		
	}
	

}
