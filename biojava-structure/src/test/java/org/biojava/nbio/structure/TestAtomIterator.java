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
package org.biojava.nbio.structure;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.NoSuchElementException;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomIterator;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.junit.Test;

public class TestAtomIterator {
	@Test
	public void test5frf() throws IOException, StructureException {
		// 5frf: 10 models; residues -2-105, binds a ZN; 1615 atoms/model
		Structure s = StructureIO.getStructure("5frf");
		assertEquals("nrModels",10,s.nrModels());
		
		Atom[] allAtomArray = StructureTools.getAllAtomArray(s);
		assertEquals("getAllAtomArray length",16150, allAtomArray.length);
		
		int atoms=0;
		AtomIterator atomIt = new AtomIterator(s);
		while(atomIt.hasNext()) {
			atoms++;
			atomIt.next();
		}
		try {
			atomIt.next();
			fail("No more elements");
		} catch( NoSuchElementException e) {}
		assertEquals("AtomIterator full length",16150, atoms);
		
		atoms=0;
		atomIt = new AtomIterator(s,0);
		while(atomIt.hasNext()) {
			atoms++;
			atomIt.next();
		}
		try {
			atomIt.next();
			fail("No more elements");
		} catch( NoSuchElementException e) {}
		assertEquals("AtomIterator single model",1615, atoms);
	}
}
