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

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

/**
 * Test for a file with some sugar molecules covalently bound to residues
 * having same residue number for sugar and residue (chain T, residue 201)
 * @author Jose Duarte
 *
 */
public class Test3th3 {

	@Test
	public void test3th3() throws StructureException, IOException {
		AtomCache cache = new AtomCache();
		StructureIO.setAtomCache(cache);

		Structure s = StructureIO.getStructure("3th3");

		// there's 2 residues with residue number 201 in chain T: LYS and a sugar BGC
		// that's more like bad annotation in the file but it's good if we can parse
		// the file without crashing and produce a good warning

		// below we make sure that we parse both residues but that we can only lookup the
		// aminoacid residue (see ChainImpl.addChain)

		// since biojava 5.0 polymer and nonpolymer chains are separated, we've modified the
		// test accordingly below
		
		Chain c = s.getPolyChainByPDB("T");

		ResidueNumber rn = ResidueNumber.fromString("201");
		rn.setChainName("T");

		Group g = c.getGroupByPDB(rn);

		// we get the aminoacid residue and not the BGC sugar residue
		assertEquals("LYS", g.getPDBName());


		// let's see if we have both the residues with that number:
		int count = 0;

		for (Group gr : c.getAtomGroups()) {
			if (gr.getResidueNumber().equals(rn)) count++;
		}

		assertEquals(1, count);
		
		c = s.getNonPolyChainsByPDB("T").get(0);
		
		g = c.getGroupByPDB(rn);
		
		assertNotNull(g);
		assertEquals("BGC", g.getPDBName());
	}


}
