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
import org.biojava.nbio.structure.io.FileParsingParameters;
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
		
		FileParsingParameters params = cache.getFileParsingParams();
		params.setUseInternalChainId(false);
		params.setCreateAtomBonds(true);
		StructureIO.setAtomCache(cache);
		
		Structure s = StructureIO.getStructure("3th3");
		
		// there's 2 residues with residue number 201 in chain T: LYS and a sugar BGC
		// that's more like bad annotation in the file but it's good if we can parse 
		// the file without crashing and produce a good warning
		
		// below we make sure that we parse both residues but that we can only lookup the 
		// aminoacid residue (see ChainImpl.addChain)
		
		Chain c = s.getChainByPDB("T");
		
		ResidueNumber rn = ResidueNumber.fromString("201");
		rn.setChainId("T");
		
		Group g = c.getGroupByPDB(rn);
				
		// we get the aminoacid residue and not the BGC sugar residue
		assertEquals("LYS", g.getPDBName());
		
		
		// let's see if we have both the residues with that number:
		int count = 0;
		
		for (Group gr : c.getAtomGroups()) {
			if (gr.getResidueNumber().equals(rn)) count++;
		}
		
		assertEquals(2, count);
	}

	
}
