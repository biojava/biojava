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
 * Test for a difficult for a file some sugar molecules covalently bound to residues
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
		
		Chain c = s.getChainByPDB("T");
		
		ResidueNumber rn = ResidueNumber.fromString("201");
		rn.setChainId("T");
		
		Group g = c.getGroupByPDB(rn);
				
		assertEquals("LYS", g.getPDBName());
		
		int count = 0;
		
		for (Group gr : c.getAtomGroups()) {
			if (gr.getResidueNumber().equals(rn)) count++;
		}
		
		assertEquals(2, count);
	}

	
}
