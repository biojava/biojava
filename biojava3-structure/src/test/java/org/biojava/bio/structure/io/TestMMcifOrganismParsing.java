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
 * created at Aug 12, 2013
 * Author: ap3 
 */

package org.biojava.bio.structure.io;





import org.biojava.bio.structure.Compound;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.StructureIO;

import junit.framework.TestCase;

public class TestMMcifOrganismParsing extends TestCase{
	
	
	
	@Override
	protected void setUp() throws Exception {
		
		super.setUp();
		
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		StructureIO.setAtomCache(cache);
	}

	public void test1STP(){
		String pdbId = "1stp";
		
		checkPDB(pdbId);
		
	}
	
	public void test1a4w(){
		String pdbId = "1a4w";
		
		checkPDB(pdbId);
		
	}
	
	public void test4hhb(){
		String pdbId = "4hhb";
		
		checkPDB(pdbId);
		
	}
	
	public void test3ZD6(){
		// a PDB ID that contains a synthetic entity
		String pdbId = "3ZD6";
		
		checkPDB(pdbId);
		
	}
	
	
	

	private void checkPDB(String pdbId) {
		try {
		Structure s = StructureIO.getStructure(pdbId);
		
		assertNotNull(s.getCompounds());
		assertTrue(s.getCompounds().size() > 0);
		
		for ( Compound c : s.getCompounds()) {
			assertNotNull(c.getOrganismTaxId());
		}
		
		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}
		
		
		
	}

}
