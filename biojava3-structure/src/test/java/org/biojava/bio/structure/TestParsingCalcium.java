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
 * Created on Nov 1, 2013
 * Author: andreas 
 *
 */

package org.biojava.bio.structure;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava3.structure.StructureIO;

import junit.framework.TestCase;

public class TestParsingCalcium extends TestCase{


	public void testCalciumParsing(){

		String pdbID = "1SU4";

		// Calcium is at position 995
		// HETATM 7673 CA    CA A 995      64.194  12.588   7.315  1.00 41.55          CA  
		
		
		try {
			AtomCache cache = new AtomCache();
			Structure s = cache.getStructure(pdbID);
			cache.setUseMmCif(true);

			Structure m = cache.getStructure(pdbID);
			
			Group g1 = s.getChainByPDB("A").getGroupByPDB(new ResidueNumber("A",995,null));
			Group g2 = m.getChainByPDB("A").getGroupByPDB(new ResidueNumber("A",995,null));
			
			// can't do that! the atom index is not the same!
			//assertEquals(g1.getAtom(0).toPDB(),g2.getAtom(0).toPDB());
			
			assertEquals(g1.getAtom(0).getFullName(),"CA  ");
			assertEquals(g1.getAtom(0).getFullName(), g2.getAtom(0).getFullName());
		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}


	}
}
