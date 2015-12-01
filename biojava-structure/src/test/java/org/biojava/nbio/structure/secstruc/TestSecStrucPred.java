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
package org.biojava.nbio.structure.secstruc;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Test the correctness of the DSSP implementation in BioJava
 * for the prediction of secondary structure in a Structure object.
 * 
 * EXAMPLES:
 * 			Big structures: 4v7r, 4V60 (use mmCif parser)
 * 			Helical: 4hhb, 4lup
 * 			Mixed small: 5pti
 * 			First sheet: 1ze3, 3k19
 * 			Insertion code: 1how
 * 
 * @author Aleix Lafita
 *
 */
public class TestSecStrucPred {

	@Test
	public void testSecStrucPred() throws StructureException, IOException {
		
		//List of names to test the DSSP prediction
		List<String> names = Arrays.asList(
				"5pti", "1tim", "4hhb", "1how", "4i4q");
		
		for (String name : names) {
			
			AtomCache cache = new AtomCache();		
			Structure s = cache.getStructure(name);
			
			//Predict with BioJava the SS
			SecStrucPred sec = new SecStrucPred();
			List<SecStrucState> biojava = sec.predict(s, true);
			
			//Download the original DSSP implementation output
			List<SecStrucState> dssp = DSSPParser.fetch(name, s, false);
			
			assertEquals("SS assignment lengths do not match",
					biojava.size(), dssp.size());
			
			for (int i=0; i<dssp.size(); i++){
				assertEquals("SS assignment position "+(i+1)+" does not match", 
						biojava.get(i), dssp.get(i));
			}
		}
	}
		
}
