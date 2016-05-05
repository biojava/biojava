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
package org.biojava.nbio.structure.io;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.junit.Test;

import static org.junit.Assert.*;

public class TestHardBioUnits {

	@Test
	public void test4A1I(){

		String pdbId = "4A1I";
		int biolAssemblyNr = 2;

		Structure bioAssembly;
		try {

			bioAssembly = StructureIO.getBiologicalAssembly(pdbId,biolAssemblyNr);

			if ( bioAssembly == null){
				System.err.println("Could not generate the biological assembly " + pdbId + " nr " + biolAssemblyNr);
			}


			/*
			 * loop_
				_pdbx_struct_assembly_gen.assembly_id
				_pdbx_struct_assembly_gen.oper_expression
				_pdbx_struct_assembly_gen.asym_id_list
				1 1 A,I,J,K,L,M,N,UA,H,PA,QA,RA,SA,TA,BB
				2 1 G,KA,LA,MA,NA,OA,AB
				2 2 B,O,P,Q,R,VA
				3 1 B,O,P,Q,R,VA
				3 3 G,KA,LA,MA,NA,OA,AB
				4 1 C,S,T,U,V,W,WA,F,FA,GA,HA,IA,JA,ZA
				5 1 D,X,Y,Z,XA,E,AA,BA,CA,DA,EA,YA
			 */

			//System.out.println(bioAssembly.toPDB());

			assertTrue(bioAssembly.nrModels() == 2);

			assertTrue(bioAssembly.getChains().size() > 0);

			Chain g = bioAssembly.getChainByPDB("G", 0);

			assertNotNull(g);

			Chain b = bioAssembly.getChainByPDB("B", 1);

			assertNotNull(b);

			assertFalse(bioAssembly.hasChain("A"));

			assertFalse(bioAssembly.hasChain("H"));

			assertEquals(1,bioAssembly.getPolyChains(0).size());
			assertEquals(1,bioAssembly.getPolyChains(1).size());
		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}

	}
}
