package org.biojava.bio.structure.io;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.quaternary.io.BioUnitDataProviderFactory;
import org.biojava3.structure.StructureIO;
import org.junit.Test;
import static org.junit.Assert.*;

public class TestHardBioUnits {
	
	@Test
	public void test4A1I(){
		
		String pdbId = "4A1I";
		int biolAssemblyNr = 2;
		
		BioUnitDataProviderFactory.setBioUnitDataProvider(BioUnitDataProviderFactory.pdbProviderClassName);
		
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
			
			
			
			assertTrue(bioAssembly.getChains().size() > 1);
			
			Chain g = bioAssembly.getChainByPDB("G");
			
			assertNotNull(g);
			
			Chain b = bioAssembly.getChainByPDB("B");
			
			assertNotNull(b);
			
			assertFalse(bioAssembly.hasChain("A"));
			
			assertFalse(bioAssembly.hasChain("H"));
			
			
			
			assertEquals(2,bioAssembly.getChains().size());
		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}

	}
}
