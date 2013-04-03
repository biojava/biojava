package org.biojava.bio.structure;

import junit.framework.TestCase;


public class TestVeryLongFileName extends TestCase{


	public void testVeryLongFilename(){

		try {
			PDBHeader header = new PDBHeader();

			header.setTitle("IAMAVERYLONGTITLEWITHOUTANYSPACECHARACTERSJUSTTOMAKESUREWECANTESTWHATISGOINGTOHAPPENIFWETRYTOWRITETHISTOAPDBFILE.");

			header.toPDB();

			String title2 = "jCE V.1.1 : file:/Users/ap3/tmp/mareike/IAMAVERYLONGTITLEWITHOUTANYSPACECHARACTERSJUSTTOMAKESUREWECANTESTWHATISGOINGTOHAPPENIFWETRYTOWRITETHISTOAPDBFILE/pdb1_chainG.pdb vs. file:/Users/ap3/tmp/mareike/IAMAVERYLONGTITLEWITHOUTANYSPACECHARACTERSJUSTTOMAKESUREWECANTESTWHATISGOINGTOHAPPENIFWETRYTOWRITETHISTOAPDBFILE/pdb2_chainB.pdb ";
			header.setTitle(title2);
			header.toPDB();
		} catch(Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
}
