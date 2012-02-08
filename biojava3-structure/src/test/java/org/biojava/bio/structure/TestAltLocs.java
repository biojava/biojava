package org.biojava.bio.structure;

import java.util.List;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileConvert;
import org.biojava.bio.structure.io.mmcif.chem.PolymerType;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;

import junit.framework.TestCase;

public class TestAltLocs extends TestCase {

	public void testAltLocParsing(){

		try {
			AtomCache cache = new AtomCache();
			Structure s = cache.getStructure("2CI1");

			//System.out.println(s);

			Chain a = s.getChainByPDB("A");
			//System.out.println(a);

			int groupCount = 0;
			List<Group> groups = a.getAtomGroups();
			for (Group g : groups){
				ChemComp cc = g.getChemComp();
				if ( ResidueType.lPeptideLinking.equals(cc.getResidueType()) ||
						PolymerType.PROTEIN_ONLY.contains(cc.getPolymerType()) ||
						PolymerType.POLYNUCLEOTIDE_ONLY.contains(cc.getPolymerType())
						){
					if (! ChainImpl.waternames.contains(g.getPDBName())) {
						//System.out.println(g);
						groupCount ++;
					}
				} else {
					// when using the Reduced Chem Comp provider
					// there are 3 groups in 2CI1 which are non-standard: SNC, KOR, CIT
					// they are not in the reduced set of standard definitions that will
					// be shipped in the .jar file.

					// if the download chem comp provider is used
					// there will be CIT, which is not a peptide, but
					// should still be counted as a valid HETATOM group...
					if (! ChainImpl.waternames.contains(g.getPDBName())) {
						//System.out.println(cc);
						//System.out.println(g);
						groupCount++;
					}
				}
			}


			ResidueNumber resNum = ResidueNumber.fromString("273");
			resNum.setChainId("A");

			Group g = a.getGroupByPDB(resNum);

			assertEquals("The residue number is not correct", resNum, g.getResidueNumber());

			assertTrue("The group does not have an altLoc ", g.hasAltLoc());

			assertTrue("The nr of altLocs is not 1, but " + g.getAltLocs().size(), g.getAltLocs().size() == 1);

			assertEquals( g.getPDBName(), "KOR");

			Group altLocG = g.getAltLocs().get(0);

			assertEquals(altLocG.getPDBName(),"K1R");

			assertEquals(276,groupCount);


			ResidueNumber resNum2 = ResidueNumber.fromString("265");

			Group g2 = a.getGroupByPDB(resNum2);
			assertTrue(g2.hasAltLoc());


		} catch (Exception e) {
			fail(e.getMessage());
		}

	}

	public void test2W72(){
		try {
			AtomCache cache = new AtomCache();
			Structure s = cache.getStructure("2W72");

			Chain a = s.getChainByPDB("A");

			Group val1 = a.getGroupByPDB(ResidueNumber.fromString("1"));
			Atom ca1 = val1.getAtom(" CA ");
			assertNotNull(ca1);

			Group lys7 = a.getGroupByPDB(ResidueNumber.fromString("7"));
			Atom ca7 = lys7.getAtom(" CA ");			
			assertNotNull(ca7);

			Atom[] caA = StructureTools.getAtomCAArray(a);

			assertEquals(caA.length,141);

		} catch(Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}

	}

	public void test1U7F(){
		try {
			AtomCache cache = new AtomCache();
			Structure s = cache.getStructure("1U7F");

			Chain c = s.getChainByPDB("B");

			Group g = c.getGroupByPDB(ResidueNumber.fromString("314"));
			//System.out.println("== original group ==");
			ensureAllAtomsSameInsCode(g);
			//System.out.println("== alternate group ==");
			for ( Group altGroup : g.getAltLocs() ) {
				ensureAllAtomsSameInsCode(altGroup);	
			}

		} catch(Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}	
	}

	public void test1JXX(){
		try {
			AtomCache cache = new AtomCache();
			Structure structure = cache.getStructure("1JXX");

			Chain chain = structure.getChain(0); // 1JXX example

			Group g = chain.getAtomGroups().get(1); // 1JXX  THR A   2
			ensureAllAtomsSameInsCode(g);
			//System.out.println("== alternate group ==");
			for ( Group altGroup : g.getAltLocs() ) {
				ensureAllAtomsSameInsCode(altGroup);	
			}
			
		} catch(Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}	

	}



	private void ensureAllAtomsSameInsCode(Group g) {

	//	System.out.println(String.format("Group size: %d", g.getAtoms().size()));

		Character defaultAltLoc = null;
		for (Atom atom : g.getAtoms()) {
			if ( defaultAltLoc == null) {
				defaultAltLoc = atom.getAltLoc();
				continue;
			}
		//	System.out.print(atom.toPDB());
			Character altLoc = atom.getAltLoc();

			assertEquals(defaultAltLoc,altLoc);
		}
	}

}
