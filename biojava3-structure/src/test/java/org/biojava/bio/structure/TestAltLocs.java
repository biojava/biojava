package org.biojava.bio.structure;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.mmcif.chem.PolymerType;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.biojava3.structure.StructureIO;
import org.junit.Test;
import static org.junit.Assert.*;

public class TestAltLocs {

	@Test
	public void testAltLocParsing() throws StructureException, IOException{ 


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
				if (! g.isWater()) {
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
				if (! g.isWater()) {
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


	}

	@Test
	public void test2W72() throws IOException, StructureException{
		
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


	}

	@Test
	public void test1U7F() throws IOException, StructureException{

		AtomCache cache = new AtomCache();
		Structure s = cache.getStructure("1U7F");

		Chain c = s.getChainByPDB("B");

		Group g = c.getGroupByPDB(ResidueNumber.fromString("314"));
		//System.out.println("== original group ==");
		ensureAllAtomsSameAltCode(g);
		//System.out.println("== alternate group ==");
		for ( Group altGroup : g.getAltLocs() ) {
			ensureAllAtomsSameAltCode(altGroup);	
		}


	}

	@Test
	public void test1JXX() throws IOException, StructureException{

		AtomCache cache = new AtomCache();
		Structure structure = cache.getStructure("1JXX");

		Chain chain = structure.getChain(0); // 1JXX example

		Group g = chain.getAtomGroups().get(1); // 1JXX  THR A   2
		ensureAllAtomsSameAltCode(g);
		//System.out.println("== alternate group ==");
		for ( Group altGroup : g.getAltLocs() ) {
			ensureAllAtomsSameAltCode(altGroup);	
		}


	}

	private void ensureAllAtomsSameAltCode(Group g) {
		
		Character defaultAltLoc = null;
		for (Atom atom : g.getAtoms()) {
		
			if ( defaultAltLoc == null) {
				defaultAltLoc = atom.getAltLoc();
		
				continue;
			}

			Character altLoc = atom.getAltLoc();

			assertEquals(defaultAltLoc,altLoc);
		}
	}

	@Test
	public void test1AAC() throws IOException, StructureException{ 

		Structure s = StructureIO.getStructure("1AAC");

		Chain a = s.getChainByPDB("A");

		Group g = a.getGroupByPDB( ResidueNumber.fromString("27"));
		testCBAtomInMainGroup(g);

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);

		Structure s1 = cache.getStructure("1AAC");
		Chain a1 = s1.getChainByPDB("A");

		Group g1 = a1.getGroupByPDB( ResidueNumber.fromString("27"));

		testCBAtomInMainGroup(g1);



		//			int pos = 0;
		//			for (Group alt: g.getAltLocs()) {
		//				pos++;
		//				System.out.println("altLoc: " + pos + " " + alt);
		//				for (Atom atom : alt.getAtoms()) {
		//					System.out.print(atom.toPDB());
		//				}
		//			}




	}

	private void testCBAtomInMainGroup(Group g) {
		// test position of C-beta

		boolean cbInMain = false;

		for (Atom atom : g.getAtoms()) {
			//System.out.print(atom.toPDB());
			if ( atom.getFullName().equals(StructureTools.caAtomName)){

				cbInMain = true;
				break;
			}
		}

		assertTrue("Did not find C beta atom in main group",cbInMain);

	}

	@Test
	public void test3PIUpdb() throws IOException, StructureException{ 

		AtomCache cache = new AtomCache("/tmp",false);
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);

		Structure structure = StructureIO.getStructure("3PIU");

		assertNotNull(structure);

		Atom[] ca = StructureTools.getAtomCAArray(structure);

		//System.out.println(structure.getPdbId() + " has # CA atoms: " + ca.length);

		List<Atom> caList = new ArrayList<Atom>();
		for ( Chain c: structure.getChains()){
			for (Group g: c.getAtomGroups()){

				for (Group altLocGroup:g.getAltLocs()) {
					ensureAllAtomsSameAltCode(altLocGroup);						
				}

				List<Atom> atoms = g.getAtoms();
				boolean caInMain = false;
				for (Atom a: atoms){

					if ( a.getFullName().equals(StructureTools.caAtomName)) {
						caList.add(a);
						caInMain = true;
						break;

					}


				}
				if (! caInMain && g.hasAtom(StructureTools.caAtomName)){
					// g.hasAtom checks altLocs
					fail("CA is not in main group, but in altLoc");
				}

			}
		}

		assertTrue(ca.length == caList.size());


	}
	
	@Test
	public void test3PIUmmcif() throws IOException, StructureException{ 

		AtomCache cache = new AtomCache("/tmp",false);
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(true);

		Structure structure = StructureIO.getStructure("3PIU");

		assertNotNull(structure);

		Atom[] ca = StructureTools.getAtomCAArray(structure);

		//System.out.println(structure.getPdbId() + " has # CA atoms: " + ca.length);

		List<Atom> caList = new ArrayList<Atom>();
		for ( Chain c: structure.getChains()){
			for (Group g: c.getAtomGroups()){

				for (Group altLocGroup:g.getAltLocs()) {
					ensureAllAtomsSameAltCode(altLocGroup);						
				}

				List<Atom> atoms = g.getAtoms();
				boolean caInMain = false;
				for (Atom a: atoms){

					if ( a.getFullName().equals(StructureTools.caAtomName)) {
						caList.add(a);
						caInMain = true;
						break;

					}


				}
				if (! caInMain && g.hasAtom(StructureTools.caAtomName)){
					// g.hasAtom checks altLocs
					fail("CA is not in main group, but in altLoc");
				}

			}
		}

		assertTrue(ca.length == caList.size());


	}

}
