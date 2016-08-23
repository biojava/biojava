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
package org.biojava.nbio.structure;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.io.mmcif.chem.ResidueType;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.structure.io.mmcif.model.ChemCompBond;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

public class TestAltLocs {

	@Test
	public void testAltLocParsing() throws StructureException, IOException{

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		Structure s = cache.getStructure("2CI1");

		Chain a = s.getPolyChainByPDB("A");

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
		resNum.setChainName("A");

		Group g = a.getGroupByPDB(resNum);

		assertEquals("The residue number is not correct", resNum, g.getResidueNumber());

		assertTrue("The group does not have an altLoc ", g.hasAltLoc());

		assertTrue("The nr of altLocs is not 1, but " + g.getAltLocs().size(), g.getAltLocs().size() == 1);

		assertEquals( g.getPDBName(), "KOR");

		Group altLocG = g.getAltLocs().get(0);

		assertEquals(altLocG.getPDBName(),"K1R");

		assertEquals(275,groupCount);

		// citric acid is now in its own chain

		Chain b = s.getChain("B");
		assertTrue(b.getAtomGroups().size() == 1);


		ResidueNumber resNum2 = ResidueNumber.fromString("265");

		Group g2 = a.getGroupByPDB(resNum2);
		assertTrue(g2.hasAltLoc());


	}

	@Test
	public void test2W72() throws IOException, StructureException{

		AtomCache cache = new AtomCache();
		Structure s = cache.getStructure("2W72");

		Chain a = s.getPolyChainByPDB("A");

		Group val1 = a.getGroupByPDB(ResidueNumber.fromString("1"));
		Atom ca1 = val1.getAtom("CA");
		assertNotNull(ca1);

		Group lys7 = a.getGroupByPDB(ResidueNumber.fromString("7"));
		Atom ca7 = lys7.getAtom("CA");
		assertNotNull(ca7);

		Atom[] caA = StructureTools.getRepresentativeAtomArray(a);

		assertEquals(caA.length,141);


	}

	@Test
	public void test1U7F() throws IOException, StructureException{

		AtomCache cache = new AtomCache();
		Structure s = cache.getStructure("1U7F");

		Chain c = s.getPolyChainByPDB("B");

		Group g = c.getGroupByPDB(ResidueNumber.fromString("314"));
		//System.out.println("== original group ==");
		ensureAllAtomsSameAltCode(g, g);
		//System.out.println("== alternate group ==");
		for ( Group altGroup : g.getAltLocs() ) {
			ensureAllAtomsSameAltCode(altGroup, g);
		}


	}

	@Test
	public void test1JXX() throws IOException, StructureException{

		AtomCache cache = new AtomCache();
		Structure structure = cache.getStructure("1JXX");

		Chain chain = structure.getChainByIndex(0); // 1JXX example

		Group g = chain.getAtomGroups().get(1); // 1JXX  THR A   2
		ensureAllAtomsSameAltCode(g, g);
		//System.out.println("== alternate group ==");
		for ( Group altGroup : g.getAltLocs() ) {
			ensureAllAtomsSameAltCode(altGroup, g);
		}


	}

	
	/**
	 * Test to check that all atoms have the same alt code (unless they're in the main group)
	 * @param groupInputAltLocGroup The input alt loc group
	 */
	private void ensureAllAtomsSameAltCode(Group groupInputAltLocGroup, Group inputMainGroup) {

		// If they're the exact same group just return
		if (groupInputAltLocGroup == inputMainGroup) {
			return;
		}
		
		// Check that the atom group is the same size as the alt loc group (as long as it's not a case of microheterogenity
		if (groupInputAltLocGroup.getPDBName().equals(inputMainGroup.getPDBName())){
			assertEquals(groupInputAltLocGroup.getAtoms().size(), inputMainGroup.getAtoms().size());
		}
		Character defaultAltLoc = null;
		for (Atom atom : groupInputAltLocGroup.getAtoms()) {
			
			// If this is in the original atom group just carry on
			if (inputMainGroup.getAtoms().contains(atom)) {
				continue;
			}
			
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

		Chain a = s.getPolyChainByPDB("A");

		Group g = a.getGroupByPDB( ResidueNumber.fromString("27"));
		testCBAtomInMainGroup(g);

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);

		Structure s1 = cache.getStructure("1AAC");
		Chain a1 = s1.getPolyChainByPDB("A");

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
			if ( atom.getName().equals(StructureTools.CA_ATOM_NAME)){

				cbInMain = true;
				break;
			}
		}

		assertTrue("Did not find C beta atom in main group",cbInMain);

	}

	@Test
	public void test3PIUpdb() throws IOException, StructureException{

		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);


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
					ensureAllAtomsSameAltCode(altLocGroup, g);

					for (Atom a:altLocGroup.getAtoms()) {
						assertNotNull(a.getGroup());
						assertNotNull(a.getGroup().getChain());
					}
				}

				List<Atom> atoms = g.getAtoms();
				boolean caInMain = false;
				for (Atom a: atoms){
					assertNotNull(a.getGroup());
					assertNotNull(a.getGroup().getChain());

					if ( a.getName().equals(StructureTools.CA_ATOM_NAME)) {
						caList.add(a);
						caInMain = true;
						break;

					}


				}
				if (! caInMain && g.hasAtom(StructureTools.CA_ATOM_NAME)){
					// g.hasAtom checks altLocs
					fail("CA is not in main group, but in altLoc");
				}

			}
		}

		assertTrue(ca.length == caList.size());

	}


	/**
	 * A test that all alternate location groups have the same number of atoms as the main group
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void testAllAltLocsSameAtomsMainGroup() throws IOException, StructureException {
		doTestAllAltLocsSamAtomsMainGroup("3nu4");
		doTestAllAltLocsSamAtomsMainGroup("3nvd");
		doTestAllAltLocsSamAtomsMainGroup("4txr");
		doTestAllAltLocsSamAtomsMainGroup("3nvd");
		doTestAllAltLocsSamAtomsMainGroup("4cup");
	}
	
	/**
	 * Actually perform the test to see all alt locs are the same size as the main group
	 * @throws StructureException 
	 * @throws IOException  
	 * 
	 */
	private void doTestAllAltLocsSamAtomsMainGroup(String pdbId) throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		params.setCreateAtomBonds(true);
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
		Structure structure = StructureIO.getStructure(pdbId);
		// Loop through the atoms
		for ( Chain c: structure.getChains()){
			for (Group g: c.getAtomGroups()){

				for (Group altLocGroup:g.getAltLocs()) {
					assertEquals(altLocGroup.size(), g.size());
				}
			}
		}
	}
	
	/**
	 * A test that adding bonds to atoms between groups - doesn't change the size of the groups
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void testAddBondsDoesntChangeGroups() throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		params.setCreateAtomBonds(true);
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
		Structure structure = StructureIO.getStructure("4CUP");
		// Loop through and find 
		for (Chain chain : structure.getChains()) {
			List<Group> groups = chain.getAtomGroups();

			for (Group mainGroup : groups) {
				// atoms with no residue number don't have atom information
				if (mainGroup.getResidueNumber() == null) {
					continue;
				}
				if (mainGroup.getAltLocs().isEmpty()) {
					continue;
				}
				int oldSize = mainGroup.size();
				// Now add support for altLocGroup
				List<Atom> atomsList = new ArrayList<>(mainGroup.getAtoms());
				for(Group altLocOne: mainGroup.getAltLocs()){
					for(Atom atomAltLocOne: altLocOne.getAtoms()){
						atomsList.add(atomAltLocOne);
					}
				}
				// Get the chem copm
				ChemComp aminoChemComp = ChemCompGroupFactory.getChemComp(mainGroup
						.getPDBName());
				// Now iterate through this list
				for(Atom atomA : atomsList){

					for (ChemCompBond chemCompBond : aminoChemComp.getBonds()) {

						//
						if(chemCompBond.getAtom_id_1().equals(atomA.getName())){
							// Get the other atom in the group
							for(Atom atomB : atomsList) { 
								if(chemCompBond.getAtom_id_2().equals(atomB.getName())){
									int bondOrder = chemCompBond.getNumericalBondOrder();
									new BondImpl(atomA, atomB, bondOrder);
								}
							}
						}
					}
				}
				assertEquals(oldSize, mainGroup.size());
			}
		}
	}

	/**
	 * A test to see that alternate location bonds are being generated
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void test4CUPBonds() throws IOException, StructureException{

		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		params.setCreateAtomBonds(true);
		cache.setFileParsingParams(params);


		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(false);

		Structure structure = StructureIO.getStructure("4CUP");

		assertNotNull(structure);

		Atom[] ca = StructureTools.getAtomCAArray(structure);


		List<Atom> caList = new ArrayList<Atom>();
		for ( Chain c: structure.getChains()){
			for (Group g: c.getAtomGroups()){

				for (Group altLocGroup:g.getAltLocs()) {
					ensureAllAtomsSameAltCode(altLocGroup, g);
					for (Atom a:altLocGroup.getAtoms()) {
						// Check the atomsall have bonds 
						assertNotEquals(a.getBonds(),null);
						assertNotEquals(a.getBonds().size(),0);

					}
				}

				List<Atom> atoms = g.getAtoms();
				boolean caInMain = false;
				for (Atom a: atoms){
					assertNotNull(a.getGroup());
					assertNotNull(a.getGroup().getChain());

					if ( a.getName().equals(StructureTools.CA_ATOM_NAME)) {
						caList.add(a);
						caInMain = true;
						break;

					}


				}
				if (! caInMain && g.hasAtom(StructureTools.CA_ATOM_NAME)){
					// g.hasAtom checks altLocs
					fail("CA is not in main group, but in altLoc");
				}

			}
		}

		assertTrue(ca.length == caList.size());


	}

	@Test
	public void test3PIUmmcif() throws IOException, StructureException{

		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);

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
					ensureAllAtomsSameAltCode(altLocGroup, g);

					for (Atom a:altLocGroup.getAtoms()) {
						assertNotNull(a.getGroup());
						assertNotNull(a.getGroup().getChain());
					}
				}

				List<Atom> atoms = g.getAtoms();
				boolean caInMain = false;
				for (Atom a: atoms){

					assertNotNull(a.getGroup());
					assertNotNull(a.getGroup().getChain());

					if ( a.getName().equals(StructureTools.CA_ATOM_NAME)) {
						caList.add(a);
						caInMain = true;
						break;

					}


				}
				if (! caInMain && g.hasAtom(StructureTools.CA_ATOM_NAME)){
					// g.hasAtom checks altLocs
					fail("CA is not in main group, but in altLoc");
				}

			}
		}

		assertEquals(ca.length, caList.size());


	}

	@Test
	public void test3U7Tmmcif() throws IOException, StructureException{

		// this test intends to check that the mmCIF parser doesn't read twice residues 22 and 25
		// which are annotated as "microheterogeneity" in the SEQRES (entity_poly_seq), see #160

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);

		Structure structure = StructureIO.getStructure("3U7T");

		assertNotNull(structure);

		Atom[] ca = StructureTools.getRepresentativeAtomArray(structure);

		//System.out.println(structure.getPdbId() + " has # CA atoms: " + ca.length);

		List<Atom> caList = new ArrayList<Atom>();
		for ( Chain c: structure.getChains()){
			// notice here we test the seqresgroups, because we want to check if microheterogeinity is treated correctly
			for (Group g: c.getSeqResGroups()){

				for (Group altLocGroup:g.getAltLocs()) {
					ensureAllAtomsSameAltCode(altLocGroup, g);
				}

				List<Atom> atoms = g.getAtoms();
				boolean caInMain = false;
				for (Atom a: atoms){

					if ( a.getName().equals(StructureTools.CA_ATOM_NAME)) {
						caList.add(a);
						caInMain = true;
						break;

					}


				}
				if (! caInMain && g.hasAtom(StructureTools.CA_ATOM_NAME)){
					// g.hasAtom checks altLocs
					fail("CA is not in main group, but in altLoc");
				}

			}
		}

		assertEquals(ca.length, caList.size());


	}

}
