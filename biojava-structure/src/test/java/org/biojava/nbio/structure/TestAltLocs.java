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
import org.biojava.nbio.structure.io.mmcif.MMCIFFileTools;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;
import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.io.mmcif.chem.ResidueType;
import org.biojava.nbio.structure.io.mmcif.model.AtomSite;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.structure.io.mmcif.model.ChemCompBond;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
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

		assertEquals(1, g.getAltLocs().size());

		assertEquals( g.getPDBName(), "KOR");

		Group altLocG = g.getAltLocs().get(0);

		assertEquals(altLocG.getPDBName(),"K1R");

		assertEquals(275, groupCount);

		// citric acid is now in its own chain

		Chain b = s.getChain("B");
		assertEquals(1, b.getAtomGroups().size());

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

		assertEquals(141, caA.length);

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

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(false);
		StructureIO.setAtomCache(cache);

		Structure s = StructureIO.getStructure("1AAC");

		Chain a = s.getPolyChainByPDB("A");

		Group g = a.getGroupByPDB( ResidueNumber.fromString("27"));
		testCBAtomInMainGroup(g);

		cache = new AtomCache();
		cache.setUseMmCif(true);
		StructureIO.setAtomCache(cache);

		Structure s1 = cache.getStructure("1AAC");
		Chain a1 = s1.getPolyChainByPDB("A");

		Group g1 = a1.getGroupByPDB( ResidueNumber.fromString("27"));

		testCBAtomInMainGroup(g1);

	}

	private void testCBAtomInMainGroup(Group g) {
		// test position of C-beta

		boolean cbInMain = false;

		for (Atom atom : g.getAtoms()) {
			//System.out.print(atom.toPDB());
			if ( atom.getName().equals(StructureTools.CB_ATOM_NAME)){

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

		assertEquals(ca.length, caList.size());

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
					assertEquals(g.size(), altLocGroup.size());
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

		assertEquals(ca.length, caList.size());


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

	@Test
	public void testMmcifConversionPartialAltlocs() throws IOException {
		String mmcifData =
				"data_test\n" +
						"loop_\n" +
						"_atom_site.group_PDB \n" +
						"_atom_site.id \n" +
						"_atom_site.type_symbol \n" +
						"_atom_site.label_atom_id \n" +
						"_atom_site.label_alt_id \n" +
						"_atom_site.label_comp_id \n" +
						"_atom_site.label_asym_id \n" +
						"_atom_site.label_entity_id \n" +
						"_atom_site.label_seq_id \n" +
						"_atom_site.pdbx_PDB_ins_code \n" +
						"_atom_site.Cartn_x \n" +
						"_atom_site.Cartn_y \n" +
						"_atom_site.Cartn_z \n" +
						"_atom_site.occupancy \n" +
						"_atom_site.B_iso_or_equiv \n" +
						"_atom_site.pdbx_formal_charge \n" +
						"_atom_site.auth_seq_id \n" +
						"_atom_site.auth_comp_id \n" +
						"_atom_site.auth_asym_id \n" +
						"_atom_site.auth_atom_id \n" +
						"_atom_site.pdbx_PDB_model_num \n" +
						"ATOM   102 N N   . ARG A 1 13 ? 9.889  23.379 13.115 1.00 6.57  ? 102  ARG A N   1\n" +
						"ATOM   103 C CA  . ARG A 1 13 ? 9.540  23.003 14.482 1.00 7.05  ? 102  ARG A CA  1\n" +
						"ATOM   104 C C   . ARG A 1 13 ? 10.407 23.758 15.489 1.00 6.88  ? 102  ARG A C   1\n" +
						"ATOM   105 O O   . ARG A 1 13 ? 9.915  24.196 16.532 1.00 7.69  ? 102  ARG A O   1\n" +
						"ATOM   106 C CB  . ARG A 1 13 ? 9.706  21.494 14.688 1.00 9.07  ? 102  ARG A CB  1\n" +
						"ATOM   107 C CG  A ARG A 1 13 ? 8.757  20.644 13.854 0.50 14.39 ? 102  ARG A CG  1\n" +
						"ATOM   108 C CG  B ARG A 1 13 ? 8.693  20.645 13.938 0.50 13.58 ? 102  ARG A CG  1\n" +
						"ATOM   109 C CD  A ARG A 1 13 ? 9.109  19.164 13.950 0.50 18.14 ? 102  ARG A CD  1\n" +
						"ATOM   110 C CD  B ARG A 1 13 ? 8.710  19.216 14.456 0.50 16.34 ? 102  ARG A CD  1\n" +
						"ATOM   111 N NE  A ARG A 1 13 ? 8.983  18.644 15.310 0.50 20.72 ? 102  ARG A NE  1\n" +
						"ATOM   112 N NE  B ARG A 1 13 ? 8.315  19.158 15.861 0.50 23.99 ? 102  ARG A NE  1\n" +
						"ATOM   113 C CZ  A ARG A 1 13 ? 7.826  18.445 15.933 0.50 23.45 ? 102  ARG A CZ  1\n" +
						"ATOM   114 C CZ  B ARG A 1 13 ? 8.404  18.072 16.622 0.50 24.56 ? 102  ARG A CZ  1\n" +
						"ATOM   115 N NH1 A ARG A 1 13 ? 6.683  18.718 15.321 0.50 23.60 ? 102  ARG A NH1 1\n" +
						"ATOM   116 N NH1 B ARG A 1 13 ? 8.881  16.942 16.118 0.50 28.42 ? 102  ARG A NH1 1\n" +
						"ATOM   117 N NH2 A ARG A 1 13 ? 7.812  17.972 17.172 0.50 24.80 ? 102  ARG A NH2 1\n" +
						"ATOM   118 N NH2 B ARG A 1 13 ? 8.013  18.115 17.888 0.50 26.52 ? 102  ARG A NH2 1\n";

		SimpleMMcifParser parser = new SimpleMMcifParser();
		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
		parser.addMMcifConsumer(consumer);

		BufferedReader buf = new BufferedReader(new StringReader(mmcifData));
		parser.parse(buf);
		buf.close();

		Structure s = consumer.getStructure();
		Chain c = s.getPolyChains().get(0);
		assertEquals(1, c.getAtomGroups().size());
		Group g = c.getAtomGroup(0);
		assertEquals(11, g.size());

		// there's the main group (. and A) plus the 2 alt loc groups (A and B). The alt locs will contain all the '.' atoms too
		assertEquals(2, g.getAltLocs().size());

		for (Atom a : g.getAtoms()) {
			if (a.getName().equals("C") || a.getName().equals("N") || a.getName().equals("O") || a.getName().equals("CA") || a.getName().equals("CB"))
				assertEquals(' ', a.getAltLoc().charValue());
			else
				assertEquals('A', a.getAltLoc().charValue());
		}

		assertEquals(11, g.getAltLocs().get(0).size());
		for (Atom a : g.getAltLocs().get(0).getAtoms()) {
			if (a.getName().equals("C") || a.getName().equals("N") || a.getName().equals("O") || a.getName().equals("CA") || a.getName().equals("CB"))
				assertEquals(' ', a.getAltLoc().charValue());
			else
				assertEquals('A', a.getAltLoc().charValue());
		}

		assertEquals(11, g.getAltLocs().get(1).size());
		for (Atom a : g.getAltLocs().get(1).getAtoms()) {
			if (a.getName().equals("C") || a.getName().equals("N") || a.getName().equals("O") || a.getName().equals("CA") || a.getName().equals("CB"))
				assertEquals(' ', a.getAltLoc().charValue());
			else
				assertEquals('B', a.getAltLoc().charValue());
		}

		List<AtomSite> atomSites = MMCIFFileTools.convertChainToAtomSites(c, 1, "A", "A");
		assertEquals(17, atomSites.size());

	}

	@Test
	public void testMmcifConversionAllAltlocs() throws IOException {
		String mmcifData =
				"data_test\n" +
						"loop_\n" +
						"_atom_site.group_PDB \n" +
						"_atom_site.id \n" +
						"_atom_site.type_symbol \n" +
						"_atom_site.label_atom_id \n" +
						"_atom_site.label_alt_id \n" +
						"_atom_site.label_comp_id \n" +
						"_atom_site.label_asym_id \n" +
						"_atom_site.label_entity_id \n" +
						"_atom_site.label_seq_id \n" +
						"_atom_site.pdbx_PDB_ins_code \n" +
						"_atom_site.Cartn_x \n" +
						"_atom_site.Cartn_y \n" +
						"_atom_site.Cartn_z \n" +
						"_atom_site.occupancy \n" +
						"_atom_site.B_iso_or_equiv \n" +
						"_atom_site.pdbx_formal_charge \n" +
						"_atom_site.auth_seq_id \n" +
						"_atom_site.auth_comp_id \n" +
						"_atom_site.auth_asym_id \n" +
						"_atom_site.auth_atom_id \n" +
						"_atom_site.pdbx_PDB_model_num \n" +
						"ATOM   204 N N   A PRO A 1 23 ? 15.057 31.425 23.772 0.50 3.09  ? 112  PRO A N   1 \n" +
						"ATOM   205 N N   B PRO A 1 23 ? 14.762 31.778 23.217 0.50 15.25 ? 112  PRO A N   1 \n" +
						"ATOM   206 C CA  A PRO A 1 23 ? 16.391 30.930 23.416 0.50 5.82  ? 112  PRO A CA  1 \n" +
						"ATOM   207 C CA  B PRO A 1 23 ? 16.049 31.406 22.622 0.50 15.44 ? 112  PRO A CA  1 \n" +
						"ATOM   208 C C   A PRO A 1 23 ? 17.360 30.580 24.546 0.50 6.73  ? 112  PRO A C   1 \n" +
						"ATOM   209 C C   B PRO A 1 23 ? 16.971 30.922 23.734 0.50 15.04 ? 112  PRO A C   1 \n" +
						"ATOM   210 O O   A PRO A 1 23 ? 18.566 30.784 24.409 0.50 10.00 ? 112  PRO A O   1 \n" +
						"ATOM   211 O O   B PRO A 1 23 ? 18.076 31.430 23.925 0.50 14.61 ? 112  PRO A O   1 \n" +
						"ATOM   212 C CB  A PRO A 1 23 ? 16.931 32.050 22.542 0.50 8.38  ? 112  PRO A CB  1 \n" +
						"ATOM   213 C CB  B PRO A 1 23 ? 16.519 32.710 21.986 0.50 14.09 ? 112  PRO A CB  1 \n" +
						"ATOM   214 C CG  A PRO A 1 23 ? 16.424 33.256 23.263 0.50 7.59  ? 112  PRO A CG  1 \n" +
						"ATOM   215 C CG  B PRO A 1 23 ? 15.960 33.743 22.908 0.50 15.66 ? 112  PRO A CG  1 \n" +
						"ATOM   216 C CD  A PRO A 1 23 ? 14.980 32.886 23.580 0.50 6.98  ? 112  PRO A CD  1 \n" +
						"ATOM   217 C CD  B PRO A 1 23 ? 14.558 33.235 23.153 0.50 14.91 ? 112  PRO A CD  1 \n";

		SimpleMMcifParser parser = new SimpleMMcifParser();
		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
		parser.addMMcifConsumer(consumer);

		BufferedReader buf = new BufferedReader(new StringReader(mmcifData));
		parser.parse(buf);
		buf.close();

		Structure s = consumer.getStructure();
		Chain c = s.getPolyChains().get(0);
		assertEquals(1, c.getAtomGroups().size());

		Group g = c.getAtomGroup(0);
		assertEquals(7, g.size());

		assertEquals(1, g.getAltLocs().size());

		for (Atom a : g.getAtoms()) {
			assertEquals('A', a.getAltLoc().charValue());
		}
		for (Atom a : g.getAltLocs().get(0).getAtoms()) {
			assertEquals('B', a.getAltLoc().charValue());
		}

		List<AtomSite> atomSites = MMCIFFileTools.convertChainToAtomSites(c, 1, "A", "A");
		assertEquals(14, atomSites.size());

	}

}
