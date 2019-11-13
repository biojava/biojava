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
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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

	/**
	 * Test that intra-residue bonds between alt locs link atoms with same altloc codes
	 * https://github.com/rcsb/mmtf/issues/44
	 */
	@Test
	public void testIntraResidueBondsBetweenAltlocs() throws IOException {
		// from 5MOO
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
						"ATOM   1405 N  N    A MET A 1 86  ? 10.748  -17.610 -6.975  0.47 16.12 ? 104 MET A N    1 \n" +
						"ATOM   1406 N  N    B MET A 1 86  ? 10.802  -17.694 -6.986  0.53 17.92 ? 104 MET A N    1 \n" +
						"ATOM   1407 C  CA   A MET A 1 86  ? 11.189  -17.392 -5.610  0.47 15.78 ? 104 MET A CA   1 \n" +
						"ATOM   1408 C  CA   B MET A 1 86  ? 11.033  -17.368 -5.587  0.53 18.29 ? 104 MET A CA   1 \n" +
						"ATOM   1409 C  C    A MET A 1 86  ? 10.952  -18.663 -4.810  0.47 15.91 ? 104 MET A C    1 \n" +
						"ATOM   1410 C  C    B MET A 1 86  ? 10.882  -18.643 -4.767  0.53 17.40 ? 104 MET A C    1 \n" +
						"ATOM   1411 O  O    A MET A 1 86  ? 10.120  -19.504 -5.154  0.47 18.21 ? 104 MET A O    1 \n" +
						"ATOM   1412 O  O    B MET A 1 86  ? 10.018  -19.474 -5.052  0.53 20.02 ? 104 MET A O    1 \n" +
						"ATOM   1413 C  CB   A MET A 1 86  ? 10.477  -16.204 -4.933  0.47 17.14 ? 104 MET A CB   1 \n" +
						"ATOM   1414 C  CB   B MET A 1 86  ? 10.001  -16.336 -5.111  0.53 18.92 ? 104 MET A CB   1 \n" +
						"ATOM   1415 C  CG   A MET A 1 86  ? 9.019   -16.476 -4.619  0.47 20.01 ? 104 MET A CG   1 \n" +
						"ATOM   1416 C  CG   B MET A 1 86  ? 10.030  -16.038 -3.634  0.53 19.12 ? 104 MET A CG   1 \n" +
						"ATOM   1417 S  SD   A MET A 1 86  ? 8.207   -15.088 -3.838  0.47 22.06 ? 104 MET A SD   1 \n" +
						"ATOM   1418 S  SD   B MET A 1 86  ? 8.874   -14.724 -3.205  0.53 20.16 ? 104 MET A SD   1 \n" +
						"ATOM   1419 C  CE   A MET A 1 86  ? 9.151   -14.973 -2.340  0.47 25.15 ? 104 MET A CE   1 \n" +
						"ATOM   1420 C  CE   B MET A 1 86  ? 7.269   -15.536 -3.380  0.53 20.38 ? 104 MET A CE   1 \n" +
						"ATOM   1421 H  H    A MET A 1 86  ? 9.931   -18.207 -7.055  0.47 15.58 ? 104 MET A H    1 \n" +
						"ATOM   1422 H  H    B MET A 1 86  ? 10.144  -18.461 -7.109  0.53 18.91 ? 104 MET A H    1 \n" +
						"ATOM   1423 H  HA   A MET A 1 86  ? 12.256  -17.182 -5.644  0.47 15.14 ? 104 MET A HA   1 \n" +
						"ATOM   1424 H  HA   B MET A 1 86  ? 12.033  -16.953 -5.465  0.53 19.55 ? 104 MET A HA   1 \n" +
						"ATOM   1425 H  HB2  A MET A 1 86  ? 10.986  -15.920 -4.008  0.47 17.68 ? 104 MET A HB2  1 \n" +
						"ATOM   1426 H  HB3  A MET A 1 86  ? 10.484  -15.364 -5.622  0.47 17.68 ? 104 MET A HB3  1 \n" +
						"ATOM   1427 H  HB3  B MET A 1 86  ? 9.001   -16.676 -5.398  0.53 20.49 ? 104 MET A HB3  1 \n" +
						"ATOM   1428 H  HG2  A MET A 1 86  ? 8.490   -16.704 -5.546  0.47 20.93 ? 104 MET A HG2  1 \n" +
						"ATOM   1429 H  HG3  A MET A 1 86  ? 8.956   -17.315 -3.927  0.47 20.93 ? 104 MET A HG3  1 \n" +
						"ATOM   1430 H  HE2  A MET A 1 86  ? 9.861   -14.153 -2.440  0.47 27.31 ? 104 MET A HE2  1 \n" +
						"ATOM   1431 H  HE2  B MET A 1 86  ? 7.346   -16.554 -2.998  0.53 23.03 ? 104 MET A HE2  1 \n" +
						"ATOM   1432 H  HE3  B MET A 1 86  ? 6.996   -15.566 -4.437  0.53 23.03 ? 104 MET A HE3  1 ";

		SimpleMMcifParser parser = new SimpleMMcifParser();
		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
		parser.addMMcifConsumer(consumer);

		FileParsingParameters params = new FileParsingParameters();
		params.setCreateAtomBonds(true);
		consumer.setFileParsingParameters(params);

		BufferedReader buf = new BufferedReader(new StringReader(mmcifData));
		parser.parse(buf);
		buf.close();

		Structure s = consumer.getStructure();
		Chain c = s.getPolyChains().get(0);
		assertEquals(1, c.getAtomGroups().size());

		Group g = c.getAtomGroup(0);

		assertEquals(1, g.getAltLocs().size());

		boolean foundCEHE3bond = false;
		for (Atom a : g.getAtoms()) {
			for (Bond b : a.getBonds()) {
//				if (b.getAtomA().getAltLoc() != b.getAtomB().getAltLoc()) {
//					System.out.println(
//							b.getAtomA().toString() + ": '" + b.getAtomA().getAltLoc() + "' --- " +
//							b.getAtomB().toString() + ": '" + b.getAtomB().getAltLoc() + "'");
//				}
				// no bonds between atoms with different alt locs
				assertEquals(b.getAtomA().toString() + " --- " + b.getAtomB().toString(),
						b.getAtomA().getAltLoc(), b.getAtomB().getAltLoc());

				// a bond should exist between CE and HE3 but only for altloc=B
				if ((b.getAtomA().getName().equals("CE") && b.getAtomB().getName().equals("HE3")) ||
						(b.getAtomA().getName().equals("HE3") && b.getAtomB().getName().equals("CE")) ) {
					foundCEHE3bond = true;
				}
			}
		}

		// there should be a bond between CE and HE3 but only for altloc=B
		assertTrue(foundCEHE3bond);

	}

	/**
	 * Test that inter-residue bonds between alt locs link atoms with same altloc codes or default alt loc to all alt locs
	 * https://github.com/rcsb/mmtf/issues/44
	 */
	@Test
	public void testInterResidueBondsBetweenAltlocs() throws IOException {
		//  from 5MOO
		String mmcifData =
				"data_test\n" +
						"# \n" +
						"loop_\n" +
						"_entity.id \n" +
						"_entity.type \n" +
						"_entity.src_method \n" +
						"_entity.pdbx_description \n" +
						"_entity.formula_weight \n" +
						"_entity.pdbx_number_of_molecules \n" +
						"_entity.pdbx_ec \n" +
						"_entity.pdbx_mutation \n" +
						"_entity.pdbx_fragment \n" +
						"_entity.details \n" +
						"1 polymer     nat 'Cationic trypsin' 23324.287 1   3.4.21.4 ? ? ? \n" +
						"# \n" +
						"loop_\n" +
						"_entity_poly_seq.entity_id \n" +
						"_entity_poly_seq.num \n" +
						"_entity_poly_seq.mon_id \n" +
						"_entity_poly_seq.hetero \n" +
						"1 1  ILE n \n" +
						"1 2  MET n \n" +
						"# \n" +
						"loop_\n" +
						"_struct_asym.id \n" +
						"_struct_asym.pdbx_blank_PDB_chainid_flag \n" +
						"_struct_asym.pdbx_modified \n" +
						"_struct_asym.entity_id \n" +
						"_struct_asym.details \n" +
						"A N N 1 ? \n" +
						"# \n" +
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
						"ATOM   1385 N  N    . ILE A 1  1  ? 10.900  -16.328 -10.274 1.00 17.47 ? 103 ILE A N    1 \n" +
						"ATOM   1386 C  CA   . ILE A 1  1  ? 10.885  -17.487 -9.388  1.00 17.76 ? 103 ILE A CA   1 \n" +
						"ATOM   1387 C  C    . ILE A 1  1  ? 11.374  -17.058 -8.011  1.00 17.35 ? 103 ILE A C    1 \n" +
						"ATOM   1388 O  O    . ILE A 1  1  ? 12.265  -16.211 -7.883  1.00 18.51 ? 103 ILE A O    1 \n" +
						"ATOM   1389 C  CB   . ILE A 1  1  ? 11.721  -18.644 -9.986  1.00 18.19 ? 103 ILE A CB   1 \n" +
						"ATOM   1390 C  CG1  . ILE A 1  1  ? 11.610  -19.916 -9.144  1.00 19.64 ? 103 ILE A CG1  1 \n" +
						"ATOM   1391 C  CG2  . ILE A 1  1  ? 13.177  -18.246 -10.209 1.00 19.73 ? 103 ILE A CG2  1 \n" +
						"ATOM   1392 C  CD1  . ILE A 1  1  ? 12.217  -21.162 -9.820  1.00 22.94 ? 103 ILE A CD1  1 \n" +
						"ATOM   1393 H  H    A ILE A 1  1  ? 11.598  -15.614 -10.041 1.00 17.71 ? 103 ILE A H    1 \n" +
						"ATOM   1394 D  D    B ILE A 1  1  ? 11.598  -15.614 -10.041 0.00 17.71 ? 103 ILE A D    1 \n" +
						"ATOM   1395 H  HA   . ILE A 1  1  ? 9.856   -17.843 -9.277  1.00 17.70 ? 103 ILE A HA   1 \n" +
						"ATOM   1396 H  HB   . ILE A 1  1  ? 11.300  -18.886 -10.957 1.00 18.93 ? 103 ILE A HB   1 \n" +
						"ATOM   1397 H  HG12 . ILE A 1  1  ? 12.149  -19.788 -8.209  1.00 20.93 ? 103 ILE A HG12 1 \n" +
						"ATOM   1398 H  HG13 . ILE A 1  1  ? 10.563  -20.127 -8.939  1.00 20.93 ? 103 ILE A HG13 1 \n" +
						"ATOM   1399 H  HG21 . ILE A 1  1  ? 13.669  -19.035 -10.776 1.00 20.97 ? 103 ILE A HG21 1 \n" +
						"ATOM   1400 H  HG22 . ILE A 1  1  ? 13.235  -17.312 -10.767 1.00 20.97 ? 103 ILE A HG22 1 \n" +
						"ATOM   1401 H  HG23 . ILE A 1  1  ? 13.683  -18.144 -9.251  1.00 20.97 ? 103 ILE A HG23 1 \n" +
						"ATOM   1402 H  HD11 . ILE A 1  1  ? 13.299  -21.078 -9.905  1.00 24.96 ? 103 ILE A HD11 1 \n" +
						"ATOM   1403 H  HD12 . ILE A 1  1  ? 11.967  -22.036 -9.223  1.00 24.96 ? 103 ILE A HD12 1 \n" +
						"ATOM   1404 H  HD13 . ILE A 1  1  ? 11.779  -21.281 -10.808 1.00 24.96 ? 103 ILE A HD13 1 \n" +
						"ATOM   1405 N  N    A MET A 1  2  ? 10.748  -17.610 -6.975  0.47 16.12 ? 104 MET A N    1 \n" +
						"ATOM   1406 N  N    B MET A 1  2  ? 10.802  -17.694 -6.986  0.53 17.92 ? 104 MET A N    1 \n" +
						"ATOM   1407 C  CA   A MET A 1  2  ? 11.189  -17.392 -5.610  0.47 15.78 ? 104 MET A CA   1 \n" +
						"ATOM   1408 C  CA   B MET A 1  2  ? 11.033  -17.368 -5.587  0.53 18.29 ? 104 MET A CA   1 \n" +
						"ATOM   1409 C  C    A MET A 1  2  ? 10.952  -18.663 -4.810  0.47 15.91 ? 104 MET A C    1 \n" +
						"ATOM   1410 C  C    B MET A 1  2  ? 10.882  -18.643 -4.767  0.53 17.40 ? 104 MET A C    1 \n" +
						"ATOM   1411 O  O    A MET A 1  2  ? 10.120  -19.504 -5.154  0.47 18.21 ? 104 MET A O    1 \n" +
						"ATOM   1412 O  O    B MET A 1  2  ? 10.018  -19.474 -5.052  0.53 20.02 ? 104 MET A O    1 \n" +
						"ATOM   1413 C  CB   A MET A 1  2  ? 10.477  -16.204 -4.933  0.47 17.14 ? 104 MET A CB   1 \n" +
						"ATOM   1414 C  CB   B MET A 1  2  ? 10.001  -16.336 -5.111  0.53 18.92 ? 104 MET A CB   1 \n" +
						"ATOM   1415 C  CG   A MET A 1  2  ? 9.019   -16.476 -4.619  0.47 20.01 ? 104 MET A CG   1 \n" +
						"ATOM   1416 C  CG   B MET A 1  2  ? 10.030  -16.038 -3.634  0.53 19.12 ? 104 MET A CG   1 \n" +
						"ATOM   1417 S  SD   A MET A 1  2  ? 8.207   -15.088 -3.838  0.47 22.06 ? 104 MET A SD   1 \n" +
						"ATOM   1418 S  SD   B MET A 1  2  ? 8.874   -14.724 -3.205  0.53 20.16 ? 104 MET A SD   1 \n" +
						"ATOM   1419 C  CE   A MET A 1  2  ? 9.151   -14.973 -2.340  0.47 25.15 ? 104 MET A CE   1 \n" +
						"ATOM   1420 C  CE   B MET A 1  2  ? 7.269   -15.536 -3.380  0.53 20.38 ? 104 MET A CE   1 \n" +
						"ATOM   1421 H  H    A MET A 1  2  ? 9.931   -18.207 -7.055  0.47 15.58 ? 104 MET A H    1 \n" +
						"ATOM   1422 H  H    B MET A 1  2  ? 10.144  -18.461 -7.109  0.53 18.91 ? 104 MET A H    1 \n" +
						"ATOM   1423 H  HA   A MET A 1  2  ? 12.256  -17.182 -5.644  0.47 15.14 ? 104 MET A HA   1 \n" +
						"ATOM   1424 H  HA   B MET A 1  2  ? 12.033  -16.953 -5.465  0.53 19.55 ? 104 MET A HA   1 \n" +
						"ATOM   1425 H  HB2  A MET A 1  2  ? 10.986  -15.920 -4.008  0.47 17.68 ? 104 MET A HB2  1 \n" +
						"ATOM   1426 H  HB3  A MET A 1  2  ? 10.484  -15.364 -5.622  0.47 17.68 ? 104 MET A HB3  1 \n" +
						"ATOM   1427 H  HB3  B MET A 1  2  ? 9.001   -16.676 -5.398  0.53 20.49 ? 104 MET A HB3  1 \n" +
						"ATOM   1428 H  HG2  A MET A 1  2  ? 8.490   -16.704 -5.546  0.47 20.93 ? 104 MET A HG2  1 \n" +
						"ATOM   1429 H  HG3  A MET A 1  2  ? 8.956   -17.315 -3.927  0.47 20.93 ? 104 MET A HG3  1 \n" +
						"ATOM   1430 H  HE2  A MET A 1  2  ? 9.861   -14.153 -2.440  0.47 27.31 ? 104 MET A HE2  1 \n" +
						"ATOM   1431 H  HE2  B MET A 1  2  ? 7.346   -16.554 -2.998  0.53 23.03 ? 104 MET A HE2  1 \n" +
						"ATOM   1432 H  HE3  B MET A 1  2  ? 6.996   -15.566 -4.437  0.53 23.03 ? 104 MET A HE3  1 ";

		SimpleMMcifParser parser = new SimpleMMcifParser();
		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
		parser.addMMcifConsumer(consumer);

		FileParsingParameters params = new FileParsingParameters();
		params.setCreateAtomBonds(true);
		consumer.setFileParsingParameters(params);

		BufferedReader buf = new BufferedReader(new StringReader(mmcifData));
		parser.parse(buf);
		buf.close();

		Structure s = consumer.getStructure();
		Chain c = s.getPolyChains().get(0);
		assertEquals(2, c.getAtomGroups().size());

		// inter residue bonds and alt locs
		// ILE-C (.) must be linked to both MET-N (A and B alt locs)
		Group g1 = c.getAtomGroup(0);

		Atom catom = g1.getAtom("C");
		List<Bond> bonds = new ArrayList<>();
		for (Bond b : catom.getBonds()) {
			if (b.getAtomA().getName().equals("N") || b.getAtomB().getName().equals("N")) {
				bonds.add(b);
			}
		}

		assertEquals(2, bonds.size());

		Set<Character> seenAltLocs = new HashSet<>();
		for (Bond b : bonds) {
			Atom aAtom = b.getAtomA();
			Atom bAtom = b.getAtomB();
			Atom nAtom;
			if (aAtom.getName().equals("N")) {
				nAtom = aAtom;
			} else {
				nAtom = bAtom;
			}
			seenAltLocs.add(nAtom.getAltLoc());
		}
		// 2 distinct N atoms: alt loc A and B
		assertEquals(2, seenAltLocs.size());
		assertTrue(seenAltLocs.contains('A'));
		assertTrue(seenAltLocs.contains('B'));

	}


}
