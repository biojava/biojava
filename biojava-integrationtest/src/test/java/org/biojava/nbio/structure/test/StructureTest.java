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
package org.biojava.nbio.structure.test;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.geometry.SuperPositions;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.biojava.nbio.structure.io.SSBondImpl;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.Set;

import javax.vecmath.Matrix4d;

import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 *
 * @author Andreas Prlic
 * @since 1.5
 */
public class StructureTest {

	private static Structure structure;

	@BeforeClass
	public static void setUp() throws IOException {
		InputStream inStream = StructureTest.class.getResourceAsStream("/5pti_old.pdb");
		assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		params.setCreateAtomBonds(true);
		pdbpars.setFileParsingParameters(params);

		structure = pdbpars.parsePDBFile(inStream) ;

		assertNotNull(structure);

		assertEquals("structure does not contain one chain ", 1 ,structure.size());
	}

	@Test
	public void testSeqResParsing() {

		// System.out.println(structure);
		List<Chain> chains = structure.getChains(0);
		// since biojava 5.0, we have 4 chains here: 1 protein, 2 non-poly (ligands), 1 water
		assertEquals(" nr of found chains not correct!",4,chains.size());
		Chain c = chains.get(0);
		//System.out.println(c);
		List<Group> seqResGroups = c.getSeqResGroups();
		assertEquals("nr of SEQRES groups not correct!",58,seqResGroups.size());

		List<Group> atomGroups = c.getAtomGroups();

		Group g3 = seqResGroups.get(2);
		int indexAtom = atomGroups.indexOf(g3);
		//System.out.println(" index in atomlist: " + indexAtom);
		assertEquals("the SEQRES group can not be found in the ATOM list",2,indexAtom);


		Group g5 = atomGroups.get(5);
		assertEquals("The ATOM group can not be fond in the SEQRES list", 5,seqResGroups.indexOf(g5));

		Chain c2 = chains.get(1);
		List<Group>atomGroups2 = c2.getAtomGroups();
		Group g58 = atomGroups2.get(0);
		assertEquals("The group is not PO4","PO4", g58.getPDBName());
		assertEquals("The group P04 should not be in the SEQRES list", -1 , seqResGroups.indexOf(g58));

	}


	/**
	 * Tests if a PDB file can be parsed
	 * @throws Exception
	 */
	@Test
	public void testReadPDBFile() throws Exception {

		assertEquals("pdb code not set!","5PTI",structure.getPDBCode());

		// since biojava 5.0, we have 4 chains here: 1 protein, 2 non-poly (ligands), 1 water
		
		Chain c = structure.getChainByIndex(0);
		assertEquals("did not find the expected 58 amino acids!",58,c.getAtomGroups(GroupType.AMINOACID).size());

		assertEquals(0 , c.getAtomGroups(GroupType.HETATM).size());

		Chain c4 = structure.getChainByIndex(3);

		// The fourth chain in the file contains 63 molecules of deutarated
		assertEquals(63, c4.getAtomGroups(GroupType.HETATM).size());
		assertEquals(0, c4.getAtomGroups(GroupType.NUCLEOTIDE).size());

		List<EntityInfo> compounds= structure.getEntityInfos();

		assertEquals(4, compounds.size());
		EntityInfo mol = compounds.get(0);
		assertTrue(mol.getDescription().startsWith("TRYPSIN INHIBITOR"));
	}

	@Test
	public void testSSBondParsing() throws Exception {
		assertNotNull(structure);

		List<Bond> ssbonds = structure.getSSBonds();
		assertEquals("did not find the correct nr of SSBonds ",3,ssbonds.size());

		String pdb1 = "SSBOND   1 CYS A    5    CYS A   55";
		String pdb2 = "SSBOND   2 CYS A   14    CYS A   38";

		Bond bond1 = ssbonds.get(0);
		assertDisulfideBond("A", "A", 5, 55, bond1);

		Bond bond2 = ssbonds.get(1);
		assertDisulfideBond("A", "A", 14, 38, bond2);

		List<SSBondImpl> list = SSBondImpl.getSsBondListFromBondList(ssbonds);

		//System.out.println(list.get(0).toPDB());
		assertEquals("PDB representation incorrect", pdb1, list.get(0).toPDB().trim());

		//System.out.println(list.get(1).toPDB());
		assertEquals("PDB representation incorrect", pdb2, list.get(1).toPDB().trim());

	}

	private void assertDisulfideBond(String expectedChainId1, String expectedChainId2, int expectedResSerial1, int expectedResSerial2, Bond bond) {
		String chainId1 = bond.getAtomA().getGroup().getChainId();
		String chainId2 = bond.getAtomB().getGroup().getChainId();
		ResidueNumber resNum1 = bond.getAtomA().getGroup().getResidueNumber();
		ResidueNumber resNum2 = bond.getAtomB().getGroup().getResidueNumber();
		assertEquals("disulfide bond first chain id failed ", expectedChainId1, chainId1);
		assertEquals("disulfide bond second chain id failed ", expectedChainId2, chainId2);
		assertEquals("disulfide bond failed first residue number failed ", new ResidueNumber(expectedChainId1, expectedResSerial1, null), resNum1);
		assertEquals("disulfide bond failed second residue number failed ", new ResidueNumber(expectedChainId2, expectedResSerial2, null), resNum2);
	}

	/**
	 * Tests that standard amino acids are working properly
	 * @throws Exception
	 */
	@Test
	public void testStandardAmino() throws Exception {

		AminoAcid arg = StandardAminoAcid.getAminoAcid("ARG");
		assertTrue(arg.size() == 11 );

		AminoAcid gly = StandardAminoAcid.getAminoAcid("G");
		assertTrue(gly.size() == 4);

	}

	@Test
	public void testPDBHeader(){

		PDBHeader header = structure.getPDBHeader();
		String classification = header.getClassification();
		assertTrue(classification.equals("PROTEINASE INHIBITOR (TRYPSIN)"));

		String idCode = header.getIdCode();
		assertEquals("the idCode in the Header is " + idCode + " and not 5PTI, as expected","5PTI",idCode);

		float resolution = header.getResolution();
		assertEquals("the resolution in the Header is " + resolution + " and not 1.0, as expected",1.0,resolution,0.0001);

		// commenting out test for deprecated method
		//String technique = header.getTechnique();
		String techShould = "X-RAY DIFFRACTION";
		//assertEquals("the technique in the Header is " + technique, techShould,technique);

		Set<ExperimentalTechnique> techniques = header.getExperimentalTechniques();
		String technique = techniques.iterator().next().getName();
		assertEquals("the technique in the Header is " + technique, techShould, technique);


		List <EntityInfo> compounds = structure.getEntityInfos();

		// from biojava 5.0 we have limited support for old pdb files with no chain identifiers
		// due to that, we don't find all compounds in this file: 1 protein, 1 PO4, 1 UNK and 1 deuterated water entity
		// thus commenting out the test
		//assertEquals("did not find the right number of compounds! ", 2, compounds.size());

		EntityInfo comp = compounds.get(0);
		assertEquals("did not get the right compounds info",true,comp.getDescription().startsWith("TRYPSIN INHIBITOR"));

		List<String> chainIds = comp.getChainIds();
		List<Chain> chains    = comp.getChains();

		assertEquals("the number of chain ids and chains did not match!",chainIds.size(),chains.size());
		assertEquals("the chain ID did not match", chainIds.get(0),chains.get(0).getId());
	}

	@Test
	public void testCreateVirtualCBAtom(){

		Group g1 = structure.getChainByIndex(0).getAtomGroup(11);

		if ( g1.getPDBName().equals("GLY")){
			if ( g1 instanceof AminoAcid){
				try {
					Atom cb = Calc.createVirtualCBAtom((AminoAcid)g1);
					g1.addAtom(cb);
				} catch (StructureException e){
					fail ("createVirtualCBAtom failed with " + e.getMessage());
				}
			}
		} else {
			fail("the group at position 11 is not a GLY!");
		}
	}

	@Test
	public void testMutation() throws Exception {

		Group g1 = (Group)structure.getChainByIndex(0).getAtomGroup(21).clone();
		assertTrue(g1 != null);


		Group g2 = (Group)structure.getChainByIndex(0).getAtomGroup(53).clone();
		assertTrue(g2 != null);


		assertEquals("The group at position 22 is not a PHE","PHE", g1.getPDBName());
		assertEquals("The group position is  not number 22","22", g1.getResidueNumber().toString());

		assertEquals("The group at position 54 is not a THR","THR", g2.getPDBName());
		assertEquals("The group position is not number 54","54", g2.getResidueNumber().toString());

		Atom[] atoms1 = new Atom[3];
		Atom[] atoms2 = new Atom[3];

		atoms1[0] = g1.getAtom("N");
		atoms1[1] = g1.getAtom("CA");
		atoms1[2] = g1.getAtom("CB");


		atoms2[0] = g2.getAtom("N");
		atoms2[1] = g2.getAtom("CA");
		atoms2[2] = g2.getAtom("CB");

		Matrix4d transform = SuperPositions.superpose(
				Calc.atomsToPoints(atoms1), Calc.atomsToPoints(atoms2));

		Group newGroup = (Group) g2.clone();

		Calc.transform(newGroup, transform);

		Atom ca1    =       g1.getAtom("CA");
		Atom oldca2 =       g2.getAtom("CA");
		Atom newca2 = newGroup.getAtom("CA");
		Element e1 = ca1.getElement();

		assertEquals(Element.C, e1);

		// this also tests the cloning ...
		double olddistance = Calc.getDistance(ca1,oldca2);
		assertTrue( olddistance > 10 );

		// final test check that the distance between the CA atoms is small ;

		double newdistance = Calc.getDistance(ca1,newca2);
		assertTrue( newdistance < 0.1);


	}

	@Test
	public void testElement() throws Exception {
		// there should be no wild card elements
		// in a structure (!= Element.R)
		for (Chain c: structure.getChains()) {
			for (Group g: c.getAtomGroups()) {
				for (Atom a: g.getAtoms()) {
					assertFalse(a.getElement().equals(Element.R));
				}
			}
		}
	}


}
