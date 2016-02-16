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
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.biojava.nbio.structure.jama.Matrix;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.Set;

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
		pdbpars.setFileParsingParameters(params);
		
		structure = pdbpars.parsePDBFile(inStream) ;

		assertNotNull(structure);

		assertEquals("structure does not contain one chain ", 2 ,structure.size());
	}

	@Test
	public void testSeqResParsing() {

		// System.out.println(structure);
		List<Chain> chains = structure.getChains(0);
		assertEquals(" nr of found chains not correct!",2,chains.size());
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

		Chain c = structure.getChain(0);
		assertEquals("did not find the expected 58 amino acids!",58,c.getAtomGroups(GroupType.AMINOACID).size());

		assertTrue(c.getAtomGroups(GroupType.HETATM).size()     == 0);

		Chain c2 = structure.getChain(1);
		assertTrue(c2.getAtomGroups(GroupType.HETATM).size()     == 65);
		assertTrue(c2.getAtomGroups(GroupType.NUCLEOTIDE).size() == 0 );

		List<Compound> compounds= structure.getCompounds();
		
		// from Biojava 4.2 on we are creating compounds whenever an entity is found to be without an assigned compound in the file
		// see issues https://github.com/biojava/biojava/issues/305 and https://github.com/biojava/biojava/pull/394
		assertEquals(2, compounds.size());
		Compound mol = compounds.get(0);
		assertTrue(mol.getMolName().startsWith("TRYPSIN INHIBITOR"));
	}

	@Test
	public void testSSBondParsing() throws Exception {
		assertNotNull(structure);

		List<SSBond> ssbonds = structure.getSSBonds();
		assertEquals("did not find the correct nr of SSBonds ",3,ssbonds.size());

		String pdb1 = "SSBOND   1 CYS A    5    CYS A   55";
		String pdb2 = "SSBOND   2 CYS A   14    CYS A   38";

		SSBond bond1 = ssbonds.get(0);

		String b1 = bond1.toPDB();

		assertTrue("PDB representation incorrect",pdb1.equals(b1.trim()));
		assertTrue("not right resnum1 " , bond1.getResnum1().equals("5"));
		assertTrue("not right resnum2 " , bond1.getResnum2().equals("55"));

		SSBond bond2 = ssbonds.get(1);
		String b2 = bond2.toPDB();
		assertTrue("not right resnum1 " , bond2.getResnum1().equals("14"));
		assertTrue("not right resnum2 " , bond2.getResnum2().equals("38"));
		assertTrue("PDB representation incorrect",pdb2.equals(b2.trim()));

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


		List <Compound> compounds = structure.getCompounds();
		
		// from Biojava 4.2 on we are creating compounds whenever an entity is found to be without an assigned compound in the file
		// see issues https://github.com/biojava/biojava/issues/305 and https://github.com/biojava/biojava/pull/394
		assertEquals("did not find the right number of compounds! ", 2, compounds.size());

		Compound comp = compounds.get(0);
		assertEquals("did not get the right compounds info",true,comp.getMolName().startsWith("TRYPSIN INHIBITOR"));

		List<String> chainIds = comp.getChainIds();
		List<Chain> chains    = comp.getChains();

		assertEquals("the number of chain ids and chains did not match!",chainIds.size(),chains.size());
		assertEquals("the chain ID did not match", chainIds.get(0),chains.get(0).getChainID());
	}

	@Test
	public void testCreateVirtualCBAtom(){

		Group g1 = structure.getChain(0).getAtomGroup(11);

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

		Group g1 = (Group)structure.getChain(0).getAtomGroup(21).clone();
		assertTrue(g1 != null);


		Group g2 = (Group)structure.getChain(0).getAtomGroup(53).clone();
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


		SVDSuperimposer svds = new SVDSuperimposer(atoms1,atoms2);


		Matrix rotMatrix = svds.getRotation();
		Atom   tran      = svds.getTranslation();

		Group newGroup = (Group)g2.clone();

		Calc.rotate(newGroup,rotMatrix);

		Calc.shift(newGroup,tran);

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
