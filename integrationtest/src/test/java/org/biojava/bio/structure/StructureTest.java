package org.biojava.bio.structure;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.Map;

import junit.framework.TestCase;

import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileParser;
import org.biojava.bio.structure.jama.Matrix;


/**
 *
 * @author Andreas Prlic
 * @since 1.5
 */


public class StructureTest extends TestCase {

	Structure structure;

	protected void setUp()
	{
		InputStream inStream = this.getClass().getResourceAsStream("/5pti_old.pdb");
		assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		pdbpars.setFileParsingParameters(params);
		
		try {
			structure = pdbpars.parsePDBFile(inStream) ;
		} catch (IOException e) {
			e.printStackTrace();
		}

		assertNotNull(structure);

		assertEquals("structure does not contain one chain ", 2 ,structure.size());
	}

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


	/** test if a PDB file can be parsed
	 * @throws Exception */
	public void testReadPDBFile() throws Exception {

		assertEquals("pdb code not set!","5PTI",structure.getPDBCode());

		Chain c = structure.getChain(0);
		assertEquals("did not find the expected 58 amino acids!",58,c.getAtomGroups("amino").size());

		assertTrue(c.getAtomGroups("hetatm").size()     == 0);

		Chain c2 = structure.getChain(1);
		assertTrue(c2.getAtomGroups("hetatm").size()     == 65);
		assertTrue(c2.getAtomGroups("nucleotide").size() == 0 );

		List<Compound> compounds= structure.getCompounds();
		assertTrue(compounds.size() == 1);
		Compound mol = compounds.get(0);
		assertTrue(mol.getMolName().startsWith("TRYPSIN INHIBITOR"));
	}


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

	/** Tests that standard amino acids are working properly
	 * @throws Exception */
	public void testStandardAmino() throws Exception {

		AminoAcid arg = StandardAminoAcid.getAminoAcid("ARG");
		assertTrue(arg.size() == 11 );

		AminoAcid gly = StandardAminoAcid.getAminoAcid("G");
		assertTrue(gly.size() == 4);

	}


	@SuppressWarnings("deprecation")
	public void testHeader() {
		Map<String, Object> m = structure.getHeader();

		assertNotNull(m);

		String classification = (String)m.get("classification");
		assertTrue(classification.equals("PROTEINASE INHIBITOR (TRYPSIN)"));

		String idCode = (String)m.get("idCode");
		assertEquals("the idCode in the Header is " + idCode + " and not 5PTI, as expected","5PTI",idCode);

		Float resolution = (Float) m.get("resolution");
		assertEquals("the resolution in the Header is " + resolution + " and not 1.0, as expected",new Float(1.0),resolution);

		String technique = (String) m.get("technique");
		String techShould = "X-RAY DIFFRACTION ";
		assertEquals("the technique in the Header is " + technique, techShould,technique);

		List <Compound> compounds = structure.getCompounds();
		assertEquals("did not find the right number of compounds! ", 1, compounds.size());

		Compound comp = compounds.get(0);
		assertEquals("did not get the right compounds info",true,comp.getMolName().startsWith("TRYPSIN INHIBITOR"));

		List<String> chainIds = comp.getChainId();
		List<Chain> chains    = comp.getChains();

		assertEquals("the number of chain ids and chains did not match!",chainIds.size(),chains.size());
		assertEquals("the chain ID did not match", chainIds.get(0),chains.get(0).getName());


	}


	@SuppressWarnings("deprecation")
	public void testPDBHeader(){
		Map<String, Object> m = structure.getHeader();
		PDBHeader header = structure.getPDBHeader();
		String classification = (String)m.get("classification");
		assertTrue(classification.equals(header.getClassification()));

		String idCode = (String)m.get("idCode");
		assertTrue(idCode.equals(header.getIdCode()));

		Float resolution = (Float) m.get("resolution");
		assertTrue(resolution.floatValue() == header.getResolution());

		String technique = (String) m.get("technique");
		assertTrue(technique.equals(header.getTechnique()));

	}

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
