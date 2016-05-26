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

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.*;


public class TestBond {


	private static AtomCache cache;

	@BeforeClass
	public static void setUp() {
		cache = new AtomCache();

		cache.setUseMmCif(true);

		FileParsingParameters params = cache.getFileParsingParams();

		params.setAlignSeqRes(true);
		params.setCreateAtomBonds(true);

		StructureIO.setAtomCache(cache);


	}

	@Test
	public void testIntraResidueBonds() throws StructureException, IOException {


		Structure s = StructureIO.getStructure("1kh9");


		Group g = s.getPolyChainByPDB("A").getSeqResGroup(274);
		Atom cg = g.getAtom("CG");

		Atom cb = g.getAtom("CB");
		Atom cd1 = g.getAtom("CD1");
		Atom cd2 = g.getAtom("CD2");

		assertEquals(3, cg.getBonds().size());
		for (Bond bond : cg.getBonds()) {
			if (bond.getOther(cg) == cb) {
				assertEquals(1, bond.getBondOrder());
			} else if (bond.getOther(cg) == cd1) {
				assertEquals(2, bond.getBondOrder());
			} else if (bond.getOther(cg) == cd2) {
				assertEquals(1, bond.getBondOrder());
			}
		}
	}

	@Test
	public void testPeptideBonds() throws StructureException, IOException {

		Structure s = StructureIO.getStructure("1kh9");


		AminoAcidImpl residue1 = (AminoAcidImpl) s.getPolyChainByPDB("A").getSeqResGroup(273);
		AminoAcidImpl residue2 = (AminoAcidImpl) s.getPolyChainByPDB("A").getSeqResGroup(274);

		Atom carboxylC = residue1.getC();
		Atom aminoN = residue2.getN();

		assertTrue(areBonded(carboxylC, aminoN));
	}

	@Test
	public void testDisulfideBonds() throws StructureException, IOException {

		Structure s = StructureIO.getStructure("1kh9");

		Atom atom1 = s.getPolyChainByPDB("A").getSeqResGroup(177).getAtom("SG");
		Atom atom2 = s.getPolyChainByPDB("A").getSeqResGroup(167).getAtom("SG");

		assertTrue(areBonded(atom1, atom2));
	}

	@Test
	public void testLigandBonds() throws StructureException, IOException {

		Structure s = StructureIO.getStructure("1kh9");

		Atom phosphateP = s.getChain("I").getAtomGroup(0).getAtom("P");
		Atom phosphateO = s.getChain("I").getAtomGroup(0).getAtom("O1");

		assertTrue(areBonded(phosphateP, phosphateO));
	}

	/**
	 * Test whether nucleotide bonds are being generated
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test 
	public void testNucleotideBonds() throws IOException, StructureException {
		Structure bio = StructureIO.getStructure("4y60");
		for( Chain c : bio.getChains()) {
			int groupCounter = 0;
			List<Group> currentGroups = c.getAtomGroups();
			for ( Group g : currentGroups) {
				if(groupCounter!=0 && groupCounter<currentGroups.size()) {
					List<Atom> atoms = g.getAtoms();
					for ( Atom a : atoms) {
						if ( a.getName().equals("P")){
							// Check to see if one of the phosphate atoms has bonding to something
							// outside of the group.
							List<Integer> indexList = new ArrayList<>();
							for (Bond b : a.getBonds()){
								indexList.add(atoms.indexOf(b.getOther(a)));
							}
							assertTrue(indexList.contains(-1));
						}
					}
				}
				groupCounter++;
			}

		}
	}

	/**
	 * Test whether these partial occupancy hydrogens are bonded to the residue.
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void testHeavyAtomBondMissing() throws IOException, StructureException {
		assertEquals(testMissingBonds("3jtm"),0);
		assertEquals(testMissingBonds("3jq8"),0);
		assertEquals(testMissingBonds("3jq9"),0);
		assertEquals(testMissingBonds("3i06"),0);
		assertEquals(testMissingBonds("3nu3"),0);
		assertEquals(testMissingBonds("3nu4"),0);
		assertEquals(testMissingBonds("3nvd"),0);
	}


	/**
	 * Test whether these partial occupancy hydrogens are bonded to the residue.
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void testHydrogenToProteinBondMissing() throws IOException, StructureException {
		assertEquals(testMissingBonds("4txr"),0);
		assertEquals(testMissingBonds("3nvd"),0);
	}
	
	/**
	 * Test whether these partial occupancy hydrogens are bonded to the residue.
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void testAltLocBondMissing() throws IOException, StructureException {
		assertEquals(testMissingBonds("4cup"),0);
	}

	/**
	 * 
	 * @throws IOException
	 * @throws StructureException
	 */
	private int testMissingBonds(String pdbId) throws IOException, StructureException { 
		Structure inputStructure = StructureIO.getStructure(pdbId);
		// Loop through the structure
		int nonBondedCounter =0;
		for(int i=0;i<inputStructure.nrModels();i++){
			for(Chain c: inputStructure.getChains(i)){
				for(Group g: c.getAtomGroups()){
					// Skip single atom groups
					if(g.size()<=1){
						continue;
					}
					// Get all the atoms
					List<Atom> atomsList = new ArrayList<>(g.getAtoms());
					for(Group altLocOne: g.getAltLocs()){
						atomsList.addAll(altLocOne.getAtoms());
					}
					// Check they all have bonds
					for(Atom a: atomsList){
						if(a.getBonds()==null){
							nonBondedCounter++;
						}
					}

				}
			}
		}
		return nonBondedCounter;

	}


	private boolean areBonded(Atom a, Atom b) {
		for (Bond bond : a.getBonds()) {
			if (bond.getOther(a) == b) {
				return true;
			}
		}

		return false;
	}

	/*
	 * Each of the following PDB IDs used to make formBonds() crash.
	 */

	@Test
	public void test145D() throws IOException, StructureException {
		StructureIO.getStructure("145D");
	}

	@Test
	public void test1APJ() throws IOException, StructureException {
		StructureIO.getStructure("1APJ");
	}

	@Test
	public void test1BDX() throws IOException, StructureException {
		StructureIO.getStructure("1BDX");
	}
	
	/**
	 * Test that all the atoms in deuterated structures are bonded.
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void testDeuterated() throws IOException, StructureException {
		// The terminal Hydrogen D3 - is missing (from the CCD)
		assertEquals(testMissingBonds("1GKT"),1);
		assertEquals(testMissingBonds("1IO5"),2);
		// All H/D2,H/D3 errors
		assertEquals(testMissingBonds("5E5J"),13);
	}
	
	/**
	 * Test this weird case - with missing Oxygen atoms, alternate locations on Deuterium 
	 * and terminal hydrogens.
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void testWeirdCase() throws IOException, StructureException {
		assertEquals(testMissingBonds("1IU6"),6);
	}

}
