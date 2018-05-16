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
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.junit.BeforeClass;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static org.junit.Assert.*;


public class TestBond {
	
	private static final Logger logger = LoggerFactory.getLogger(TestBond.class);


	private static AtomCache cache;

	@BeforeClass
	public static void setUp() {
		
		// important: without this the tests can fail when running in maven (but not in IDE)
		// that's because it depends on the order on how tests were run - JD 2018-03-10
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider()); 
		
		cache = new AtomCache();

		cache.setUseMmCif(true);

		FileParsingParameters params = cache.getFileParsingParams();

		params.setAlignSeqRes(true);
		params.setCreateAtomBonds(true);

		StructureIO.setAtomCache(cache);


	}

	@Test
	public void testStructConnModels() throws IOException, StructureException {
		Structure s = StructureIO.getStructure("1cdr");
		Group groupOne = s.getPolyChain("A",1).getGroupByPDB(new ResidueNumber("A", 18, ' '));
		Group groupTwo = s.getNonPolyChain("B",1).getGroupByPDB(new ResidueNumber("A", 78, ' '));
		Atom atomOne = groupOne.getAtom("ND2");
		Atom atomTwo = groupTwo.getAtom("C1");
		assertTrue(areBonded(atomOne, atomTwo));
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
		assertEquals(0, countAtomsWithoutBonds("3jtm"));
		assertEquals(0, countAtomsWithoutBonds("3jq8"));
		assertEquals(0, countAtomsWithoutBonds("3jq9"));
		assertEquals(0, countAtomsWithoutBonds("3i06"));
		assertEquals(0, countAtomsWithoutBonds("3nu3"));
		assertEquals(0, countAtomsWithoutBonds("3nu4"));
		assertEquals(0, countAtomsWithoutBonds("3nvd"));
	}


	/**
	 * Test whether these partial occupancy hydrogens are bonded to the residue.
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void testHydrogenToProteinBondMissing() throws IOException, StructureException {
//		assertEquals(0, countAtomsWithoutBonds("4txr"));
//      This test has been commented out because of a data error in 4txr.cif:
//	    On 2018-03-29 the hydrogen names for EDO were modified, however, the PDB structures were not updated 
//	    (see snippets below). Note, the same issue may affect other ligands where atom names were updated.
//      PWR filed a ticket with RCSB PDB on 2018-05-16.
//	    TODO: uncomment test after data issue has been fixed.
	    
//	    edo.cif
//	    _chem_comp.pdbx_modified_date 2018-03-29
//	    ...
//	    EDO H1 H1 H 0 1 N N N 40.269 18.578 -36.179 -0.626 0.552 -1.391 H1 EDO 5
//	    EDO H2 H2 H 0 1 N N N 38.845 17.541 -35.829 -1.224 1.508 -0.013 H2 EDO 6
//	    EDO H3 H3 H 0 1 N N N 40.451 16.981 -38.299 1.224 1.508 0.014 H3 EDO 7
//	    EDO H4 H4 H 0 1 N N N 38.587 16.690 -39.453 2.328 -0.604 0.176 H4 EDO 8
//	    EDO H5 H5 H 0 1 N N N 40.908 16.123 -36.247 -2.328 -0.604 -0.177 H5 EDO 9
//	    EDO H6 H6 H 0 1 N N N 39.425 18.451 -38.410 0.626 0.551 1.391 H6 EDO 10
//
//	    4txr.cif
//	    HETATM 4078 H H11 . EDO E 5 . ? -2.936 -17.609 18.748 1.00 14.13 ? 202 EDO A H11 1
//	    HETATM 4079 H H12 . EDO E 5 . ? -1.441 -18.486 18.410 1.00 14.13 ? 202 EDO A H12 1
//	    HETATM 4080 H HO1 . EDO E 5 . ? -3.099 -19.869 19.278 1.00 16.64 ? 202 EDO A HO1 1
//	    HETATM 4081 H H21 . EDO E 5 . ? -2.090 -18.733 16.039 1.00 15.08 ? 202 EDO A H21 1
//	    HETATM 4082 H H22 . EDO E 5 . ? -2.056 -17.033 16.508 1.00 15.08 ? 202 EDO A H22 1
//	    HETATM 4083 H HO2 . EDO E 5 . ? -3.968 -17.575 15.297 1.00 15.81 ? 202 EDO A HO2 1
		assertEquals(0, countAtomsWithoutBonds("3nvd"));
	}

	/**
	 * Test whether these partial occupancy hydrogens are bonded to the residue.
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void testAltLocBondMissing() throws IOException, StructureException {
		assertEquals(0, countAtomsWithoutBonds("4cup"));
	}

	/**
	 * Loops through whole structure counting all atoms (in groups larger than 1 atom)
	 * that have no bonds.
	 * @throws IOException
	 * @throws StructureException
	 */
	private int countAtomsWithoutBonds(String pdbId) throws IOException, StructureException { 
		Structure inputStructure = StructureIO.getStructure(pdbId);
		// Loop through the structure
		int nonBondedCounter = 0;
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
							logger.debug("Atom {}-{} has no bonds", a.getPDBserial(), a.getName());
							nonBondedCounter++;
						}
					}

				}
			}
		}
		return nonBondedCounter;
	}

	private int countBondedToSelf(String pdbId) throws IOException, StructureException {
		Structure inputStructure = StructureIO.getStructure(pdbId);
		int bondedToSelf =0;
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
						if(a.getBonds()!=null){
							for(Bond b: a.getBonds()){
								if(b.getAtomA().equals(b.getAtomB())){
									bondedToSelf+=1;
								}
							}
						}
					}
				}
			}
		}
		return bondedToSelf;
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
	 * @throws IOException an error getting the required file
	 * @throws StructureException an error parsing the required file
	 */
	@Test
	public void testDeuterated() throws IOException, StructureException {
		// The terminal Hydrogen D3 - is missing (from the CCD)
		assertEquals(2, countAtomsWithoutBonds("1GKT"));
		assertEquals(2, countAtomsWithoutBonds("1IO5"));
		// All H/D2,H/D3 errors
		assertEquals(13, countAtomsWithoutBonds("5E5J"));
	}

	/**
	 * Test this weird case - with missing Oxygen atoms, alternate locations on Deuterium 
	 * and terminal hydrogens.
	 * @throws IOException an error getting the required file
	 * @throws StructureException an error parsing the required file
	 */
	@Test
	public void testWeirdCase() throws IOException, StructureException {
		assertEquals(6, countAtomsWithoutBonds("1IU6"));
	}


	/**
	 * Test that Sulphur atoms are not found to be bonded to themselves
	 * @throws IOException an error getting the required file
	 * @throws StructureException an error parsing the required file
	 */
	@Test
	public void testSSBonds() throws IOException, StructureException {
		for(String pdbCode : new String[]{"3ZXW","1NTY", "4H2I", "2K6D", "2MLM"}){
			assertEquals(0, countBondedToSelf(pdbCode));
		}
	}

}
