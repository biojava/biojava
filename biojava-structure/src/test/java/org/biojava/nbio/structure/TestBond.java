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

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.*;


public class TestBond {
	

	@BeforeClass
	public static void setUp() throws IOException, StructureException {		
		AtomCache cache = new AtomCache();

		cache.setUseMmCif(false);

		FileParsingParameters params = cache.getFileParsingParams();
		
		params.setAlignSeqRes(true);
		params.setLoadChemCompInfo(true);
		params.setCreateAtomBonds(true);

		StructureIO.setAtomCache(cache);


	}

	@Test
	public void testIntraResidueBonds() throws StructureException, IOException {

		
		Structure s = StructureIO.getStructure("1kh9");
		

		Group g = s.getChainByPDB("A").getSeqResGroup(274);
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
		

		AminoAcidImpl residue1 = (AminoAcidImpl) s.getChainByPDB("A").getSeqResGroup(273);
		AminoAcidImpl residue2 = (AminoAcidImpl) s.getChainByPDB("A").getSeqResGroup(274);

		Atom carboxylC = residue1.getC();
		Atom aminoN = residue2.getN();

		assertTrue(areBonded(carboxylC, aminoN));
	}

	@Test
	public void testLINKBonds() throws StructureException, IOException {

		Structure s = StructureIO.getStructure("1kh9");
		
		Group g1 = s.getChainByPDB("A").getSeqResGroup(50);
		assertNotNull(g1);
		
		assertTrue(g1 instanceof AminoAcid);

		AminoAcid aa = (AminoAcid)g1;
		assertTrue(aa.getRecordType().equals(AminoAcid.ATOMRECORD));

		Atom atom1 = g1.getAtom("OD1");
		Atom atom2 = s.getChainByPDB("A").getAtomGroup(446).getAtom("MG");
		assertNotNull(atom1);
		assertNotNull(atom2);
		assertTrue(areBonded(atom1, atom2));
	}

	@Test
	public void testDisulfideBonds() throws StructureException, IOException {
		
		Structure s = StructureIO.getStructure("1kh9");
		
		Atom atom1 = s.getChainByPDB("A").getSeqResGroup(177).getAtom("SG");
		Atom atom2 = s.getChainByPDB("A").getSeqResGroup(167).getAtom("SG");

		assertTrue(areBonded(atom1, atom2));
	}

	@Test
	public void testLigandBonds() throws StructureException, IOException {
		
		Structure s = StructureIO.getStructure("1kh9");
		
		Atom phosphateP = s.getChainByPDB("A").getAtomGroup(447).getAtom("P");
		Atom phosphateO = s.getChainByPDB("A").getAtomGroup(447).getAtom("O1");

		assertTrue(areBonded(phosphateP, phosphateO));
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
}
