package org.biojava.bio.structure;

import java.io.IOException;

import junit.framework.TestCase;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava3.structure.StructureIO;
import org.junit.Before;

public class TestBond extends TestCase {
	private Structure s;

	@Before
	public void setUp() throws IOException, StructureException {		
		AtomCache cache = new AtomCache();

		cache.setUseMmCif(false);

		FileParsingParameters params = cache.getFileParsingParams();

		params.setStoreEmptySeqRes(true);
		params.setAlignSeqRes(true);
		params.setLoadChemCompInfo(true);
		params.setCreateAtomBonds(true);

		StructureIO.setAtomCache(cache);


	}

	public void testIntraResidueBonds() throws StructureException, IOException {

		if ( s == null) {
			setUp();
			s = StructureIO.getStructure("1kh9");
		}

	


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

	public void testPeptideBonds() throws StructureException, IOException {
		if ( s == null) {
			setUp();
			s = StructureIO.getStructure("1kh9");
		}

		AminoAcidImpl residue1 = (AminoAcidImpl) s.getChainByPDB("A").getSeqResGroup(273);
		AminoAcidImpl residue2 = (AminoAcidImpl) s.getChainByPDB("A").getSeqResGroup(274);

		Atom carboxylC = residue1.getC();
		Atom aminoN = residue2.getN();

		assertTrue(areBonded(carboxylC, aminoN));
	}

	public void testLINKBonds() throws StructureException, IOException {
		if ( s == null) {
			setUp();
			s = StructureIO.getStructure("1kh9");
		}
		
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

	public void testDisulfideBonds() throws StructureException, IOException {
		if ( s == null) {
			setUp();
			s = StructureIO.getStructure("1kh9");
		}
		Atom atom1 = s.getChainByPDB("A").getSeqResGroup(177).getAtom("SG");
		Atom atom2 = s.getChainByPDB("A").getSeqResGroup(167).getAtom("SG");

		assertTrue(areBonded(atom1, atom2));
	}

	public void testLigandBonds() throws StructureException, IOException {
		if ( s == null) {
			setUp();
			s = StructureIO.getStructure("1kh9");
		}
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

	public void test145D() throws IOException, StructureException {
		StructureIO.getStructure("145D");
	}

	public void test1APJ() throws IOException, StructureException {
		StructureIO.getStructure("1APJ");
	}

	public void test1BDX() throws IOException, StructureException {
		StructureIO.getStructure("1BDX");
	}
}
