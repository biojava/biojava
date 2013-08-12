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
		FileParsingParameters params = cache.getFileParsingParams();
		
		params.setStoreEmptySeqRes(true);
		params.setAlignSeqRes(true);
		params.setLoadChemCompInfo(true);
		params.setCreateAtomBonds(true);
		
		StructureIO.setAtomCache(cache);
		s = StructureIO.getStructure("1KH9");
	}
	
	public void testIntraResidueBonds() throws StructureException {
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
	
	public void testPeptideBonds() throws StructureException {
		AminoAcidImpl residue1 = (AminoAcidImpl) s.getChainByPDB("A").getSeqResGroup(273);
		AminoAcidImpl residue2 = (AminoAcidImpl) s.getChainByPDB("A").getSeqResGroup(274);
		
		Atom carboxylC = residue1.getC();
		Atom aminoN = residue2.getN();
		
		assertTrue(areBonded(carboxylC, aminoN));
	}
	
	public void testLINKBonds() throws StructureException {
		Atom atom1 = s.getChainByPDB("A").getSeqResGroup(50).getAtom("OD1");
		Atom atom2 = s.getChainByPDB("A").getAtomGroup(446).getAtom("MG");
		
		assertTrue(areBonded(atom1, atom2));
	}
	
	public void testDisulfideBonds() throws StructureException {
		Atom atom1 = s.getChainByPDB("A").getSeqResGroup(177).getAtom("SG");
		Atom atom2 = s.getChainByPDB("A").getSeqResGroup(167).getAtom("SG");
		
		assertTrue(areBonded(atom1, atom2));
	}
	
	public void testLigandBonds() throws StructureException {
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
}
