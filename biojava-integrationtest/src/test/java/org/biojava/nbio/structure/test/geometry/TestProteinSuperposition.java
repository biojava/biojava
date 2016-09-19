package org.biojava.nbio.structure.test.geometry;

import static org.junit.Assert.*;

import java.io.IOException;

import javax.vecmath.Point3d;

import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.geometry.SuperPosition;
import org.biojava.nbio.structure.geometry.SuperPositionQCP;
import org.biojava.nbio.structure.geometry.SuperPositionQuat;
import org.biojava.nbio.structure.geometry.SuperPositionSVD;
import org.junit.BeforeClass;
import org.junit.Test;

public class TestProteinSuperposition {
	
	private static Point3d[] chain1;
	private static Point3d[] chain2;
	
	@BeforeClass
	public static void setUpBeforeClass() throws StructureException, IOException {
		Structure s = StructureIO.getStructure("1smt");
		Chain origChainA = s.getPolyChainByPDB("A");
		Chain clonedChainA = (Chain) origChainA.clone();
		
		chain1 = Calc.atomsToPoints(StructureTools.getAtomCAArray(origChainA));
		chain2 = Calc.atomsToPoints(StructureTools.getAtomCAArray(clonedChainA));
		
	}

	@Test
	public void testSuperpositionSVD()  {
		
		SuperPosition sup = new SuperPositionSVD(false);
		
		double rmsd = sup.getRmsd(chain1, chain2);
		
		assertEquals(0.0, rmsd, 0.0001);
	}

	@Test
	public void testSuperpositionQCP() {
		
		SuperPosition sup = new SuperPositionQCP(false);
		
		double rmsd = sup.getRmsd(chain1, chain2);
		
		assertEquals(0.0, rmsd, 0.0001);
	}
	
	@Test
	public void testSuperpositionQuat() {
		
		SuperPosition sup = new SuperPositionQuat(false);
		
		double rmsd = sup.getRmsd(chain1, chain2);
		
		assertEquals(0.0, rmsd, 0.0001);
	}
}
