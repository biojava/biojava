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
package org.biojava.nbio.structure.io;

import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.geometry.CalcPoint;
import org.biojava.nbio.structure.geometry.SuperPosition;
import org.biojava.nbio.structure.geometry.SuperPositionQCP;
import org.junit.Test;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

public class TestHardBioUnits {

	/**
	 * This tests that the biounit and operator ids are the right ones when parsing from mmcif
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void test4A1Immcif() throws IOException, StructureException {

		String pdbId = "4A1I";
		int biolAssemblyNr = 2;

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		cache.setUseMmtf(false);
		StructureIO.setAtomCache(cache);
		
		Structure bioAssembly = StructureIO.getBiologicalAssembly(pdbId,biolAssemblyNr);

		if ( bioAssembly == null){
			System.err.println("Could not generate the biological assembly " + pdbId + " nr " + biolAssemblyNr);
		}


		/*
		 * loop_
				_pdbx_struct_assembly_gen.assembly_id
				_pdbx_struct_assembly_gen.oper_expression
				_pdbx_struct_assembly_gen.asym_id_list
				1 1 A,I,J,K,L,M,N,UA,H,PA,QA,RA,SA,TA,BB
				2 1 G,KA,LA,MA,NA,OA,AB
				2 2 B,O,P,Q,R,VA
				3 1 B,O,P,Q,R,VA
				3 3 G,KA,LA,MA,NA,OA,AB
				4 1 C,S,T,U,V,W,WA,F,FA,GA,HA,IA,JA,ZA
				5 1 D,X,Y,Z,XA,E,AA,BA,CA,DA,EA,YA
		 */

		//System.out.println(bioAssembly.toPDB());

		assertEquals(1, bioAssembly.nrModels());

		assertEquals(2, bioAssembly.getPolyChains().size());

		// this tests checks that the operator ids are exactly those read from mmcif, it doesn't necessarily work in mmtf where there are no ids
		Chain g = bioAssembly.getPolyChainByPDB("G_1");
		Chain b = bioAssembly.getPolyChainByPDB("B_2");
		
		assertNotNull(g);

		assertNotNull(b);

	}
	
	/**
	 * This tests that the biounit is correctly represented (should work from all sources mmcif, pdb, mmtf)
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void test4A1I() throws IOException, StructureException {

		String pdbId = "4A1I";
		int biolAssemblyNr = 2;

		Structure bioAssembly = StructureIO.getBiologicalAssembly(pdbId,biolAssemblyNr);

		if ( bioAssembly == null){
			System.err.println("Could not generate the biological assembly " + pdbId + " nr " + biolAssemblyNr);
		}


		/*
		 * loop_
				_pdbx_struct_assembly_gen.assembly_id
				_pdbx_struct_assembly_gen.oper_expression
				_pdbx_struct_assembly_gen.asym_id_list
				1 1 A,I,J,K,L,M,N,UA,H,PA,QA,RA,SA,TA,BB
				2 1 G,KA,LA,MA,NA,OA,AB
				2 2 B,O,P,Q,R,VA
				3 1 B,O,P,Q,R,VA
				3 3 G,KA,LA,MA,NA,OA,AB
				4 1 C,S,T,U,V,W,WA,F,FA,GA,HA,IA,JA,ZA
				5 1 D,X,Y,Z,XA,E,AA,BA,CA,DA,EA,YA
		 */

		//System.out.println(bioAssembly.toPDB());

		assertEquals(1, bioAssembly.nrModels());

		assertEquals(2, bioAssembly.getPolyChains().size());
		
		// here we'll store all author chain ids without the operator id part
		Set<String> chainIdsNoOps = new HashSet<String>();
		
		for (Chain poly:bioAssembly.getPolyChains()) {
			chainIdsNoOps.add(poly.getName().split("_")[0]);			
		}

		assertEquals(2, chainIdsNoOps.size());
		
		// we should have B and G only
		assertTrue(chainIdsNoOps.contains("B"));
		assertTrue(chainIdsNoOps.contains("G"));
		assertFalse(chainIdsNoOps.contains("A"));
		assertFalse(chainIdsNoOps.contains("H"));
		
		// now let's check that the right operators were applied to the right chains
		
		// first we need the original structure
		Structure original = StructureIO.getStructure(pdbId);
		
		
		Point3d[] atomsOrigChainG = Calc.atomsToPoints(StructureTools.getAtomCAArray(original.getPolyChainByPDB("G"))); 
		Point3d[] atomsOrigChainB = Calc.atomsToPoints(StructureTools.getAtomCAArray(original.getPolyChainByPDB("B")));
		
		List<Chain> bioAssemblyChains = bioAssembly.getPolyChains();
		Chain transfChainB = null;
		Chain transfChainG = null;
		// get the bioassembly's equivalent chains B and G
		for (Chain c : bioAssemblyChains) {
			if (c.getName().startsWith("B")) transfChainB = c;
			if (c.getName().startsWith("G")) transfChainG = c;
		}
		
		assertNotNull(transfChainB);
		assertNotNull(transfChainG);
		
		Point3d[] atomsTransfChainG = Calc.atomsToPoints(StructureTools.getAtomCAArray(transfChainG));
		Point3d[] atomsTransfChainB = Calc.atomsToPoints(StructureTools.getAtomCAArray(transfChainB)); 
		
		SuperPosition sqcp = new SuperPositionQCP(false);
		
		// operator 1 is the identity, trace should be == 3
		Matrix4d m1 = sqcp.superposeAndTransform(atomsOrigChainG, atomsTransfChainG);		
		assertEquals(3.0, m1.m00 + m1.m11 + m1.m22, 0.00001);	
		assertEquals(0.0, CalcPoint.rmsd(atomsOrigChainG, atomsTransfChainG), 0.00001);
		
		
		// operator 2 is a 2-fold, trace should be == -1
		Matrix4d m2 = sqcp.superposeAndTransform(atomsOrigChainB, atomsTransfChainB);
		assertEquals(-1.0, m2.m00 + m2.m11 + m2.m22, 0.00001);
		assertEquals(0.0, CalcPoint.rmsd(atomsOrigChainB, atomsTransfChainB), 0.00001);

		
	}
	
}
