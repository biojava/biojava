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

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.junit.Test;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class TestQuaternaryStructureProviders {

	@Test
	public void test1STP() throws IOException, StructureException{
		comparePdbVsMmcif("1stp",1, 4);
	}

	@Test
	public void test3FAD() throws IOException, StructureException{
		comparePdbVsMmcif("3fad",1, 1);
		comparePdbVsMmcif("3fad",2, 2);
	}

	@Test
	public void test5LDH() throws IOException, StructureException{
		comparePdbVsMmcif("5LDH",1, 4);
		
		// the pdb file of 5ldh contains only 1 bioassembly, whilst the mmcif contains 2,
		// thus we can't test here the comparison between the 2
		//testID("5LDH",2, 2); 
		
		// since v5 remediation there's 4 bioassemblies with numerical ids for 5ldh, no more PAU and XAU
		boolean gotException = false;
		try {
			AtomCache cache = new AtomCache();
			cache.setUseMmCif(true);
			StructureIO.setAtomCache(cache);
			StructureIO.getBiologicalAssembly("5LDH",3);
		} catch (StructureException e) {
			gotException = true;
		}

		assertTrue("Bioassembly 3 for PDB id 5LDH should fail with a StructureException!", !gotException);

		// bioassembly 2 does exist in mmcif file, let's check that
		gotException = false;
		try {
			AtomCache cache = new AtomCache();
			cache.setUseMmCif(true);
			StructureIO.setAtomCache(cache);
			StructureIO.getBiologicalAssembly("5LDH",2);
		} catch (StructureException e) {
			gotException = true;
		}
		assertTrue("Bioassembly 2 for PDB id 5LDH should not fail with a StructureException!", !gotException);

	}

	@Test
	public void test3NTU() throws IOException, StructureException{
		comparePdbVsMmcif("3NTU",1, 6);
	}

	@Test
	public void test1A29() throws IOException, StructureException{
		comparePdbVsMmcif("1A29",1, 1);
	}

	@Test
	public void test1EI7() throws IOException, StructureException {

		comparePdbVsMmcif("1ei7",1, 68);
		
	}

	@Test
	public void testGetNrBioAssemblies5LDH() throws IOException, StructureException {
		assertEquals("There should be 4 bioassemblies for 5LDH, see github issue #230", 4, StructureIO.getBiologicalAssemblies("5LDH").size());
	}


	/**
	 * Bioassembly tests for a single PDB entry 
	 * @param pdbId
	 * @param bioMolecule the bio assembly identifier to test
	 * @param mmSize the expected mmSize of given bioMolecule number
	 * @throws IOException
	 * @throws StructureException
	 */
	private void comparePdbVsMmcif(String pdbId, int bioMolecule, int mmSize) throws IOException, StructureException{

			
		Structure pdbS = getPdbBioAssembly(pdbId, bioMolecule, true);

		Structure mmcifS = getMmcifBioAssembly(pdbId, bioMolecule, true);

		PDBHeader pHeader = pdbS.getPDBHeader();
		PDBHeader mHeader = mmcifS.getPDBHeader();

		assertTrue("not correct nr of bioassemblies " + pHeader.getNrBioAssemblies() + " " , pHeader.getNrBioAssemblies() >= bioMolecule);
		assertTrue("not correct nr of bioassemblies " + mHeader.getNrBioAssemblies() + " " , mHeader.getNrBioAssemblies() >= bioMolecule);

		// mmcif files contain sometimes partial virus assemblies, so they can contain more info than pdb
		assertTrue(pHeader.getNrBioAssemblies() <= mHeader.getNrBioAssemblies());


		Map<Integer, BioAssemblyInfo> pMap = pHeader.getBioAssemblies();
		Map<Integer, BioAssemblyInfo> mMap = mHeader.getBioAssemblies();


		assertTrue(pMap.keySet().size()<= mMap.keySet().size());
		
		assertEquals(mmSize, mMap.get(bioMolecule).getMacromolecularSize());


		for ( int k : pMap.keySet()) {
			assertTrue(mMap.containsKey(k));
			
			BioAssemblyInfo pBioAssemb = pMap.get(k);
			BioAssemblyInfo mBioAssemb = mMap.get(k);

			assertEquals("Macromolecular sizes don't coincide!",pBioAssemb.getMacromolecularSize(), mBioAssemb.getMacromolecularSize());
			
			List<BiologicalAssemblyTransformation> pL = pBioAssemb.getTransforms();

			// mmcif list can be longer due to the use of internal chain IDs
			List<BiologicalAssemblyTransformation> mL = mBioAssemb.getTransforms();

			//assertEquals(pL.size(), mL.size());


			for (BiologicalAssemblyTransformation m1 : pL){

				boolean found = false;
				for ( BiologicalAssemblyTransformation m2 : mL){

					if  (! m1.getChainId().equals(m2.getChainId()))
						continue;

					if ( ! m1.getTransformationMatrix().epsilonEquals(m2.getTransformationMatrix(), 0.0001))
						continue;

					found = true;

				}

				if ( ! found ){
					System.err.println("did not find matching matrix " + m1);
					System.err.println(mL);
				}
				assertTrue(found);

			}
		}


		assertEquals("Not the same number of chains!" , pdbS.size(),mmcifS.size());

		Atom[] pdbA = StructureTools.getAllAtomArray(pdbS);

		Atom[] mmcifA = StructureTools.getAllAtomArray(mmcifS);

		assertEquals(pdbA.length, mmcifA.length);

		assertEquals(pdbA[0].toPDB(), mmcifA[0].toPDB());


		

	}

	private Structure getPdbBioAssembly(String pdbId, int bioMolecule, boolean multiModel) throws IOException, StructureException {
		// get bio assembly from PDB file
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(false); 
		StructureIO.setAtomCache(cache);		
		Structure pdbS = StructureIO.getBiologicalAssembly(pdbId, bioMolecule, multiModel);
		return pdbS;
	}
	
	private Structure getMmcifBioAssembly(String pdbId, int bioMolecule, boolean multiModel) throws IOException, StructureException {
		// get bio assembly from mmcif file
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true); 
		StructureIO.setAtomCache(cache);		
		Structure mmcifS = StructureIO.getBiologicalAssembly(pdbId, bioMolecule, multiModel);
		return mmcifS;
	}
	


}
