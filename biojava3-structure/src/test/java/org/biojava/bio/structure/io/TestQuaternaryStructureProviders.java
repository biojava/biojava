package org.biojava.bio.structure.io;

import static org.junit.Assert.*;

import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.PDBHeader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.quaternary.ModelTransformationMatrix;
import org.biojava.bio.structure.quaternary.io.BioUnitDataProviderFactory;
import org.biojava.bio.structure.quaternary.io.MmCifBiolAssemblyProvider;
import org.biojava.bio.structure.quaternary.io.PDBBioUnitDataProvider;
import org.biojava3.structure.StructureIO;
import org.junit.Test;

public class TestQuaternaryStructureProviders {

	@Test
	public void test1STP() {
		
		
		testID("1stp",1);
	}
	@Test
	public void test3FAD(){
		testID("3fad",1);
		testID("3fad",2);
	}
	
	@Test
	public void test5LDH(){
		testID("5LDH",1);

	}
	
	@Test
	public void test2TBV(){
		testID("2TBV",1);

	}
//	@Test
//	public void test1EI7(){
//		testID("1ei7",1);
//
//	}
	
		
	private void testID(String pdbId, int bioMolecule){
		
		
		try {
			// get bio assembly from PDB file
			PDBBioUnitDataProvider pdbProvider = new PDBBioUnitDataProvider();
			BioUnitDataProviderFactory.setBioUnitDataProvider(pdbProvider.getClass().getCanonicalName());
			Structure pdbS = StructureIO.getBiologicalAssembly(pdbId, bioMolecule);

			// get bio assembly from mmcif file
			MmCifBiolAssemblyProvider mmcifProvider = new MmCifBiolAssemblyProvider();
			BioUnitDataProviderFactory.setBioUnitDataProvider(mmcifProvider.getClass().getCanonicalName());			
			Structure mmcifS = StructureIO.getBiologicalAssembly(pdbId, bioMolecule);
		
			BioUnitDataProviderFactory.setBioUnitDataProvider(BioUnitDataProviderFactory.DEFAULT_PROVIDER_CLASSNAME);
			
			
			
			PDBHeader pHeader = pdbS.getPDBHeader();
			PDBHeader mHeader = mmcifS.getPDBHeader();
			//PDBHeader fHeader = flatFileS.getPDBHeader();
			
			assertTrue("not correct nr of bioassemblies " + pHeader.getNrBioAssemblies() + " " , pHeader.getNrBioAssemblies() >= bioMolecule);
			assertTrue("not correct nr of bioassemblies " + mHeader.getNrBioAssemblies() + " " , mHeader.getNrBioAssemblies() >= bioMolecule);
			//assertTrue("not correct nr of bioassemblies " + fHeader.getNrBioAssemblies() + " " , fHeader.getNrBioAssemblies() >= bioMolecule);
			
			// mmcif files contain sometimes partial virus assemblies, so they can contain more info than pdb
			assertTrue(pHeader.getNrBioAssemblies() <= mHeader.getNrBioAssemblies());
			
			
			Map<Integer, List<ModelTransformationMatrix>> pMap = pHeader.getBioUnitTranformationMap();
			Map<Integer, List<ModelTransformationMatrix>> mMap = mHeader.getBioUnitTranformationMap();
			
			//System.out.println("PDB: " + pMap);
			
			//System.out.println("Mmcif: " + mMap);
			
			assertTrue(pMap.keySet().size()<= mMap.keySet().size());
			
			for ( Integer k : pMap.keySet()) {
				assertTrue(mMap.containsKey(k));
				
				List<ModelTransformationMatrix> pL = pMap.get(k);
				
				// mmcif list can be longer due to the use of internal chain IDs
				List<ModelTransformationMatrix> mL = mMap.get(k);
				
				//assertEquals(pL.size(), mL.size());
				
				
				for (ModelTransformationMatrix m1 : pL){
					
					boolean found = false;
					for ( ModelTransformationMatrix m2 : mL){
						
						if  (! m1.getNdbChainId().equals(m2.getNdbChainId()))
								continue;
						if ( ! m1.getMatrix().toString().equals(m2.getMatrix().toString()))
								continue;
						if ( ! equalVectors(m1.getVector(),m2.getVector()))
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
			
			
			// compare with flat file version:
			AtomCache cache = new AtomCache();
			FileParsingParameters params = cache.getFileParsingParams();
			params.setAlignSeqRes(true);
			params.setParseCAOnly(false);

			Structure flatFileS = cache.getBiologicalAssembly(pdbId, bioMolecule, false);
			
			Atom[] fileA = StructureTools.getAllAtomArray(flatFileS);
			
			assertEquals(pdbA.length, fileA.length);
			
			assertEquals(pdbS.nrModels(),flatFileS.nrModels());
						
		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}
		
		
	
		
		
	}
	private boolean equalVectors(double[] vector, double[] vector2) {
		
		String s1 = String.format("%.5f %.5f %.5f", vector[0], vector[1], vector[2]);
		String s2 = String.format("%.5f %.5f %.5f", vector2[0], vector2[1], vector2[2]);
		//System.out.println(s1 + " " + s2);
		return s1.equals(s2);
		
	}

}
