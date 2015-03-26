package org.biojava.nbio.structure;

import static org.junit.Assert.*;

import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.junit.Test;

public class TestCompoundResIndexMapping {

	private static final String PATH_TO_TEST_FILES = "/org/biojava/nbio/structure/io/";
	
	
	@Test
	public void test1B8G() throws IOException, StructureException { 

		
		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure s = StructureIO.getStructure("1B8G");
		
		Chain chainA = s.getChainByPDB("A");
		int i = chainA.getCompound().getAlignedResIndex(chainA.getAtomGroup(0),chainA);
		assertEquals("First residue in 1b8gA "+chainA.getAtomGroup(0).toString()+" should map to 1 in SEQRES",1,i);
		
		checkAllResidues(s);		
		
	}

	@Test
	public void test1SMT() throws IOException, StructureException { 

		
		AtomCache cache = new AtomCache();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);

		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure s = StructureIO.getStructure("1SMT");
		
		Chain chainA = s.getChainByPDB("A");
		int i = chainA.getCompound().getAlignedResIndex(chainA.getAtomGroup(0),chainA);
		assertEquals("First residue in 1smtA "+chainA.getAtomGroup(0).toString()+" should map to 24 in SEQRES",24,i);
		Chain chainB = s.getChainByPDB("B");
		i = chainB.getCompound().getAlignedResIndex(chainB.getAtomGroup(0),chainB);
		assertEquals("First residue in 1smtB "+chainA.getAtomGroup(0).toString()+" should map to 20 in SEQRES",20,i);
		
		checkAllResidues(s);
	}
	
	private void checkAllResidues(Structure s) {
		for (Chain c:s.getChains()) {
			for (Group g:c.getAtomGroups()) {
				if (g.getType() != GroupType.AMINOACID) continue;
				
				int resIndex = c.getCompound().getAlignedResIndex(g, c);
				assertNotEquals("Residue "+g.getResidueNumber()+" should not map to -1 in SEQRES",-1, resIndex);
				assertNotEquals("Residue "+g.getResidueNumber()+" should not map to 0 in SEQRES",0, resIndex);
				assertTrue(resIndex <= c.getSeqResLength());
			}
		}
	}
	
	//@Test
	public void test3ddoRawNoSeqres() throws IOException, StructureException { 
		
		// 3ddo has 6 chains in 1 entity, all of them with different residue numbering (chain A is 1000+, chain B 2000+ ...)
		Structure s = getStructure("3ddo_raw_noseqres.pdb.gz", true);
		
		assertEquals(1,s.getCompounds().size());
		
		Chain chainA = s.getChainByPDB("A");
		Chain chainB = s.getChainByPDB("B");
		Chain chainC = s.getChainByPDB("C");

		// the last 3 residues of each chain are all the same: they should map to the same index
		for (int resNum=251;resNum<=253;resNum++) {
			Group groupInChainA = chainA.getGroupByPDB(new ResidueNumber("A", resNum+1000, null));
			int indexInChainA = chainA.getCompound().getAlignedResIndex(groupInChainA, chainA);
			
			Group groupInChainB = chainB.getGroupByPDB(new ResidueNumber("B", resNum+2000, null));
			int indexInChainB = chainB.getCompound().getAlignedResIndex(groupInChainB, chainB);
			
			Group groupInChainC = chainC.getGroupByPDB(new ResidueNumber("C", resNum+3000, null));
			int indexInChainC = chainC.getCompound().getAlignedResIndex(groupInChainC, chainC);

			assertNotEquals(-1, indexInChainA);
			assertNotEquals(-1, indexInChainB);
			assertNotEquals(-1, indexInChainC);
			
			assertEquals(indexInChainA,indexInChainB);
			assertEquals(indexInChainA,indexInChainC);
		}

		
		
		// this should work either with or without setAlignSeqRes, since the mapping happens in CompoundFinder
		s = getStructure("3ddo_raw_noseqres.pdb.gz", false);
		
		assertEquals(1,s.getCompounds().size());
		
		chainA = s.getChainByPDB("A");
		chainB = s.getChainByPDB("B");
		chainC = s.getChainByPDB("C");

		// the last 3 residues of each chain are all the same: they should map to the same index
		for (int resNum=251;resNum<=253;resNum++) {
			Group groupInChainA = chainA.getGroupByPDB(new ResidueNumber("A", resNum+1000, null));
			int indexInChainA = chainA.getCompound().getAlignedResIndex(groupInChainA, chainA);
			
			Group groupInChainB = chainB.getGroupByPDB(new ResidueNumber("B", resNum+2000, null));
			int indexInChainB = chainB.getCompound().getAlignedResIndex(groupInChainB, chainB);
			
			Group groupInChainC = chainC.getGroupByPDB(new ResidueNumber("C", resNum+3000, null));
			int indexInChainC = chainC.getCompound().getAlignedResIndex(groupInChainC, chainC);

			assertNotEquals(-1, indexInChainA);
			assertNotEquals(-1, indexInChainB);
			assertNotEquals(-1, indexInChainC);

			assertEquals(indexInChainA,indexInChainB);
			assertEquals(indexInChainA,indexInChainC);
		}


	}
	
	@Test
	public void test3ddoRawSeqres() throws IOException, StructureException { 

		// 3ddo has 6 chains in 1 entity, all of them with different residue numbering (chain A is 1000+, chain B 2000+ ...)
		Structure s = getStructure("3ddo_raw_trunc_seqres.pdb.gz", true);

		assertEquals(1,s.getCompounds().size());

		Chain chainA = s.getChainByPDB("A");
		Chain chainB = s.getChainByPDB("B");
		Chain chainC = s.getChainByPDB("C");

		// the last 3 residues of each chain are all the same: they should map to the same index
		for (int resNum=28;resNum<=30;resNum++) {
			Group groupInChainA = chainA.getGroupByPDB(new ResidueNumber("A", resNum+1000, null));
			int indexInChainA = chainA.getCompound().getAlignedResIndex(groupInChainA, chainA);

			Group groupInChainB = chainB.getGroupByPDB(new ResidueNumber("B", resNum+2000, null));
			int indexInChainB = chainB.getCompound().getAlignedResIndex(groupInChainB, chainB);

			Group groupInChainC = chainC.getGroupByPDB(new ResidueNumber("C", resNum+3000, null));
			int indexInChainC = chainC.getCompound().getAlignedResIndex(groupInChainC, chainC);

			assertNotEquals(-1, indexInChainA);
			assertNotEquals(-1, indexInChainB);
			assertNotEquals(-1, indexInChainC);

			assertEquals(indexInChainA,indexInChainB);
			assertEquals(indexInChainA,indexInChainC);
		}



		// this will not work without setAlignSeqRes, since the mapping happens in SeqRes2AtomAligner
		//s = getStructure("3ddo_raw_trunc_seqres.pdb.gz", false);
		



	}
	
	private Structure getStructure(String fileName, boolean setAlignSeqRes) throws IOException {
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream(PATH_TO_TEST_FILES+fileName));
		
		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(setAlignSeqRes);
		pdbpars.setFileParsingParameters(params);

		return pdbpars.parsePDBFile(inStream) ;
	}
}
