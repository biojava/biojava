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

import static org.junit.Assert.*;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.junit.Ignore;
import org.junit.Test;

/**
 * Various tests for functionality in {@link EntityInfo} and {@link org.biojava.nbio.structure.io.EntityFinder}
 * @author Jose Duarte
 */
public class TestEntityResIndexMapping {

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

		Chain chainA = s.getPolyChainByPDB("A");
		int i = chainA.getEntityInfo().getAlignedResIndex(chainA.getAtomGroup(0),chainA);
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

		Chain chainA = s.getPolyChainByPDB("A");
		int i = chainA.getEntityInfo().getAlignedResIndex(chainA.getAtomGroup(0),chainA);
		assertEquals("First residue in 1smtA "+chainA.getAtomGroup(0).toString()+" should map to 24 in SEQRES",24,i);
		Chain chainB = s.getPolyChainByPDB("B");
		i = chainB.getEntityInfo().getAlignedResIndex(chainB.getAtomGroup(0),chainB);
		assertEquals("First residue in 1smtB "+chainB.getAtomGroup(0).toString()+" should map to 20 in SEQRES",20,i);

		// group with seqres index 19 is observed in chain B but not in chain A, we should still get the index back from getAlignedResIndex
		i = chainA.getEntityInfo().getAlignedResIndex(chainA.getSeqResGroup(19),chainA);
		assertEquals("Seqres residue 20 in 1smtA "+chainA.getSeqResGroup(19).toString()+" should map to 20 in SEQRES",20,i);

		checkAllResidues(s);
	}

	private void checkAllResidues(Structure s) {
		for (Chain c:s.getChains()) {
			for (Group g:c.getAtomGroups()) {
				if (g.getType() != GroupType.AMINOACID) continue;

				int resIndex = c.getEntityInfo().getAlignedResIndex(g, c);
				assertNotEquals("Residue "+g.getResidueNumber()+" should not map to -1 in SEQRES",-1, resIndex);
				assertNotEquals("Residue "+g.getResidueNumber()+" should not map to 0 in SEQRES",0, resIndex);
				assertTrue(resIndex <= c.getSeqResLength());
			}
		}
	}

	// This doesn't work yet, since for raw files without a SEQRES, the seqres groups are not populated. Instead
	// in that case EntityInfo.getAlignedResIndex() returns residue numbers as given (without insertion codes) and
	// thus in general residues will not be correctly aligned between different chains of same entity. This breaks
	// cases like 3ddo (with no SEQRES records) where residue numbering is different in every chain of the one entity.
	// see https://github.com/eppic-team/eppic/issues/39
	@Ignore
	@Test
	public void test3ddoRawNoSeqres() throws IOException, StructureException {

		// 3ddo has 6 chains in 1 entity, all of them with different residue numbering (chain A is 1000+, chain B 2000+ ...)
		Structure s = getStructure("3ddo_raw_noseqres.pdb.gz", true);

		List<EntityInfo> polyEntities = new ArrayList<>();
		for (EntityInfo entityInfo : s.getEntityInfos()) {
			if (entityInfo.getType() == EntityType.POLYMER) {
				polyEntities.add(entityInfo);
			}
		}

		assertEquals(1, polyEntities.size());

		Chain chainA = s.getPolyChainByPDB("A");
		Chain chainB = s.getPolyChainByPDB("B");
		Chain chainC = s.getPolyChainByPDB("C");

		// the last 3 residues of each chain are all the same: they should map to the same index
		for (int resNum=251;resNum<=253;resNum++) {
			Group groupInChainA = chainA.getGroupByPDB(new ResidueNumber("A", resNum+1000, null));
			int indexInChainA = chainA.getEntityInfo().getAlignedResIndex(groupInChainA, chainA);

			Group groupInChainB = chainB.getGroupByPDB(new ResidueNumber("B", resNum+2000, null));
			int indexInChainB = chainB.getEntityInfo().getAlignedResIndex(groupInChainB, chainB);

			Group groupInChainC = chainC.getGroupByPDB(new ResidueNumber("C", resNum+3000, null));
			int indexInChainC = chainC.getEntityInfo().getAlignedResIndex(groupInChainC, chainC);

			assertNotEquals(-1, indexInChainA);
			assertNotEquals(-1, indexInChainB);
			assertNotEquals(-1, indexInChainC);

			assertEquals(indexInChainA,indexInChainB);
			assertEquals(indexInChainA,indexInChainC);
		}



		// this should work either with or without setAlignSeqRes, since the mapping happens in EntityFinder
		s = getStructure("3ddo_raw_noseqres.pdb.gz", false);

		polyEntities = new ArrayList<>();
		for (EntityInfo entityInfo : s.getEntityInfos()) {
			if (entityInfo.getType() == EntityType.POLYMER) {
				polyEntities.add(entityInfo);
			}
		}

		assertEquals(1, polyEntities.size());

		chainA = s.getPolyChainByPDB("A");
		chainB = s.getPolyChainByPDB("B");
		chainC = s.getPolyChainByPDB("C");

		// the last 3 residues of each chain are all the same: they should map to the same index
		for (int resNum=251;resNum<=253;resNum++) {
			Group groupInChainA = chainA.getGroupByPDB(new ResidueNumber("A", resNum+1000, null));
			int indexInChainA = chainA.getEntityInfo().getAlignedResIndex(groupInChainA, chainA);

			Group groupInChainB = chainB.getGroupByPDB(new ResidueNumber("B", resNum+2000, null));
			int indexInChainB = chainB.getEntityInfo().getAlignedResIndex(groupInChainB, chainB);

			Group groupInChainC = chainC.getGroupByPDB(new ResidueNumber("C", resNum+3000, null));
			int indexInChainC = chainC.getEntityInfo().getAlignedResIndex(groupInChainC, chainC);

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

		assertEquals(1,s.getEntityInfos().size());

		Chain chainA = s.getPolyChainByPDB("A");
		Chain chainB = s.getPolyChainByPDB("B");
		Chain chainC = s.getPolyChainByPDB("C");

		// the last 3 residues of each chain are all the same: they should map to the same index
		for (int resNum=28;resNum<=30;resNum++) {
			Group groupInChainA = chainA.getGroupByPDB(new ResidueNumber("A", resNum+1000, null));
			int indexInChainA = chainA.getEntityInfo().getAlignedResIndex(groupInChainA, chainA);

			Group groupInChainB = chainB.getGroupByPDB(new ResidueNumber("B", resNum+2000, null));
			int indexInChainB = chainB.getEntityInfo().getAlignedResIndex(groupInChainB, chainB);

			Group groupInChainC = chainC.getGroupByPDB(new ResidueNumber("C", resNum+3000, null));
			int indexInChainC = chainC.getEntityInfo().getAlignedResIndex(groupInChainC, chainC);

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
