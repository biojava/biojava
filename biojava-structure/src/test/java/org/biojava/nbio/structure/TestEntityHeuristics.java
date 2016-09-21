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

import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.junit.Test;

import java.io.IOException;
import java.io.InputStream;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.*;

public class TestEntityHeuristics {

	private static final String PATH_TO_TEST_FILES = "/org/biojava/nbio/structure/io/";

	//private static final String[] FILES_TO_TEST =
	//	{"1b8g_raw.pdb.gz", "2m7y_raw.pdb.gz", "3c5f_raw.pdb.gz", "4b19_raw.pdb.gz"};



	//public void testCompoundHeuristics() throws IOException {
	//	for (String fileName:FILES_TO_TEST) {
	//		testFile(fileName);
	//	}
	//}

	@Test
	public void test1b8gRaw() throws IOException, StructureException {

		Structure s = getStructure("1b8g_raw.pdb.gz", true);

		// 1 prot, 1 PLP, 1 water: 3 entities
		assertEquals(3, s.getEntityInfos().size());

		Chain chainA = s.getPolyChainByPDB("A");

		assertEquals(2, chainA.getEntityInfo().getChains().size());

		assertEquals(chainA,chainA.getEntityInfo().getRepresentative());

		checkEntitiesNumbered(s.getEntityInfos());


		s = getStructure("1b8g_raw.pdb.gz", false);

		// 1 prot, 1 PLP, 1 water: 3 entities
		assertEquals(3, s.getEntityInfos().size());

		chainA = s.getPolyChainByPDB("A");

		assertEquals(2, chainA.getEntityInfo().getChains().size());

		assertEquals(chainA,chainA.getEntityInfo().getRepresentative());

		checkEntitiesNumbered(s.getEntityInfos());
	}


	@Test
	public void test2m7yRaw() throws IOException, StructureException {

		Structure s = getStructure("2m7y_raw.pdb.gz", true);

		// 2 entities: 1 polymer (protein) 1 non-polymer (ZN)
		assertEquals(2, s.getEntityInfos().size());

		Chain chainA = s.getPolyChainByPDB("A");

		assertEquals(1, chainA.getEntityInfo().getChains().size());

		assertEquals(chainA,chainA.getEntityInfo().getRepresentative());

		checkEntitiesNumbered(s.getEntityInfos());


		s = getStructure("2m7y_raw.pdb.gz", false);

		// 2 entities: 1 polymer (protein) 1 non-polymer (ZN)
		assertEquals(2,s.getEntityInfos().size());

		chainA = s.getPolyChainByPDB("A");

		assertEquals(1, chainA.getEntityInfo().getChains().size());

		assertEquals(chainA,chainA.getEntityInfo().getRepresentative());

		checkEntitiesNumbered(s.getEntityInfos());
	}

	@Test
	public void test3c5fRaw() throws IOException, StructureException {

		Structure s = getStructure("3c5f_raw.pdb.gz", true);

		int polyEntities = 0;
		for (EntityInfo e:s.getEntityInfos()) {
			if (e.getType()==EntityType.POLYMER) polyEntities++;
		}
		
		assertEquals(4, polyEntities);

		Chain chainA = s.getPolyChainByPDB("A");

		// there's 2 models in file, thus for the protien polymeric entity there's 4 chains (2 from each model)
		assertEquals(4, chainA.getEntityInfo().getChains().size());

		assertEquals(chainA,chainA.getEntityInfo().getRepresentative());

		checkEntitiesNumbered(s.getEntityInfos());


		s = getStructure("3c5f_raw.pdb.gz", false);

		polyEntities = 0;
		for (EntityInfo e:s.getEntityInfos()) {
			if (e.getType()==EntityType.POLYMER) polyEntities++;
		}
		
		assertEquals(4,polyEntities);

		chainA = s.getPolyChainByPDB("A");

		// there's 2 models in file, thus for the protien polymeric entity there's 4 chains (2 from each model)
		assertEquals(4, chainA.getEntityInfo().getChains().size());

		assertEquals(chainA,chainA.getEntityInfo().getRepresentative());

		checkEntitiesNumbered(s.getEntityInfos());
	}

	@Test
	public void test4b19Raw() throws IOException, StructureException {

		Structure s = getStructure("4b19_raw.pdb.gz", true);

		assertEquals(1,s.getEntityInfos().size());

		Chain chainA = s.getPolyChainByPDB("A");

		// only 1 protein entity with 1 chain, but 5 models
		assertEquals(5, chainA.getEntityInfo().getChains().size());

		assertEquals(chainA,chainA.getEntityInfo().getRepresentative());

		checkEntitiesNumbered(s.getEntityInfos());


		s = getStructure("4b19_raw.pdb.gz", false);

		assertEquals(1,s.getEntityInfos().size());

		chainA = s.getPolyChainByPDB("A");

		assertEquals(5, chainA.getEntityInfo().getChains().size());

		assertEquals(chainA,chainA.getEntityInfo().getRepresentative());

		checkEntitiesNumbered(s.getEntityInfos());

	}

	@Test
	public void test3ddoNoseqres() throws IOException, StructureException {

		// 3ddo contains 6 chains in 1 entity, with residue numbering completely different in each of the chains

		Structure s = getStructure("3ddo_raw_noseqres.pdb.gz", true);

		assertNotNull(s);

		assertEquals(6, s.getPolyChains().size());

		// checking that heuristics in CompoundFinder work. We should have a 2 entities: 1 polymeric (prot) and 1 nonpolymeric
		assertEquals(2, s.getEntityInfos().size());

		// trying without seqAlignSeqRes
		s = getStructure("3ddo_raw_noseqres.pdb.gz", false);
		assertNotNull(s);

		assertEquals(6, s.getPolyChains().size());

		assertEquals(2, s.getEntityInfos().size());
	}

	@Test
	public void test3ddoSeqres() throws IOException, StructureException {

		// 3ddo contains 6 chains in 1 entity, with residue numbering completely different in each of the chains

		Structure s = getStructure("3ddo_raw_trunc_seqres.pdb.gz", true);

		assertNotNull(s);

		assertEquals(6, s.getChains().size());

		// checking that heuristics in CompoundFinder work. We should have a single entity (compound)
		assertEquals(1, s.getEntityInfos().size());

		// trying without seqAlignSeqRes
		s = getStructure("3ddo_raw_trunc_seqres.pdb.gz", true);
		assertNotNull(s);

		assertEquals(6, s.getChains().size());

		assertEquals(1, s.getEntityInfos().size());
	}


	private Structure getStructure(String fileName, boolean setAlignSeqRes) throws IOException {
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream(PATH_TO_TEST_FILES+fileName));

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(setAlignSeqRes);
		pdbpars.setFileParsingParameters(params);

		Structure s = pdbpars.parsePDBFile(inStream) ;

		return s;
	}

	private void checkEntitiesNumbered(List<EntityInfo> entities) {

		Collections.sort(entities, new Comparator<EntityInfo>() {

			@Override
			public int compare(EntityInfo o1, EntityInfo o2) {
				return new Integer(o1.getMolId()).compareTo(o2.getMolId());
			}
		});

		int id = 1;
		for (EntityInfo compound:entities) {
			assertEquals(id,compound.getMolId());
			id++;
		}


	}
}
