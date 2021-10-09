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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assume.assumeNotNull;
import static org.junit.Assume.assumeTrue;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.cif.CifStructureConverter;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.junit.Test;

/**
 * Testing parsing of some difficult mmCIF files.
 * For instance those containing multi-line quoting using ";\n" as delimiters
 * Feel free to add any other difficult case here
 *
 *
 * @author duarte_j
 *
 */
public class TestDifficultMmCIFFiles {

	/**
	 * The 2KSA mmCIF contains a 5 model NMR structure.  The first residue of the sequence is not visible
	 * and the models should all begin with Asp indexed as residue #2.
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void test2KSA() throws IOException, StructureException {
		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		FileParsingParameters params = cache.getFileParsingParams();
		params.setParseBioAssembly(true);
		params.setAlignSeqRes(true);
		StructureIO.setAtomCache(cache);


		cache.setFiletype(StructureFiletype.CIF);
		Structure sCif = StructureIO.getStructure("2KSA");

		assertNotNull(sCif);

		// Unit test for each of the chains to show they begin with the correct first residue.
		for (int i = 0; i < sCif.nrModels(); i++) {
			List<Chain> chains = sCif.getModel(i);

			// Chain A first residue should start at ASP 2..
			// but if replaceGroupSeqPos(PdbxPolySeqScheme ppss) is used, this is incorrect and will be 1.
			assertEquals(2, chains.get(0).getAtomGroup(0).getResidueNumber().getSeqNum().intValue());
		}
	}

	@Test
	public void test2BI6() throws IOException, StructureException {

		// In this entry _struct_conf contains multiline quoting (quoting with "\n;" ) in a non-loop field

		// It seems that at the moment the field is not parsed by the mmCIF parser, anyway let's
		// keep this here if in the future it is


		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setFiletype(StructureFiletype.CIF);
		Structure sCif = StructureIO.getStructure("2BI6");

		assertNotNull(sCif);

		// an NMR entry
		assertFalse(sCif.isCrystallographic());

		assertTrue(sCif.isNmr());

		// no more revision records since v5 - JD 2017-07-13
		//assertTrue(sCif.getPDBHeader().getRevisionRecords().size() > 1);


	}

	@Test
	public void test1GQO() throws IOException, StructureException {

		// In this entry _pdbx_struct_assembly_gen contains multiline quoting (quoting with "\n;" ) in loop field

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		FileParsingParameters params = cache.getFileParsingParams();
		params.setParseBioAssembly(true);
		StructureIO.setAtomCache(cache);

		cache.setFiletype(StructureFiletype.PDB);
		Structure sPdb = StructureIO.getStructure("1GQO");

		cache.setFiletype(StructureFiletype.CIF);
		Structure sCif = StructureIO.getStructure("1GQO");

		assertNotNull(sCif);

		assertNotNull(sPdb.getPDBHeader().getBioAssemblies());
		assertNotNull(sCif.getPDBHeader().getBioAssemblies());

		Map<Integer,BioAssemblyInfo> mapPdb = sPdb.getPDBHeader().getBioAssemblies();
		Map<Integer,BioAssemblyInfo> mapCif = sCif.getPDBHeader().getBioAssemblies();



		assertEquals(mapPdb.size(),mapCif.size());

		assertEquals(60, mapCif.get(1).getTransforms().size());
		assertEquals(60, mapCif.get(2).getTransforms().size());

		// an X-RAY entry
		assertTrue(sPdb.isCrystallographic());
		assertTrue(sCif.isCrystallographic());

		assertFalse(sPdb.isNmr());
		assertFalse(sCif.isNmr());


	}

	@Test
	public void testResidueNumbers() throws IOException, StructureException {
		AtomCache cache = new AtomCache();
		cache.setFiletype(StructureFiletype.CIF);

		Structure s = cache.getStructure("2PTC");
		Chain c = s.getChainByIndex(0);
		System.out.println(c);
		assertEquals("Wrong first chain",c.getName(),"E");

		Group res = c.getAtomGroup(0);
		ResidueNumber resNum = res.getResidueNumber();

		assertEquals("Groups have wrong chain in resnum",resNum.getChainName(),"E");
	}

	@Test
	public void test4letterChains() throws IOException, StructureException, URISyntaxException {
		String filename = "/1hh0_4char.cif.gz";
		URL url = getClass().getResource(filename);
		assumeNotNull("Can't find resource "+filename,url);

		File file = new File(url.toURI());
		assumeNotNull(file);
		assumeTrue(file.exists());

		CifFileReader reader = new CifFileReader();
		Structure s = reader.getStructure(file);

		assertNotNull("Failed to load structure from jar",s);

		List<Chain> chains = s.getChains();
		assertEquals("Wrong number of chains",chains.size(), 1);

		Chain chain = chains.get(0);
		assertEquals("Wrong chain ID",chain.getId(),"ABCD");

		Chain chain2 = s.getPolyChainByPDB("ABCD");
		assertNotNull(chain2);
		assertEquals(chain2, chain);
	}

	/**
	 * This is to test the issue discussed here:
	 * http://www.globalphasing.com/startools/
	 * Essentially single quote characters (') are valid not only for quoting, but also as parts of
	 * data values as long as some rules of the STAR format are followed.
	 * For instance Phenix produces mmCIF files with non-quoted strings containing single quote characters
	 * @throws IOException
	 */
	@Test
	public void testQuotingCornerCase () throws IOException {
		InputStream inStream = this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/difficult_mmcif_quoting.cif");

		FileParsingParameters fileParsingParams = new FileParsingParameters();
		fileParsingParams.setAlignSeqRes(true);
		Structure s = CifStructureConverter.fromInputStream(inStream, fileParsingParams);

		assertNotNull(s);


	}

	/**
	 * The last category in 2KLI mmCIF file is _pdbx_struct_oper_list, which is needed for
	 * the biounit annotation.
	 * This tests makes sure that the last category in a mmCIF file is not missed because
	 * of its position as last one in file.
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void test2KLI() throws IOException, StructureException {

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		FileParsingParameters params = cache.getFileParsingParams();
		params.setParseBioAssembly(true);
		StructureIO.setAtomCache(cache);


		cache.setFiletype(StructureFiletype.CIF);
		Structure sCif = StructureIO.getStructure("2KLI");

		assertNotNull(sCif);

		assertNotNull(sCif.getPDBHeader().getBioAssemblies());

		Map<Integer,BioAssemblyInfo> mapCif = sCif.getPDBHeader().getBioAssemblies();

		assertNotNull(mapCif);

	}
}
