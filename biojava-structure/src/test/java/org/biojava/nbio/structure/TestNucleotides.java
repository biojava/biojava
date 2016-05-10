/**
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
 * Created on Apr 20, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.structure;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.ChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.*;


/** This class tests the correct loading of Nucleotides
 *
 * @author Andreas Prlic
 * @since 3.0.3
 */
public class TestNucleotides {

	private static AtomCache cache;

	@BeforeClass
	public static void beforeClass() {
		cache = new AtomCache();
	}

	@Test
	public void test3T5N() throws IOException, StructureException{

		String pdbId = "3T5N";
		Structure s = getStructure(pdbId);


		assertEquals(2,s.getPolyChains().size());

		Chain c = s.getChains().get(1);
		System.out.println(c);
		assertEquals("C", c.getName());
		List<Group> ngr = c.getAtomGroups(GroupType.NUCLEOTIDE);
		assertEquals(6,ngr.size());


		// now test if we download all definitions correctly for this one...
		PDBFileReader reader = new PDBFileReader();
		FileParsingParameters params = new FileParsingParameters();
		params.setParseSecStruc(true);
		params.setAlignSeqRes(true);
		params.setParseCAOnly(false);
		reader.setFileParsingParameters(params);

		ChemCompProvider chemProv = ChemCompGroupFactory.getChemCompProvider();

		DownloadChemCompProvider download = new DownloadChemCompProvider();

		ChemCompGroupFactory.setChemCompProvider(download);

		Structure s1 = reader.getStructureById(pdbId);

		assertNotNull(s1);

		assertEquals(2,s1.getPolyChains().size());

		Chain c1 = s1.getChains().get(1);

		assertEquals("C", c1.getName());

		Group g = c1.getAtomGroup(0);
		assertNotNull(g);
		assertNotNull(g.getChemComp());
		assertNotNull(g.getChemComp().getPolymerType());
		assertNotNull(g.getChemComp().getPolymerType().name());

		assertTrue("Found an unknown polymertype!", (! g.getChemComp().getPolymerType().equals(PolymerType.unknown)));
		//System.out.println(g.getChemComp().getPolymerType());

		List<Group> ngr1 = c1.getAtomGroups(GroupType.NUCLEOTIDE);
		assertEquals(6,ngr1.size());


		ChemCompGroupFactory.setChemCompProvider(chemProv);


	}

	@Test
	public void test1OFX() throws StructureException, IOException {
		Structure s = getStructure("1OFX");

		assertEquals(2,s.getPolyChains().size());

		Chain a = s.getChains().get(0);
		assertEquals("A", a.getId());
		List<Group> ngrA = a.getAtomGroups(GroupType.NUCLEOTIDE);
		assertEquals(10,ngrA.size());

		Chain b = s.getChains().get(1);
		assertEquals("B", b.getId());
		List<Group> ngrB = b.getAtomGroups(GroupType.NUCLEOTIDE);
		assertEquals(10,ngrB.size());
	}

	private Structure getStructure(String pdbId) throws IOException, StructureException {
		//System.out.println("cache: " + ChemCompGroupFactory.getChemCompProvider().getClass().getName());

		//System.out.println("cache: download chem comp:" + cache.getFileParsingParams().isLoadChemCompInfo());
		return cache.getStructure(pdbId);
	}

	@Test
	public void test1REP() throws StructureException, IOException{

		PDBFileReader reader = new PDBFileReader();
		FileParsingParameters params = new FileParsingParameters();
		params.setParseSecStruc(true);
		params.setAlignSeqRes(true);
		params.setParseCAOnly(false);
		reader.setFileParsingParameters(params);



		Structure s = reader.getStructureById("1REP");
		//System.out.println(s);
		//System.out.println(s.toPDB());
		Chain b = s.getPolyChainByPDB("B");

		assertEquals(22,b.getSeqResGroups().size());
		assertEquals(21,b.getAtomGroups().size());

		Group n1 = b.getSeqResGroup(0);
		Group n2 = b.getAtomGroup(0);
		//System.out.println(n1);
		//System.out.println(n2);
		//System.out.println(n1.getChemComp());


		assertNotNull("Could not acces Chem Comp file!" , n1.getChemComp());
		assertTrue("ChemComp is not DC",n1.getChemComp().getId().equals("DC"));
		assertNotNull("Could not determine polymer type " , n1.getChemComp().getPolymerType());
		//System.out.println(n1.getChemComp().getPolymerType());
		assertTrue(n1.getChemComp().getPolymerType().equals(PolymerType.dna));

		assertNotNull(n1.getPDBName());
		assertNotNull(n1.getResidueNumber());
		assertNotNull(n2.getResidueNumber());
		assertEquals("23", n2.getResidueNumber().toString());
		assertTrue(n1.getResidueNumber().equals(n2.getResidueNumber()));


	}
}
