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
 * Created on 2012-11-20
 * Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.bio.structure.rcsb;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.List;

import org.junit.Test;

/**
 * Tests {@link RCSBDescriptionFactory}.
 * @author dmyerstu
 * @since 3.0.6
 */
public class RCSBDescriptionFactoryTest {

	private static final String TEST_DIR = "src/test/resources/";
	
	/**
	 * Opens the file as a {@link FileInputStream}. Copied from ResourceList, which is not in biojava.
	 */
	private FileInputStream openStream(String filename) {
		File file = new File(TEST_DIR + filename);
		FileInputStream fis;
		try {
			fis = new FileInputStream(file);
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}
		return fis;
	}

	/**
	 * Tests on the live database. Just makes sure the resource can be found.
	 * If this test fails, it may be because the database changed.
	 */
	@Test
	public void testUrl() {
		RCSBDescriptionFactory.get("1w0p"); // just make sure it doesn't throw an exception
	}
	
	/**
	 * Covers all the basic features, including EC numbers. Does not cover multiple polymers or multiple chains.
	 */
	@Test
	public void test1() {
		RCSBDescription description = RCSBDescriptionFactory.get(openStream("describeMol/1w0p.xml"));
		
		assertEquals("1W0P", description.getPdbId());
		List<RCSBPolymer> polymers = description.getPolymers();
		assertEquals(1, polymers.size());
		
		RCSBPolymer polymer = polymers.get(0);
		assertEquals("protein", polymer.getType());
		assertEquals(1, polymer.getIndex().intValue());
		assertEquals("SIALIDASE", polymer.getDescription());
		assertEquals("3.2.1.18", polymer.getEnzClass());
		assertEquals(781, polymer.getLength().intValue());
		assertEquals(85675.5, polymer.getWeight(), 0);
		
		List<Character> chains = polymer.getChains();
		assertEquals(1, chains.size());
		assertEquals('A', (char) chains.get(0));
		
		List<String> synonyms = polymer.getSynonyms();
		assertEquals(2, synonyms.size());
		assertEquals("NEURAMINIDASE", synonyms.get(0));
		assertEquals("NANASE", synonyms.get(1));
		
		RCSBTaxonomy tax = polymer.getTaxonomy();
		assertEquals(666, tax.getId());
		assertEquals("Vibrio cholerae", tax.getName());
		
		RCSBMacromolecule mol = polymer.getMolecule();
		assertEquals("Sialidase", mol.getName());
		List<String> accessions = mol.getAccessions();
		assertEquals(4, accessions.size());
		assertEquals("A5F7A4", accessions.get(0));
		assertEquals("C3M1H8", accessions.get(1));
		assertEquals("P37060", accessions.get(2));
		assertEquals("Q9KR59", accessions.get(3));
	}

	/**
	 * What if we have a structureId but no polymers?
	 */
	@Test
	public void testEmpty() {
		RCSBDescription description = RCSBDescriptionFactory.get(openStream("describeMol/empty.xml"));
		assertEquals("empty", description.getPdbId());
		List<RCSBPolymer> polymers = description.getPolymers();
		assertEquals(0, polymers.size());
	}

	/**
	 * What if we have polymers but no macroMolecule or chains?
	 * And what if a polymer contains no attributes?
	 */
	@Test
	public void testAlmostEmpty() {
		
		RCSBDescription description = RCSBDescriptionFactory.get(openStream("describeMol/almost_empty.xml"));
		assertEquals("almost_empty", description.getPdbId());
		List<RCSBPolymer> polymers = description.getPolymers();
		assertEquals(2, polymers.size());

		RCSBPolymer polymer = polymers.get(0);
		assertEquals("notype", polymer.getType());
		assertEquals(1, polymer.getIndex().intValue());
		assertEquals("really close to empty", polymer.getDescription());
		assertEquals(null, polymer.getEnzClass());
		assertEquals(10, polymer.getLength().intValue());
		assertEquals(0, polymer.getWeight(), 0);

		polymer = polymers.get(1);
		assertEquals(null, polymer.getType()); // make sure these are null and not ""
		assertEquals(null, polymer.getIndex());
		assertEquals(null, polymer.getDescription());
		assertEquals(null, polymer.getEnzClass());
		assertEquals(null, polymer.getLength());
		assertEquals(null, polymer.getWeight());
		
	}
	
	/**
	 * Covers multiple polymers and multiple chains.
	 */
	@Test
	public void test2() {
		RCSBDescription description = RCSBDescriptionFactory.get(openStream("describeMol/4hhb.xml"));
		assertEquals("4HHB", description.getPdbId());
		List<RCSBPolymer> polymers = description.getPolymers();
		assertEquals(2, polymers.size());
		
		// first polymer
		RCSBPolymer polymer = polymers.get(0);
		assertEquals("protein", polymer.getType());
		assertEquals(1, polymer.getIndex().intValue());
		assertEquals("HEMOGLOBIN (DEOXY) (ALPHA CHAIN)", polymer.getDescription());
		assertEquals(null, polymer.getEnzClass());
		assertEquals(141, polymer.getLength().intValue());
		assertEquals(15150.5, polymer.getWeight(), 0);
		
		List<Character> chains = polymer.getChains();
		assertEquals(2, chains.size());
		assertEquals('A', (char) chains.get(0));
		assertEquals('C', (char) chains.get(1));
		
		List<String> synonyms = polymer.getSynonyms();
		assertEquals(0, synonyms.size());
		
		RCSBTaxonomy tax = polymer.getTaxonomy();
		assertEquals(9606, tax.getId());
		assertEquals("Homo sapiens", tax.getName());
		
		RCSBMacromolecule mol = polymer.getMolecule();
		assertEquals("Hemoglobin subunit alpha", mol.getName());
		List<String> accessions = mol.getAccessions();
		assertEquals(8, accessions.size());
		assertEquals("P69905", accessions.get(0));
		assertEquals("P01922", accessions.get(1));
		assertEquals("Q1HDT5", accessions.get(2));
		assertEquals("Q3MIF5", accessions.get(3));
		assertEquals("Q53F97", accessions.get(4));
		assertEquals("Q96KF1", accessions.get(5));
		assertEquals("Q9NYR7", accessions.get(6));
		assertEquals("Q9UCM0", accessions.get(7));
		
		// second polymer
		polymer = polymers.get(1);
		assertEquals("protein", polymer.getType());
		assertEquals(2, polymer.getIndex().intValue());
		assertEquals("HEMOGLOBIN (DEOXY) (BETA CHAIN)", polymer.getDescription());
		assertEquals(null, polymer.getEnzClass());
		assertEquals(146, polymer.getLength().intValue());
		assertEquals(15890.4, polymer.getWeight(), 0);
		
		chains = polymer.getChains();
		assertEquals(2, chains.size());
		assertEquals('B', (char) chains.get(0));
		assertEquals('D', (char) chains.get(1));
		
		synonyms = polymer.getSynonyms();
		assertEquals(0, synonyms.size());
		
		tax = polymer.getTaxonomy();
		assertEquals(9606, tax.getId());
		assertEquals("Homo sapiens", tax.getName());
		
		mol = polymer.getMolecule();
		assertEquals("Hemoglobin subunit beta", mol.getName());
		accessions = mol.getAccessions();
		assertEquals(16, accessions.size());
		assertEquals("P68871", accessions.get(0));
		assertEquals("A4GX73", accessions.get(1));
		assertEquals("B2ZUE0", accessions.get(2));
		assertEquals("P02023", accessions.get(3));
		assertEquals("Q13852", accessions.get(4));
		assertEquals("Q14481", accessions.get(5));
		assertEquals("Q14510", accessions.get(6));
		assertEquals("Q45KT0", accessions.get(7));
		assertEquals("Q549N7", accessions.get(8));
		assertEquals("Q6FI08", accessions.get(9));
		assertEquals("Q6R7N2", accessions.get(10));
		assertEquals("Q8IZI1", accessions.get(11));
		assertEquals("Q9BX96", accessions.get(12));
		assertEquals("Q9UCD6", accessions.get(13));
		assertEquals("Q9UCP8", accessions.get(14));
		assertEquals("Q9UCP9", accessions.get(15));
		
	}

}
