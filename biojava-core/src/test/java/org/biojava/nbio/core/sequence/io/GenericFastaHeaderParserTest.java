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
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DataSource;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.junit.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static org.junit.Assert.assertEquals;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GenericFastaHeaderParserTest {

	private final static Logger logger = LoggerFactory.getLogger(GenericFastaHeaderParserTest.class);

	public GenericFastaHeaderParserTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@AfterClass
	public static void tearDownClass() throws Exception {
	}

	@Before
	public void setUp() {
	}

	@After
	public void tearDown() {
	}

	/**
	 * GenBank                           gi|gi-number|gb|accession|locus
	 * ENA Data Library                 gi|gi-number|emb|accession|locus
	 * DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
	 * NBRF PIR                          pir||entry
	 * Protein Research Foundation       prf||name
	 * SWISS-PROT                        sp|accession|name
	 * Brookhaven Protein Data Bank (1)  pdb|entry|chain
	 * Brookhaven Protein Data Bank (2)  entry:chain|PDBID|CHAIN|SEQUENCE
	 * PDB EBI                           PDB:1ECY_A mol:protein length:142  ECOTIN
	 * Patents                           pat|country|number
	 * GenInfo Backbone Id               bbs|number
	 * General database identifier       gnl|database|identifier
	 * NCBI Reference Sequence           ref|accession|locus
	 * Local Sequence identifier         lcl|identifier
	 *
	 * @author Scooter Willis <willishf at gmail dot com>
	 */
	@Test
	public void testParseHeader() throws CompoundNotFoundException {
		logger.info("parseHeader");
		String header = "";
		ProteinSequence sequence = new ProteinSequence("");
		GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound> instance = new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>();

		header = "gi|gi-number|gb|accession|locus";
		instance.parseHeader(header, sequence);
		assertEquals("accession", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.GENBANK);

		header = "gi|gi-number|emb|accession|locus";
		instance.parseHeader(header, sequence);
		assertEquals("accession", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.ENA);

		header = "gi|gi-number|dbj|accession|locus";
		instance.parseHeader(header, sequence);
		assertEquals("accession", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.DDBJ);

		header = "pir||entry";
		instance.parseHeader(header, sequence);
		assertEquals("entry", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.NBRF);

		header = "prf||name";
		instance.parseHeader(header, sequence);
		assertEquals("name", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.PRF);

		header = "sp|accession|name";
		instance.parseHeader(header, sequence);
		assertEquals("accession", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.UNIPROT);

		header = "pdb|entry|chain";
		instance.parseHeader(header, sequence);
		assertEquals("entry:chain", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.PDB1);

		header = "entry:chain|PDBID|CHAIN|SEQUENCE";
		instance.parseHeader(header, sequence);
		assertEquals("entry:chain", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.PDB2);
		header = "PDB:1ECY_A mol:protein length:142  ECOTIN";
		instance.parseHeader(header, sequence);
		assertEquals("1ECY_A", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.PDBe);

		header = "pat|country|number";
		instance.parseHeader(header, sequence);
		assertEquals("number", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.PATENTS);

		header = "bbs|number";
		instance.parseHeader(header, sequence);
		assertEquals("number", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.GENINFO);

		header = "gnl|database|identifier";
		instance.parseHeader(header, sequence);
		assertEquals("identifier", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.GENERAL);

		header = "ref|accession|locus";

		instance.parseHeader(header, sequence);
		assertEquals("accession", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.NCBI);

		header = "lcl|identifier";
		instance.parseHeader(header, sequence);
		assertEquals("identifier", sequence.getAccession().getID());
		assertEquals(sequence.getAccession().getDataSource(), DataSource.LOCAL);
		// TODO review the generated test code and remove the default call to fail.
		//fail("The test case is a prototype.");
	}
}
