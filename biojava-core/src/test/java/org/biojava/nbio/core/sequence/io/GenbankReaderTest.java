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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.junit.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.InputStream;
import java.util.LinkedHashMap;

import static org.junit.Assert.assertNotNull;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GenbankReaderTest {

	private final static Logger logger = LoggerFactory.getLogger(GenbankReaderTest.class);

	public GenbankReaderTest() {
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
	 * Test of process method, of class GenbankReader.
	 */
	@Test
	public void testProcess() throws Exception {

		logger.info("process protein");
		InputStream inStream = this.getClass().getResourceAsStream("/BondFeature.gb");
		assertNotNull(inStream);
		
		GenbankReader<ProteinSequence,AminoAcidCompound> genbankProtein = 
				new GenbankReader<ProteinSequence,AminoAcidCompound>(
						inStream, 
						new GenericGenbankHeaderParser<ProteinSequence,AminoAcidCompound>(), 
						new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet())
						);
		@SuppressWarnings("unused")
		LinkedHashMap<String,ProteinSequence> proteinSequences = genbankProtein.process();
		inStream.close();

		logger.info("process DNA");
		inStream = this.getClass().getResourceAsStream("/NM_000266.gb");
		assertNotNull(inStream);

		GenbankReader<DNASequence,NucleotideCompound> genbankDNA = 
				new GenbankReader<DNASequence,NucleotideCompound>(
						inStream,
						new GenericGenbankHeaderParser<DNASequence,NucleotideCompound>(), 
						new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
						);
		@SuppressWarnings("unused")
		LinkedHashMap<String,DNASequence> oneDNASequence = genbankDNA.process();
		inStream.close();

		logger.info("process 2 DNAs");
		inStream = this.getClass().getResourceAsStream("/two-dnaseqs.gb");
		assertNotNull(inStream);

		GenbankReader<DNASequence,NucleotideCompound> readTwoDNAs = 
				new GenbankReader<DNASequence,NucleotideCompound>(
						inStream,
						new GenericGenbankHeaderParser<DNASequence,NucleotideCompound>(), 
						new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
						);
		@SuppressWarnings("unused")
		LinkedHashMap<String,DNASequence> twoDNAs = readTwoDNAs.process();
		inStream.close();

		logger.info("process coli genome");
		inStream = this.getClass().getResourceAsStream("/NC_000913.gb");
		assertNotNull(inStream);

		
		GenbankReader<DNASequence,NucleotideCompound> readColiGenome = 
				new GenbankReader<DNASequence,NucleotideCompound>(
						inStream,
						new GenericGenbankHeaderParser<DNASequence,NucleotideCompound>(), 
						new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
						);
		@SuppressWarnings("unused")
		LinkedHashMap<String,DNASequence> coliGenome = readColiGenome.process();
		assertNotNull(inStream);
		inStream.close();	
	}
	
}
