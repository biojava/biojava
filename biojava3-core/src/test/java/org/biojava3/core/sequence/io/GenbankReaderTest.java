/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import static org.junit.Assert.assertNotNull;

import java.io.InputStream;
import java.util.LinkedHashMap;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GenbankReaderTest {

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

		System.out.println("process protein");
		InputStream inStream = this.getClass().getResourceAsStream("/BondFeature.gb");
		assertNotNull(inStream);
		
		GenbankReader<ProteinSequence,AminoAcidCompound> GenbankProtein = 
				new GenbankReader<ProteinSequence,AminoAcidCompound>(
						inStream, 
						new GenericGenbankHeaderParser<ProteinSequence,AminoAcidCompound>(), 
						new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet())
						);
		LinkedHashMap<String,ProteinSequence> proteinSequences = GenbankProtein.process();
		inStream.close();

		System.out.println("process DNA");
		inStream = this.getClass().getResourceAsStream("/NM_000266.gb");
		assertNotNull(inStream);

		GenbankReader<DNASequence,NucleotideCompound> GenbankDNA = 
				new GenbankReader<DNASequence,NucleotideCompound>(
						inStream,
						new GenericGenbankHeaderParser<DNASequence,NucleotideCompound>(), 
						new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
						);
		LinkedHashMap<String,DNASequence> dnaSequences = GenbankDNA.process();
		inStream.close();
	}
	
}
