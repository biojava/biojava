/**
 * 
 */
package org.biojava3.core.sequence.io;

import static org.junit.Assert.assertEquals;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.biojava3.core.sequence.DNASequence;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * @author mckeee1
 * 
 */
public class GenbankWriterTest {

	public GenbankWriterTest() {}
	
	@BeforeClass
	public static void setUpClass() throws Exception {}

	@AfterClass
	public static void tearDownClass() throws Exception {}

	@Before
	public void setUp() {}

	@After
	public void tearDown() {}
	
	@Test
	public void testProcess() throws Exception {
		File dnaFile = new File("src/test/resources/NM_000266.gb");
		LinkedHashMap<String, DNASequence> dnaSequences = GenbankReaderHelper.readGenbankDNASequence( dnaFile );
		ByteArrayOutputStream fragwriter = new ByteArrayOutputStream();
		ArrayList<DNASequence> seqs = new ArrayList<DNASequence>();
		for(DNASequence seq : dnaSequences.values()) {
			seqs.add(seq);
		}
		GenbankWriterHelper.writeNucleotideSequence(fragwriter, seqs,
				GenbankWriterHelper.LINEAR_DNA);
		//System.out.println(fragwriter.toString());
		ByteArrayInputStream fragreader = new ByteArrayInputStream(fragwriter.toByteArray());
		dnaSequences = GenbankReaderHelper.readGenbankDNASequence( fragreader );
		fragwriter.close();
		assertEquals(seqs.get(0).getSequenceAsString(),dnaSequences.values().iterator().next().getSequenceAsString());
	}
}
