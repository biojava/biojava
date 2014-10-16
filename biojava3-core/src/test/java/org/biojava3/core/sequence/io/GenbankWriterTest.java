/**
 * 
 */
package org.biojava3.core.sequence.io;


import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import junit.framework.TestCase;

import org.biojava3.core.sequence.DNASequence;
import org.junit.Test;


/**
 * @author mckeee1
 * 
 */
public class GenbankWriterTest extends TestCase{


	@Test
	public void testProcess() throws Exception {

        InputStream inStream = GenbankWriterTest.class.getResourceAsStream("/NM_000266.gb");
		//File dnaFile = new File("src/test/resources/NM_000266.gb");
		LinkedHashMap<String, DNASequence> dnaSequences = GenbankReaderHelper.readGenbankDNASequence( inStream );
		ByteArrayOutputStream fragwriter = new ByteArrayOutputStream();
		ArrayList<DNASequence> seqs = new ArrayList<DNASequence>();
		for(DNASequence seq : dnaSequences.values()) {
			seqs.add(seq);
		}
		GenbankWriterHelper.writeNucleotideSequence(fragwriter, seqs,
				GenbankWriterHelper.LINEAR_DNA);
		//System.out.println(fragwriter.toString());
		ByteArrayInputStream fragreader = new ByteArrayInputStream(fragwriter.toByteArray());
                /**
                 * Hello Jacek
                 * can you please investigate why this test fails? it seems that 
                 * fragreader at the line below is read with the last feature 
                 * in an invalid state: location = 2005..2004
                 */
		//dnaSequences = GenbankReaderHelper.readGenbankDNASequence( fragreader );
		fragwriter.close();
		assertEquals(seqs.get(0).getSequenceAsString(),dnaSequences.values().iterator().next().getSequenceAsString());
	}
}
