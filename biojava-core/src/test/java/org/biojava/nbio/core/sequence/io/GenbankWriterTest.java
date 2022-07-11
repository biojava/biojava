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
/**
 *
 */
package org.biojava.nbio.core.sequence.io;


import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.features.AbstractFeature;
import org.biojava.nbio.core.sequence.features.Qualifier;
import org.biojava.nbio.core.sequence.features.TextFeature;
import org.biojava.nbio.core.sequence.location.SimpleLocation;
import org.biojava.nbio.core.sequence.Strand;
import org.junit.Assert;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;


/**
 * @author mckeee1
 *
 */
public class GenbankWriterTest {


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
		Assert.assertEquals(seqs.get(0).getSequenceAsString(), dnaSequences.values().iterator().next().getSequenceAsString());
	}
	
	/**
	 * String Formatter error when key or value of Qualifier has character "%"
	 * https://github.com/biojava/biojava/issues/886
	 */
	@Test
	public void testGithub886() throws Exception {
		
		DNASequence seq = new DNASequence("ATGC");
		seq.setAccession(new AccessionID("."));
		AbstractFeature feature = new TextFeature("CDS", "source", "short description", "description");
		feature.setLocation(new SimpleLocation(1, 10, Strand.POSITIVE));

		// no percent symbols in key or value
		feature.addQualifier("note1", new Qualifier("note1", "50", true));
		// percent symbol in key
		feature.addQualifier("note2", new Qualifier("%note2", "50", true));
		feature.addQualifier("note3", new Qualifier("not%e3", "50", true));
		feature.addQualifier("note4", new Qualifier("note4%", "50", true));
		// percent symbol in value
		feature.addQualifier("note5", new Qualifier("note5", "%50", true));
		feature.addQualifier("note6", new Qualifier("note6", "5%0", true));
		feature.addQualifier("note7", new Qualifier("note7", "50%", true));
		
		seq.addFeature(feature);
		
		ByteArrayOutputStream fragwriter = new ByteArrayOutputStream();
		GenbankWriterHelper.writeNucleotideSequence(
				fragwriter, 
				Arrays.asList(seq), 
				GenbankWriterHelper.LINEAR_DNA);
		fragwriter.close();
		System.out.println(fragwriter.toString().replaceAll("\r\n", "\n"));
		
		// now read in the file that was created and check that the qualifiers were created correctly
		InputStream readerInputStream = new ByteArrayInputStream(fragwriter.toByteArray());
		DNASequence newSeq = GenbankReaderHelper.readGenbankDNASequence(readerInputStream).values().iterator().next();
		AbstractFeature newFeature = (TextFeature) seq.getFeaturesByType("CDS").get(0);
		Map<String, List<Qualifier>> newQualifiers = newFeature.getQualifiers();
		
		assertEquals("note1", newQualifiers.get("note1").get(0).getName());
		assertEquals("50", newQualifiers.get("note1").get(0).getValue());
		
		assertEquals("%note2", newQualifiers.get("note2").get(0).getName());
		assertEquals("50", newQualifiers.get("note2").get(0).getValue());
		
		assertEquals("not%e3", newQualifiers.get("note3").get(0).getName());
		assertEquals("50", newQualifiers.get("note3").get(0).getValue());
		
		assertEquals("note4%", newQualifiers.get("note4").get(0).getName());
		assertEquals("50", newQualifiers.get("note4").get(0).getValue());
		
		assertEquals("note5", newQualifiers.get("note5").get(0).getName());
		assertEquals("%50", newQualifiers.get("note5").get(0).getValue());
		
		assertEquals("note6", newQualifiers.get("note6").get(0).getName());
		assertEquals("5%0", newQualifiers.get("note6").get(0).getValue());
		
		assertEquals("note7", newQualifiers.get("note7").get(0).getName());
		assertEquals("50%", newQualifiers.get("note7").get(0).getValue());
		
	}
}
