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
 * Created on Oct 5, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava3.alignment;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

import org.biojava3.alignment.template.Profile;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.util.ConcurrencyTools;

import junit.framework.TestCase;

public class TestDNAAlignment extends TestCase{


	public void testDNAAlignment(){

		try {
			List<DNASequence> lst = getDNAFASTAFile();

			Profile<DNASequence, NucleotideCompound> profile = Alignments.getMultipleSequenceAlignment(lst);

			assertTrue(profile.getSize() == 10);

			assertTrue(profile.getAlignedSequence(1).getSequenceAsString().length() > 50);
			
			
			// here how to print the MSA:
			
			//System.out.printf("MSA:%n%s%n", profile);
		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}
		ConcurrencyTools.shutdown();
	}

	private static List<DNASequence> getDNAFASTAFile() throws Exception {

		InputStream inStream = TestDNAAlignment.class.getResourceAsStream(String.format("/dna-fasta.txt")); 
		LinkedHashMap<String, DNASequence> fastas = FastaReaderHelper.readFastaDNASequence(inStream);

		List<DNASequence> sequences = new ArrayList<DNASequence>();

		for ( String key: fastas.keySet()){
			DNASequence seq = fastas.get(key);
			sequences.add(seq);
		}

		return sequences;
	}


}
