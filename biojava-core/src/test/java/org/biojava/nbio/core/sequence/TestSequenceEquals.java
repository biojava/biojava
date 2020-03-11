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
package org.biojava.nbio.core.sequence;

import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

/**
 * Created by yana on 3/27/17.
 */
public class TestSequenceEquals {

	@Test
	public void testSameCompounds() throws Exception{



		ProteinSequence seq1 = new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX");

		ProteinSequence seq2 = new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJXARNDCEQGHILKMFPSTWYVBZJX");


        assertNotEquals(seq1, seq2);

        assertEquals(seq1, seq1);

        assertEquals(seq2, seq2);


		ProteinSequence seq3 = new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX");

        assertEquals(seq3, seq1);

        assertNotEquals(seq2, seq3);


		DNASequence dnaSeq = new DNASequence("ATGGCGGCGCTGAGCGGT");

        assertNotEquals(seq1, dnaSeq);



	}
}
