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
package org.biojava.nbio.alignment;

import static org.junit.Assert.*;

import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.util.ConcurrencyTools;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

public class TestSubOptimalMSA {

	private List<DNASequence> sequences = new ArrayList<DNASequence>();

	@Before
	public void setUp() {
		try {
			sequences.add(new DNASequence("TTGGGGCCTCTAAACGGGGTCTT"));
			sequences.add(new DNASequence("TTGGGGCCTCTAAACGGGTCTT"));
			sequences.add(new DNASequence("TTGGGGCTCTAACGGGTCTT"));
		} catch (CompoundNotFoundException e) {
			e.printStackTrace();
		}
	}

	@Test
	public void gapPenalty52() {
		SimpleGapPenalty gapP = new SimpleGapPenalty((short) 5, (short) 2);
		Profile<DNASequence, NucleotideCompound> msa = Alignments
				.getMultipleSequenceAlignment(sequences, gapP);

		assertEquals("TTGGGGCCTCTAAACGGGGTCTT" + System.lineSeparator()
				+ "TTGGGGCCTCTAAACGGG-TCTT"    + System.lineSeparator()
				+ "TTGGGGC-TCTAA-CGGG-TCTT"    + System.lineSeparator(),
				msa.toString());

		ConcurrencyTools.shutdown();
	}

	@Test @Ignore
	public void gapPenaltyDefault() {
		// Default is currently 10-1
		SimpleGapPenalty gapP = new SimpleGapPenalty((short) 10, (short) 1);
		Profile<DNASequence, NucleotideCompound> msa = Alignments
				.getMultipleSequenceAlignment(sequences, gapP);

		// TODO test not passing (see issue 288 in github) - Aleix 03.2016
		assertEquals("TTGGGGCCTCTAAACGGGGTCTT" + System.lineSeparator()
				+ "TTGGGGCCTCTAAACGGG-TCTT"    + System.lineSeparator()
				+ "TTGGGGC-TCTAA-CGGG-TCTT"    + System.lineSeparator(),
				msa.toString());

		ConcurrencyTools.shutdown();
	}
}
