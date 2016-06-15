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

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.storage.JoiningSequenceReader;
import org.biojava.nbio.core.sequence.template.SequenceMixin;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class JoiningSequenceReaderTest {

	@SuppressWarnings("unchecked")
	@Test
	public void canScan() throws CompoundNotFoundException {
		JoiningSequenceReader<NucleotideCompound> seq =
			new JoiningSequenceReader<NucleotideCompound>(
					new DNASequence("AAAA"),
					new DNASequence("GGG"),
					new JoiningSequenceReader<NucleotideCompound>(new DNASequence("A"), new DNASequence("C")),
					new DNASequence("TT"),
					new DNASequence("C")
		);

		String expected = "AAAAGGGACTTC";

		StringBuilder builderByIndex = new StringBuilder();
		for(int i = 1; i <= seq.getLength(); i++) {
			builderByIndex.append(seq.getCompoundAt(i));
		}

		StringBuilder builderByIterator = SequenceMixin.toStringBuilder(seq);

		assertEquals("Index builder", expected, builderByIndex.toString());
		assertEquals("Iterator builder", expected, builderByIterator.toString());
	}

	@SuppressWarnings("unchecked")
	@Test
	public void empty() throws CompoundNotFoundException {
		JoiningSequenceReader<NucleotideCompound> seq =
			new JoiningSequenceReader<NucleotideCompound>(
					new DNASequence(""),
					new DNASequence(""),
					new DNASequence("A"),
					new DNASequence("")
			);
		assertEquals("Testing empty sequences", "A", seq.getSequenceAsString());
	}
}
