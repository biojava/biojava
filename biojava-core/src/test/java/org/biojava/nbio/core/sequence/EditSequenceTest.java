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
import org.biojava.nbio.core.sequence.edits.Edit;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class EditSequenceTest {

	@Test
	public void substitute() throws CompoundNotFoundException {
		DNASequence seq = new DNASequence("ACGT");
		assertSeq(new Edit.Substitute<NucleotideCompound>("T", 2).edit(seq), "ATGT");
		assertSeq(new Edit.Substitute<NucleotideCompound>("TT", 2).edit(seq), "ATTT");
		assertSeq(new Edit.Substitute<NucleotideCompound>("T", 1).edit(seq), "TCGT");
		assertSeq(new Edit.Substitute<NucleotideCompound>("TTC", 2).edit(seq), "ATTC");
	}

	@Test(expected=IndexOutOfBoundsException.class)
	public void badSubstitute() throws CompoundNotFoundException {
		new Edit.Substitute<NucleotideCompound>("AAAA", 4).edit(new DNASequence("ACGT"));
	}

	@Test
	public void delete() throws CompoundNotFoundException {
		DNASequence seq = new DNASequence("ACGT");
		assertSeq(new Edit.Delete<NucleotideCompound>(1).edit(seq), "CGT");
		assertSeq(new Edit.Delete<NucleotideCompound>(4).edit(seq), "ACG");
		assertSeq(new Edit.Delete<NucleotideCompound>(2,3).edit(seq), "AT");

		//disabling this test, because we can't create a CompoundSet if we have no sequences...
		// assertSeq(new Edit.Delete<NucleotideCompound>(1,4).edit(seq), "");
	}

	@Test
	public void insert() throws CompoundNotFoundException {
		DNASequence seq = new DNASequence("ACGT");
		assertSeq(new Edit.Insert<NucleotideCompound>("TT", 1).edit(seq), "TTACGT");
		assertSeq(new Edit.Insert<NucleotideCompound>("TT", 2,3).edit(seq), "ACTTGT");
		assertSeq(new Edit.Insert<NucleotideCompound>("TT", 3,4).edit(seq), "ACGTTT");
		assertSeq(new Edit.Insert<NucleotideCompound>("A", 4).edit(seq), "ACGTA");

		//Original BioJava example
		assertSeq(
				new Edit.Insert<NucleotideCompound>("atgga", 3,4).edit(new DNASequence("gataca")),
				"gatatggaaca"
		);
	}

	private void assertSeq(Sequence<? extends Compound> seq, String expected) {
		assertEquals("Asserting sequence "+expected, expected, seq.getSequenceAsString());
	}
}
