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
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.SequenceView;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class SequenceViewTest {

	@Test
	public void testGetCompoundAt() throws CompoundNotFoundException {
		SequenceView<NucleotideCompound> s = new DNASequence("ATGC").getSubSequence(1, 4);
		assertEquals("Compound @ 1", s.getCompoundAt(1).toString(), "A");
		assertEquals("Compound @ 3", s.getCompoundAt(3).toString(), "G");
		assertEquals("Compound @ 3", s.getSubSequence(2,3).getCompoundAt(1).toString(), "T");
	}

	@Test
	public void testLastIndexOf() throws CompoundNotFoundException {
		SequenceView<NucleotideCompound> s = new DNASequence("ATGC").getSubSequence(1, 4);
		CompoundSet<NucleotideCompound> cs = s.getCompoundSet();
		assertEquals("Last index of ", 4, s.getLastIndexOf(cs.getCompoundForString("C")));

		s = new DNASequence("GAAAAAAAAG").getSubSequence(1, 10);
		assertEquals("Last index of G is 10", 10,
				s.getLastIndexOf(cs.getCompoundForString("G")));
		assertEquals("Last index of G is 5", 6,
				s.getSubSequence(5, 10).getLastIndexOf(cs.getCompoundForString("G")));
	}

	@Test
	public void testInverse() throws CompoundNotFoundException {
		SequenceView<NucleotideCompound> s = new DNASequence("ATGC").getSubSequence(2, 3).getInverse();
		assertEquals("Reversed complementing view", s.getSequenceAsString(), "CA");
	}
}
