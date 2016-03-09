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
package org.biojava.nbio.core.sequence.views;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.SequenceView;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.*;

public class WindowViewTests {

	@Test
	public void basicWindow() throws CompoundNotFoundException {
		RNASequence rna = new RNASequence("AUGCCU");
		WindowedSequence<NucleotideCompound> window = new WindowedSequence<NucleotideCompound>(rna, 3);

		Iterator<SequenceView<NucleotideCompound>> iter = window.iterator();
		assertTrue("hasNext() returns true", iter.hasNext());

		int count = 0;
		for(SequenceView<NucleotideCompound> c: window) {
			count++;
			if(count == 0) {
				String extracted = c.getCompoundAt(0).toString() + c.getCompoundAt(1).toString() + c.getCompoundAt(2).toString();
				assertEquals("Checking codon string", "AUG", extracted);
			}
		}
		assertThat("Windowed iterator broken", count, is(2));
	}

	@Test
	public void reaminderWindow() throws CompoundNotFoundException {
		RNASequence rna = new RNASequence("AUGCC");
		WindowedSequence<NucleotideCompound> window = new WindowedSequence<NucleotideCompound>(rna, 3);
		List<SequenceView<NucleotideCompound>> list = new ArrayList<SequenceView<NucleotideCompound>>();
		for(SequenceView<NucleotideCompound> c: window) {
			list.add(c);
		}
		assertThat("First window is size 3", list.get(0).getLength(), is(3));
		assertThat("Only 1 window", list.size(), is(1));
	}
}
