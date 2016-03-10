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
package org.biojava.nbio.core.sequence.io;


import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.DNASequence;
import org.junit.Test;

import java.io.ByteArrayOutputStream;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

public class FastaWriterTest {

	@Test
	public void writeBasicFasta() throws Exception {
		String id         = "Example";
		String dnaLineOne = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		String dnaLineTwo = "T";

		DNASequence s = new DNASequence(dnaLineOne+dnaLineTwo);
		s.setAccession(new AccessionID(id));
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		FastaWriterHelper.writeSequence(baos, s);

		String actual = new String(baos.toByteArray());
		String expected = String.format(">%s%n%s%n%s%n", id, dnaLineOne, dnaLineTwo);

		assertThat("Writer not as expected", actual, is(expected));
	}

	@Test
	public void writeFastaEqualToLineLength() throws Exception {
		String id  = "Example";
		String dna = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT";

		DNASequence s = new DNASequence(dna);
		s.setAccession(new AccessionID(id));
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		FastaWriterHelper.writeSequence(baos, s);

		String actual = new String(baos.toByteArray());
		String expected = String.format(">%s%n%s%n", id, dna);

		assertThat("Writer not as expected", actual, is(expected));
	}

}
