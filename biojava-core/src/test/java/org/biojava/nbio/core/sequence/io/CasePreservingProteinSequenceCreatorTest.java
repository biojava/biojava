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

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Collection;
import java.util.Iterator;

class CasePreservingProteinSequenceCreatorTest {

	@Test
	void testConstructor() throws CompoundNotFoundException {
		CasePreservingProteinSequenceCreator creator = new CasePreservingProteinSequenceCreator(
				AminoAcidCompoundSet.getAminoAcidCompoundSet());

		String seq = "aCDEfgHI-Jkl";
		ProteinSequence prot = (ProteinSequence) creator.getSequence(seq, 0);
		Collection<Object> uppercase = prot.getUserCollection();

		// test some assumptions. Hopefully work on non-english locals too?
		assertFalse(Character.isUpperCase('-'));
		assertFalse(Character.isUpperCase('.'));

		assertEquals(seq.length(), uppercase.size(), "Lengths differ");

		int i = 0;
		for (Object obj : uppercase) {
			assertTrue(obj instanceof Boolean, "Not a Boolean");
			Boolean bool = (Boolean) obj;
			assertEquals(Character.isUpperCase(seq.charAt(i)), bool, "Doesn't match case of " + seq.charAt(i));
			i++;
		}
	}

	@Test
	void booleanConversion() throws CompoundNotFoundException {
		CasePreservingProteinSequenceCreator creator = new CasePreservingProteinSequenceCreator(
				AminoAcidCompoundSet.getAminoAcidCompoundSet());
		AbstractSequence<AminoAcidCompound> seq = creator.getSequence("aaAA", 0);
		assertEquals("AAAA", seq.getSequenceAsString());
		Boolean[] expected = new Boolean[] { Boolean.FALSE, Boolean.FALSE, Boolean.TRUE, Boolean.TRUE };
		Iterator<Object> userCollection = seq.getUserCollection().iterator();
		for (int i = 0; i < seq.getLength(); i++) {
			assertEquals(expected[i], userCollection.next());
		}
	}

}
