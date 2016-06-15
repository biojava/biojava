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
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.junit.Test;

import java.util.Collection;

import static org.junit.Assert.*;

public class CasePreservingProteinSequenceCreatorTest {

	@Test
	public void testConstructor() throws CompoundNotFoundException {
		CasePreservingProteinSequenceCreator creator = new CasePreservingProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet());

		String seq = "aCDEfgHI-Jkl";
		ProteinSequence prot = (ProteinSequence) creator.getSequence(seq, 0);
		Collection<Object> uppercase = prot.getUserCollection();

		//test some assumptions. Hopefully work on non-english locals too?
		assertFalse(Character.isUpperCase('-'));
		assertFalse(Character.isUpperCase('.'));

		assertEquals("Lengths differ",seq.length(),uppercase.size());

		int i=0;
		for(Object obj : uppercase) {
			assertTrue("Not a Boolean",obj instanceof Boolean);
			Boolean bool = (Boolean)obj;
			assertEquals("Doesn't match case of "+seq.charAt(i),Character.isUpperCase(seq.charAt(i)),bool);
			i++;
		}
	}
}
