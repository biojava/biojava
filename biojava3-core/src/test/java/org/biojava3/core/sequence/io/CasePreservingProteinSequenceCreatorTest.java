package org.biojava3.core.sequence.io;

import java.util.Collection;

import junit.framework.TestCase;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.CasePreservingProteinSequenceCreator;

public class CasePreservingProteinSequenceCreatorTest extends TestCase {

	public void testConstructor() {
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
			assertEquals("Doesn't match case of "+seq.charAt(i),(Boolean)Character.isUpperCase(seq.charAt(i)),bool);
			i++;
		}
	}
}
