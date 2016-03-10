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
package org.biojava.nbio.core.sequence.compound;

import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class AmbiguityDNACompoundTest {

	private AmbiguityDNACompoundSet set = AmbiguityDNACompoundSet.getDNACompoundSet();

	@Test
	public void testAmbiguity() {
		NucleotideCompound actual = set.getAmbiguity(getCompounds("M","V"));
		assertEquals("Checking M & G = V", getCompounds("V")[0], actual);
	}

	@Test
	public void testBasicAmbiguity() {
		NucleotideCompound actual = set.getAmbiguity(getCompounds("A","C"));
		assertEquals("Checking A & C = M", getCompounds("M")[0], actual);
	}

	private NucleotideCompound[] getCompounds(String... compoundStrings) {
		List<NucleotideCompound> c = new ArrayList<NucleotideCompound>();
		for(String s: compoundStrings) {
			c.add(set.getCompoundForString(s));
		}
		return c.toArray(new NucleotideCompound[0]);
	}
}
