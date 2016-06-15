/**
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
 * Created on Dec 7, 2013
 * Created by Douglas Myers-Turnbull
 *
 */
package org.biojava.nbio.structure.io.sifts;

import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertEquals;

/**
 * Tests {@link SiftsChainToUniprotMapping}.
 * @author dmyersturnbull
 * @since 3.0.7
 */
public class SiftsChainToUniprotMappingTest {

	@Test
	public void test() throws IOException {
		SiftsChainToUniprotMapping sifts = SiftsChainToUniprotMapping.load();
		SiftsChainEntry entry = sifts.getByChainId("1hiv", "A");
		assertEquals("1hiv", entry.getPdbId());
		assertEquals("A", entry.getChainId());
		assertEquals("P04585", entry.getUniProtId());
	}

}
