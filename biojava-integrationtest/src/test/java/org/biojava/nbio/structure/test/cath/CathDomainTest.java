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
package org.biojava.nbio.structure.test.cath;

import org.biojava.nbio.structure.cath.CathDomain;
import org.biojava.nbio.structure.cath.CathFactory;
import org.junit.Test;

import static org.junit.Assert.assertEquals;


/**
 * A test for {@link CathDomain}.
 * @author dmyersturnbull
 */
public class CathDomainTest {
	@Test
	public void test() {
		String id = "1qvrC03";
		CathDomain domain = CathFactory.getCathDatabase().getDomainByCathId(id);
		assertEquals("1qvr.C_332-400,C_514-540", domain.toCanonical().getIdentifier());
	}
}
