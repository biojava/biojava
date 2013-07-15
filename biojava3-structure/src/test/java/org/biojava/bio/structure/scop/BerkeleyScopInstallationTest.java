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
 */
package org.biojava.bio.structure.scop;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;

/**
 * Tests {@link BerkeleyScopInstallation}.
 * @author dmyerstu
 * @since 3.0.6
 */
public class BerkeleyScopInstallationTest {

	@Test
	public void testComments() {
		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75B);
		List<String> comments;
		
		// d3ueea_ was added in 1.75B update
		// domain
		comments = scop.getComments(190700);
		assertEquals("Wrong number of comments", 1, comments.size());
		assertEquals("Wrong comment", "not a true protein", comments.get(0));

		// fold
		comments = scop.getComments(57923);
		assertEquals("Wrong number of comments", 1, comments.size());
		assertEquals("Wrong comment", "metal(zinc)-bound alpha+beta fold", comments.get(0));
	}
	
}
