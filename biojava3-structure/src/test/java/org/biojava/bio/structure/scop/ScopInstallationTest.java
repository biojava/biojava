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
 * Created on Jul 1, 2013
 * Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.bio.structure.scop;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;



public class ScopInstallationTest {

	@Test
	public void testComments() {
		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75);
		List<String> comments;
		
		comments = scop.getComments(127355);
		assertEquals("Wrong number of comments", 2, comments.size());
		assertEquals("Wrong comment", "automatically matched to d2hbia_", comments.get(0));
		assertEquals("Wrong comment", "complexed with hem; mutant", comments.get(1));
		
		comments = scop.getComments(160555);
		assertEquals("Wrong number of comments", 1, comments.size());
		assertEquals("Wrong comment", "<a href=\"http://pfam.sanger.ac.uk/family?acc=PF06262\">PF06262</a>; DUF1025; minimal zincin fold that retains 3-stranded mixed beta-sheet and contains HExxH motif in the C-terminal helix; no metal ion bound to this motif is observed in the first determined structures", comments.get(0));
	}
	
}
