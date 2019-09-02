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
 * Created on May 11, 2010
 * Author: Andreas Prlic
 *
 */

package org.biojava.nbio.structure;

import java.io.IOException;

import static org.junit.Assert.*;

import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

/** make sure DNA alignments fail gracefully
 *
 * @author Andreas Prlic
 *
 */
public class TestDNAAlignment
{

	@Test
	public void test1() throws IOException {
		String name1="1l3s.A";
		String name2="1t7p.P";

		AtomCache cache = new AtomCache();
		try {
			Atom[] ca1 = cache.getAtoms(name1);
			Atom[] ca2 = cache.getAtoms(name2);
			CeMain ce = new CeMain();
			AFPChain afpChain = ce.align(ca1,ca2);
			assertNotNull(afpChain);

		  String txt = afpChain.toFatcat(ca1, ca2);

		  assertNotNull(txt);

		} catch (StructureException e){
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
}
