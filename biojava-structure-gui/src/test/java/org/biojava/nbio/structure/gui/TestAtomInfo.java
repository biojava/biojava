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
 * Created on Jan 21, 2010
 *
 */
package org.biojava.nbio.structure.gui;

import org.biojava.nbio.structure.align.gui.jmol.AtomInfo;
import org.junit.Assert;
import org.junit.Test;


public class TestAtomInfo {

	@Test
	public void testAtomInfoConversion(){
		String s1 = "[MET]508:A.CA/1 #3918";

		AtomInfo aa = AtomInfo.fromString(s1);
		Assert.assertTrue(aa.getAtomName().equals("CA"));
		Assert.assertEquals(aa.getChainId(), "A");
		Assert.assertEquals(aa.getModelNumber(), 1);
		Assert.assertEquals(aa.getResidueName(), "MET");
		Assert.assertEquals(aa.getResidueNumber(), "508");

		String s1New = aa.toString();
		// atom nr not supported yet
		Assert.assertEquals("[MET]508:A.CA/1", s1New);

	}

	@Test
	public void testInsertionCode(){
		String s1 = "[ASP]1^A:A.CA/2 #2";

		AtomInfo aa = AtomInfo.fromString(s1);

		Assert.assertEquals(aa.getAtomName(), "CA");
		Assert.assertEquals(aa.getChainId(), "A");
		Assert.assertEquals(aa.getModelNumber(), 2);
		Assert.assertEquals(aa.getResidueName(), "ASP");
		Assert.assertEquals(aa.getResidueNumber(), "1A");
		Assert.assertEquals("[ASP]1^A:A.CA/2", aa.toString());

	}
}
