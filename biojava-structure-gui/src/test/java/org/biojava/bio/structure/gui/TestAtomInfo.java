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
package org.biojava.bio.structure.gui;

import junit.framework.TestCase;

import org.biojava.bio.structure.align.gui.jmol.AtomInfo;


public class TestAtomInfo extends TestCase {

	public void testAtomInfoConversion(){
		String s1 = "[MET]508:A.CA/1 #3918";

		AtomInfo aa = AtomInfo.fromString(s1);
		assertTrue(aa.getAtomName().equals("CA"));
		assertEquals(aa.getChainId(),"A");
		assertEquals(aa.getModelNumber(),1);
		assertEquals(aa.getResidueName(),"MET");
		assertEquals(aa.getResidueNumber(), "508");

		String s1New = aa.toString();
		// atom nr not supported yet
		assertEquals("[MET]508:A.CA/1",s1New);

	}

	public void testInsertionCode(){
		String s1 = "[ASP]1^A:A.CA/2 #2";

		AtomInfo aa = AtomInfo.fromString(s1);

		assertEquals(aa.getAtomName(), "CA");
		assertEquals(aa.getChainId(),"A");
		assertEquals(aa.getModelNumber(),2);
		assertEquals(aa.getResidueName(),"ASP");		
		assertEquals(aa.getResidueNumber(), "1A");
		assertEquals("[ASP]1^A:A.CA/2",aa.toString());

	}
}
