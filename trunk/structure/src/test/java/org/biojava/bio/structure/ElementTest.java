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

package org.biojava.bio.structure;

import junit.framework.TestCase;

/**
 * Tests for Element class.
 *
 * @author Peter Rose
 * @since 3.0
 * @version %I% %G%
 */
public class ElementTest extends TestCase {

	public void testIsMetal() {
		Element h = Element.H;
		assertFalse(h.isMetal());
		Element he = Element.He;
		assertFalse(he.isMetal());
		Element li = Element.Li;
		assertTrue(li.isMetal());
		Element be = Element.Be;
		assertTrue(be.isMetal());
		Element b = Element.B;
		assertFalse(b.isMetal());
		Element c = Element.C;
		assertFalse(c.isMetal());
		Element f = Element.F;
		assertFalse(f.isMetal());
		Element al = Element.Al;
		assertTrue(al.isMetal());
		Element sc = Element.Sc;
		assertTrue(sc.isMetal());
		Element la = Element.La;
		assertTrue(la.isMetal());
		Element ac = Element.Ac;
		assertTrue(ac.isMetal());
	};

	public void testIsMetalloid() {
		Element h = Element.H;
		assertFalse(h.isMetalloid());
		Element he = Element.He;
		assertFalse(he.isMetalloid());
		Element li = Element.Li;
		assertFalse(li.isMetalloid());
		Element be = Element.Be;
		assertFalse(be.isMetalloid());
		Element b = Element.B;
		assertTrue(b.isMetalloid());
		Element c = Element.C;
		assertFalse(c.isMetalloid());
		Element f = Element.F;
		assertFalse(f.isMetalloid());
		Element al = Element.Al;
		assertFalse(al.isMetalloid());
		Element sc = Element.Sc;
		assertFalse(sc.isMetalloid());
		Element la = Element.La;
		assertFalse(la.isMetalloid());
		Element ac = Element.Ac;
		assertFalse(ac.isMetalloid());
	}

	public void testIsNonMetal() {
		Element h = Element.H;
		assertTrue(h.isNonMetal());
		Element he = Element.He;
		assertTrue(he.isNonMetal());
		Element li = Element.Li;
		assertFalse(li.isNonMetal());
		Element be = Element.Be;
		assertFalse(be.isNonMetal());
		Element b = Element.B;
		assertFalse(b.isNonMetal());
		Element c = Element.C;
		assertTrue(c.isNonMetal());
		Element f = Element.F;
		assertTrue(f.isNonMetal());
		Element al = Element.Al;
		assertFalse(al.isNonMetal());
		Element sc = Element.Sc;
		assertFalse(sc.isNonMetal());
		Element la = Element.La;
		assertFalse(la.isNonMetal());
		Element ac = Element.Ac;
		assertFalse(ac.isNonMetal());
	}

}
