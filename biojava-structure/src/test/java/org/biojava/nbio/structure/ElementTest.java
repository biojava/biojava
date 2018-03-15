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

package org.biojava.nbio.structure;

import junit.framework.TestCase;
import org.junit.Assert;
import org.junit.Test;

/**
 * Tests for Element class.
 *
 * @author Peter Rose
 * @since 3.0
 */
public class ElementTest {

	@Test
	public void testIsMetal() {
		Element h = Element.H;
		Assert.assertFalse(h.isMetal());
		Element he = Element.He;
		Assert.assertFalse(he.isMetal());
		Element li = Element.Li;
		Assert.assertTrue(li.isMetal());
		Element be = Element.Be;
		Assert.assertTrue(be.isMetal());
		Element b = Element.B;
		Assert.assertFalse(b.isMetal());
		Element c = Element.C;
		Assert.assertFalse(c.isMetal());
		Element f = Element.F;
		Assert.assertFalse(f.isMetal());
		Element al = Element.Al;
		Assert.assertTrue(al.isMetal());
		Element sc = Element.Sc;
		Assert.assertTrue(sc.isMetal());
		Element la = Element.La;
		Assert.assertTrue(la.isMetal());
		Element ac = Element.Ac;
		Assert.assertTrue(ac.isMetal());
	};

	@Test
	public void testIsMetalloid() {
		Element h = Element.H;
		Assert.assertFalse(h.isMetalloid());
		Element he = Element.He;
		Assert.assertFalse(he.isMetalloid());
		Element li = Element.Li;
		Assert.assertFalse(li.isMetalloid());
		Element be = Element.Be;
		Assert.assertFalse(be.isMetalloid());
		Element b = Element.B;
		Assert.assertTrue(b.isMetalloid());
		Element c = Element.C;
		Assert.assertFalse(c.isMetalloid());
		Element f = Element.F;
		Assert.assertFalse(f.isMetalloid());
		Element al = Element.Al;
		Assert.assertFalse(al.isMetalloid());
		Element sc = Element.Sc;
		Assert.assertFalse(sc.isMetalloid());
		Element la = Element.La;
		Assert.assertFalse(la.isMetalloid());
		Element ac = Element.Ac;
		Assert.assertFalse(ac.isMetalloid());
	}

	@Test
	public void testIsNonMetal() {
		Element h = Element.H;
		Assert.assertTrue(h.isNonMetal());
		Element he = Element.He;
		Assert.assertTrue(he.isNonMetal());
		Element li = Element.Li;
		Assert.assertFalse(li.isNonMetal());
		Element be = Element.Be;
		Assert.assertFalse(be.isNonMetal());
		Element b = Element.B;
		Assert.assertFalse(b.isNonMetal());
		Element c = Element.C;
		Assert.assertTrue(c.isNonMetal());
		Element f = Element.F;
		Assert.assertTrue(f.isNonMetal());
		Element al = Element.Al;
		Assert.assertFalse(al.isNonMetal());
		Element sc = Element.Sc;
		Assert.assertFalse(sc.isNonMetal());
		Element la = Element.La;
		Assert.assertFalse(la.isNonMetal());
		Element ac = Element.Ac;
		Assert.assertFalse(ac.isNonMetal());
	}

}
