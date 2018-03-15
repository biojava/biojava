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
/**
 *
 */
package org.biojava.nbio.structure.test.scop;

import org.biojava.nbio.structure.scop.ScopCategory;
import org.biojava.nbio.structure.scop.ScopDescription;
import org.junit.Assert;
import org.junit.Test;


/**
 * @author Spencer Bliven <sbliven@ucsd.edu>
 *
 */
public class ScopDescriptionTest {
	@Test
	public void testClassification() {
		ScopDescription s = new ScopDescription();
		s.setClassificationId("b.12.1.7");

		String scopClass;

		try {
			scopClass = s.getClassificationId(ScopCategory.Px);
			Assert.fail("Illegal category. Should have thrown IllegalArgumentException");
		} catch( IllegalArgumentException e) {
			//expected
		}

		scopClass = s.getClassificationId(ScopCategory.Class);
		Assert.assertEquals("b", scopClass);

		scopClass = s.getClassificationId(ScopCategory.Fold);
		Assert.assertEquals("b.12", scopClass);

		scopClass = s.getClassificationId(ScopCategory.Superfamily);
		Assert.assertEquals("b.12.1", scopClass);

		scopClass = s.getClassificationId(ScopCategory.Family);
		Assert.assertEquals("b.12.1.7", scopClass);


		s.setClassificationId("k.14");

		scopClass = s.getClassificationId(ScopCategory.Class);
		Assert.assertEquals("k", scopClass);

		scopClass = s.getClassificationId(ScopCategory.Fold);
		Assert.assertEquals("k.14", scopClass);

		scopClass = s.getClassificationId(ScopCategory.Superfamily);
		Assert.assertNull(scopClass);

		scopClass = s.getClassificationId(ScopCategory.Family);
		Assert.assertNull(scopClass);
	}
}
