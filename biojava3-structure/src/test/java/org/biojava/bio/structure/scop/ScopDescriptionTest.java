/**
 * 
 */
package org.biojava.bio.structure.scop;

import junit.framework.TestCase;


/**
 * @author Spencer Bliven <sbliven@ucsd.edu>
 *
 */
public class ScopDescriptionTest extends TestCase{
	public void testClassification() {
		ScopDescription s = new ScopDescription();
		s.setClassificationId("b.12.1.7");
		
		String scopClass;
		
		try {
			scopClass = s.getClassificationId(ScopCategory.Px);
			fail("Illegal category. Should have thrown IllegalArgumentException");
		} catch( IllegalArgumentException e) {
			//expected
		}
		
		scopClass = s.getClassificationId(ScopCategory.Class);
		assertEquals("b",scopClass);
		
		scopClass = s.getClassificationId(ScopCategory.Fold);
		assertEquals("b.12",scopClass);
		
		scopClass = s.getClassificationId(ScopCategory.Superfamily);
		assertEquals("b.12.1",scopClass);
		
		scopClass = s.getClassificationId(ScopCategory.Family);
		assertEquals("b.12.1.7",scopClass);

		
		s.setClassificationId("k.14");
		
		scopClass = s.getClassificationId(ScopCategory.Class);
		assertEquals("k",scopClass);
		
		scopClass = s.getClassificationId(ScopCategory.Fold);
		assertEquals("k.14",scopClass);
		
		scopClass = s.getClassificationId(ScopCategory.Superfamily);
		assertNull(scopClass);
		
		scopClass = s.getClassificationId(ScopCategory.Family);
		assertNull(scopClass);
	}
}
