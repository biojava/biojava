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
package org.biojava.nbio.structure.test.util;

import org.biojava.nbio.core.util.StringManipulationHelper;
import org.junit.Assert;
import org.junit.Test;

public class StringManipulationTestsHelper {

	/**
	 * Asserts that two strings are equal, line by line, ignoring any difference
	 * of end line delimiters contained within the 2 Strings. This method should
	 * be used if and only if two XMLs are considered identical when all nodes
	 * are identical including their relative order. Generally useful when
	 * asserting identity of <b>automatically regenerated</b> XML.
	 *
	 * @param expected
	 * @param actual
	 */
	public static void assertEqualsIgnoreEndline(String expected, String actual) {
		Assert.assertTrue(StringManipulationHelper.equalsToIgnoreEndline(expected, actual));
	}

	/**
	 * Asserts that two XML-representing strings are equal, by recursively
	 * comparing each node's set of properties & children nodes. This method
	 * should be used when two XMLs are considered identical when all nodes are
	 * identical regardless to their order
	 *
	 * @param expectedXml
	 * @param actualXml
	 */
	public static void assertEqualsXml(String expectedXml, String actualXml) {
		Assert.assertTrue(StringManipulationHelper.equalsToXml(expectedXml, actualXml));
	}

	/**
	 * @param s
	 * @param pdb
	 */
	public static void compareString(String s, String pdb) {
		for (int i = 0; i < s.length(); i++) {
			System.out.println("@" + i + "\t>" + s.charAt(i) + ":"
					+ pdb.charAt(i) + "<\t" + Integer.toHexString(s.charAt(i))
					+ ":" + Integer.toHexString(pdb.charAt(i)));
			if (Character.toUpperCase(s.charAt(i)) != Character.toUpperCase(pdb.charAt(i))) {
				break;
			}
		}
	}

	@Test
	public void testNothing(){
		Assert.assertTrue(true);
	}
}
