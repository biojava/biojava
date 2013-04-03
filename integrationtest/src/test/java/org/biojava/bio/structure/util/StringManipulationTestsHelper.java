package org.biojava.bio.structure.util;

import junit.framework.TestCase;

import org.biojava3.core.util.StringManipulationHelper;

public class StringManipulationTestsHelper extends TestCase {

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
		assertTrue(StringManipulationHelper.equalsToIgnoreEndline(expected, actual));
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
		assertTrue(StringManipulationHelper.equalsToXml(expectedXml, actualXml));
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
	
	public void testNothing(){
		assertTrue(true);
	}
}
