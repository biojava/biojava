package org.biojava.bio.structure.util;

import java.util.Scanner;

import junit.framework.TestCase;

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
		Scanner scanner1 = new Scanner(expected);
		Scanner scanner2 = new Scanner(actual);
		String line1, line2;
		while (scanner1.hasNextLine()) {
			line1 = scanner1.nextLine();
			line2 = scanner2.nextLine();
			assertEquals(line1, line2);
		}
		if (scanner2.hasNextLine()) {
			// force fail
			assertEquals(expected, actual);
		}
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
		throw new UnsupportedOperationException("not yet implemented");
	}
	
	public static void compareString(String t, String pdb) {
		for (int i = 0; i < t.length(); i++) {
			System.out.println("@" + i + "\t>" + t.charAt(i) + ":"
					+ pdb.charAt(i) + "<\t" + Integer.toHexString(t.charAt(i))
					+ ":" + Integer.toHexString(pdb.charAt(i)));
			if (Character.toUpperCase(t.charAt(i)) != Character.toUpperCase(pdb
					.charAt(i))) {
				break;
			}
		}
	}
}
