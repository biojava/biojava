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
 * Created on Sep 14, 2011
 * Author: Amr AL-Hossary 
 *
 */
package org.biojava3.core.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Scanner;

import junit.framework.TestCase;

/**
 * A utility class for common {@link String} manipulation tasks.
 * All functions are static methods.
 * 
 * @author Amr AL-Hossary
 */
public abstract class StringManipulationTestsHelper extends TestCase {

	/**
	 * we are using Unix endline here, since this is used for testing XML and it
	 * is part of the XML recommendations: <a href
	 * ="http://www.w3.org/TR/REC-xml/#sec-line-ends"
	 * >http://www.w3.org/TR/REC-xml/#sec-line-ends</a>
	 */
	private static final String UNIX_NEWLINE = "\n";

	private StringManipulationTestsHelper() {
		// to prevent instantiation
	}

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

	public static String convertStreamToString(InputStream stream) {
		BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
		StringBuilder sb = new StringBuilder();

		String line = null;
		try {
			while ((line = reader.readLine()) != null) {

				sb.append(line + UNIX_NEWLINE);
			}
		} catch (IOException e) {
			// e.printStackTrace();
		} finally {
			try {
				stream.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return sb.toString();
	}

}
