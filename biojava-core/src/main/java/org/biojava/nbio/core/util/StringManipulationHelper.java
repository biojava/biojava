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
 * Author: Amr ALHOSSARY, Richard Adams
 *
 */
package org.biojava.nbio.core.util;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Scanner;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.DocumentType;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.xml.sax.SAXException;


/**
 * A utility class for common {@link String} manipulation tasks.
 * All functions are static methods.
 *
 * @author Amr ALHOSSARY
 * @author Richard Adams
 */
public class StringManipulationHelper  {

	private final static Logger logger = LoggerFactory.getLogger(StringManipulationHelper.class);

	/**
	 * we are using Unix endline here, since this is used for testing XML and it
	 * is part of the XML recommendations: <a href
	 * ="http://www.w3.org/TR/REC-xml/#sec-line-ends"
	 * >http://www.w3.org/TR/REC-xml/#sec-line-ends</a>
	 */
	private static final String UNIX_NEWLINE = "\n";

	private StringManipulationHelper() {
		// to prevent instantiation
	}

	/**
	 * Converts an InputStream of text to a String, closing the stream
	 * before returning.
	 * <ul>
	 * <li> Newlines are converted to Unix newlines (\n)
	 * <li> Default charset encoding is used to read the stream.
	 * <li> Any IOException reading the stream is 'squashed' and not made
	 *   available to caller
	 * <li> An additional newline is appended at the end of the string.
	 * <ul>
	 * @author andreas
	 * @param stream
	 * @return a possibly empty but non-null String
	 */
	public static String convertStreamToString(InputStream stream) {
		BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
		StringBuilder sb = new StringBuilder();

		String line = null;
		try {
			while ((line = reader.readLine()) != null) {
			    sb.append(line).append(UNIX_NEWLINE);
			}
		} catch (IOException e) {
			// logger.error("Exception: ", e);
		} finally {
			try {
				stream.close();
			} catch (IOException e) {
				logger.error("Exception: ", e);
			}
		}
		return sb.toString();
	}

	/**
	 * Compares two strings in a case-sensitive manner for equality, line by line, ignoring any difference
	 * of end line delimiters contained within the 2 Strings.
	 * <br/>
	 * This method should
	 * be used if and only if two Strings are considered identical when all nodes
	 * are identical including their relative order. Generally useful when
	 * asserting identity of <b>automatically regenerated</b> XML or PDB.
	 *
	 * @param expected
	 * @param actual
	 */
	public static boolean equalsToIgnoreEndline(String expected, String actual) {
		if (expected == null && actual == null) {
			return true;
		}
		if (expected != null ^ actual != null) {
			return false;
		}
		Scanner scanner1 = new Scanner(expected);
		Scanner scanner2 = new Scanner(actual);
		String line1, line2;
		while (scanner1.hasNextLine()) {
			line1 = scanner1.nextLine();
			if(scanner2.hasNextLine()) {
				line2 = scanner2.nextLine();
				if (! line1.equals(line2)) {
					closeScanners(scanner1, scanner2);
					return false;
				}
			} else {
				closeScanners(scanner1, scanner2);
				return false;
			}
		}
		if (scanner2.hasNextLine()) {
			closeScanners(scanner1, scanner2);
			return false;
		}

		closeScanners(scanner1, scanner2);
		return true;
	}

	private static void closeScanners(Scanner s1, Scanner s2) {
		s1.close();
		s2.close();
	}

	/**
	 * This method is not implemented or used, never returns true
	 * and should probably be removed.
	 * @param expected
	 * @param actual
	 * @return
	 * @throws UnsupportedOperationException in most cases
	 */
	public static boolean equalsToXml(String expected, String actual) {
		Document expectedDocument=null;
		Document actualDocument=null;
		try {
			DocumentBuilderFactory documentBuilderFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder documentBuilder = documentBuilderFactory.newDocumentBuilder();
			expectedDocument = documentBuilder.parse(new ByteArrayInputStream(expected.getBytes()));
			actualDocument = documentBuilder.parse(new ByteArrayInputStream(actual.getBytes()));
		} catch (ParserConfigurationException e) {
			logger.error("Exception: ", e);
			throw new RuntimeException("Couldn't Parse XML", e);
		} catch (SAXException e) {
			logger.error("Exception: ", e);
			throw new RuntimeException("Couldn't Parse XML", e);
		} catch (IOException e) {
			logger.error("Exception: ", e);
			throw new RuntimeException("Couldn't Parse XML", e);
		}
		final DocumentType doctype1 = expectedDocument.getDoctype();
		final DocumentType doctype2 = actualDocument.getDoctype();
		if (doctype1==null ^ doctype2 == null) {
			return false;
		}else if (doctype1!= null /*&& doctype2 != null*/) {
			NamedNodeMap expectedNotations = doctype1.getNotations();
			NamedNodeMap actualNotations = doctype2.getNotations();
			if (expectedNotations.getLength() == actualNotations.getLength()) {
				for (int i = 0; i < expectedNotations.getLength(); i++) {
					Node node= expectedNotations.item(i);
					node.isEqualNode(null);
				}
			}else{
				return false;
			}

		}

		throw new UnsupportedOperationException("not yet implemented");
	}

	/**
	 * Adds padding to left of supplied string
	 * @param s The String to pad
	 * @param n an integer >= 1
	 * @return The left-padded string. 
	 * @throws IllegalArgumentException if n <= 0
	 */
	public static String padLeft(String s, int n) {
		validatePadding(n);
	    return String.format("%1$" + n + "s", s);
	}

	/**
	 * Adds padding to right of supplied string
	 * @param s The String to pad
	 * @param n an integer >= 1
	 * @return The right-padded string. 
	 * @throws IllegalArgumentException if n <= 0
	 */
	public static String padRight(String s, int n) {
		validatePadding(n);
	    return String.format("%1$-" + n + "s", s);
	}

	private static void validatePadding(int n) {
		if (n <=0 ) {
			throw new IllegalArgumentException("padding must be >= 1");
		}
	}

	/**
	 * Joins Strings together with a delimiter to a single
	 * @param s An {@link Iterable} of Strings
	 * @param delimiter
	 * @return
	 */
	public static String join(Collection<String> s, String delimiter) {
		if (s==null) return "";
		return s.stream().collect(Collectors.joining(delimiter));
	}

}
