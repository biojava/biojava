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
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Scanner;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.DocumentType;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.xml.sax.SAXException;


/**
 * A utility class for common {@link String} manipulation tasks.
 * All functions are static methods.
 * 
 * @author Amr AL-Hossary
 */
public class StringManipulationHelper  {

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
	 * @author andreas
	 * @param stream
	 * @return
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
	
	/**
	 * compares two strings for equality, line by line, ignoring any difference
	 * of end line delimiters contained within the 2 Strings. This method should
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
			line2 = scanner2.nextLine();
			if (! line1.equals(line2))
				return false;
		}
		if (scanner2.hasNextLine()) {
			return false;
		}

		return true;
	}
	
	
	public static boolean equalsToXml(String expected, String actual) {
		Document expectedDocument=null;
		Document actualDocument=null;
		try {
			DocumentBuilderFactory documentBuilderFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder documentBuilder = documentBuilderFactory.newDocumentBuilder();
			expectedDocument = documentBuilder.parse(new ByteArrayInputStream(expected.getBytes()));
			actualDocument = documentBuilder.parse(new ByteArrayInputStream(actual.getBytes()));
		} catch (ParserConfigurationException e) {
			e.printStackTrace();
			throw new RuntimeException("Couldn't Parse XML", e);
		} catch (SAXException e) {
			e.printStackTrace();
			throw new RuntimeException("Couldn't Parse XML", e);
		} catch (IOException e) {
			e.printStackTrace();
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
	
}
