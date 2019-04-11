/**
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
 * Created on 2013-06-13
 * Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.nbio.structure.rcsb;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.io.InputStream;

/**
 * Package-level static utilities for parsing XML.
 * @author dmyerstu
 */
public class ReadUtils {

	private static final Logger logger = LoggerFactory.getLogger(ReadUtils.class);

	// this IS needed
	private static boolean documentBuilderFactorySet = false;

	/**
	 * @param s
	 * @return {@code s}, or null if {@code s} is the empty string
	 */
	static String toStr(String s) {
		if (s.isEmpty()) return null;
		return s;
	}

	/**
	 * @param stream
	 * @return A {@link NodeList} of top-level {@link Node Nodes} in {@code stream}.
	 * @throws IOException
	 */
	static NodeList getNodes(InputStream stream) throws IOException {

		if (!documentBuilderFactorySet) { // it's really stupid, but we have to do this
			System.setProperty("javax.xml.parsers.DocumentBuilderFactory",
					"com.sun.org.apache.xerces.internal.jaxp.DocumentBuilderFactoryImpl");
			documentBuilderFactorySet = true;
		}
		DocumentBuilderFactory builderFactory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = null;
		Document document = null;
		try {
			builder = builderFactory.newDocumentBuilder();
		} catch (ParserConfigurationException e) {
			logger.warn("Couldn't configure parser", e);
			stream.close();
			throw new IOException(e);
		}
		try {
			document = builder.parse(stream);
		} catch (SAXException e) {
			stream.close();
			throw new IOException(e);
		}
		Node root = document.getDocumentElement();
		return root.getChildNodes();
	}

	static Double toDouble(String s) {
		if (s.isEmpty()) return null;
		try {
			return Double.parseDouble(s);
		} catch (NumberFormatException e) {
			logger.warn(s + " is not a floating-point number", e);
		}
		return null;
	}

	static Integer toInt(String s) {
		if (s.isEmpty()) return null;
		try {
			return Integer.parseInt(s);
		} catch (NumberFormatException e) {
			logger.warn(s + " is not an integer", e);
		}
		return null;
	}

}
