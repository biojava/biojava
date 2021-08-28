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
package org.biojava.nbio.core.util;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;
import java.io.*;
import java.util.ArrayList;

import static org.biojava.nbio.core.sequence.io.util.IOUtils.close;
import static org.biojava.nbio.core.sequence.io.util.IOUtils.openFile;

/**
 * Helper methods to simplify boilerplate XML parsing code for  {@code}org.w3c.dom{@code} XML objects
 * @author Scooter
 */
public class XMLHelper {

	/**
	 * Creates a new element called {@code}elementName{@code} and adds it to {@code}parentElement{@code}
	 * @param parentElement
	 * @param elementName
	 * @return the new child element
	 */
	public static Element addChildElement(Element parentElement, String elementName) {
		Element childElement = parentElement.getOwnerDocument().createElement(elementName);
		parentElement.appendChild(childElement);
		return childElement;
	}

	/**
	 * Create a new, empty {@code}org.w3c.dom.Document{@code}
	 * @return a new {@code}org.w3c.dom.Document{@code}
	 * @throws ParserConfigurationException
	 */
	public static Document getNewDocument() throws ParserConfigurationException  {

		//Create instance of DocumentBuilderFactory
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		//Get the DocumentBuilder
		DocumentBuilder docBuilder = factory.newDocumentBuilder();
		//Create blank DOM Document
		Document doc = docBuilder.newDocument();
		return doc;
	}

	/**
	 * Given a path to an XML file, parses into an {@code}org.w3c.dom.Document{@code} 
	 * @param fileName path to a readable XML file
	 * @return
	 * @throws SAXException
	 * @throws IOException
	 * @throws ParserConfigurationException
	 */
	public static Document loadXML(String fileName) throws SAXException, IOException, ParserConfigurationException  {
		InputStream is = openFile(new File(fileName));
		Document doc = inputStreamToDocument(new BufferedInputStream(is));
		close(is);
		return doc;
	}

	/**
	 * Creates an {@code}org.w3c.dom.Document{@code} from the content of the {@code}inputStream{@code}
	 * @param inputStream
	 * @return a {@code}Document{@code}
	 * @throws SAXException
	 * @throws IOException
	 * @throws ParserConfigurationException
	 */
	public static Document inputStreamToDocument(InputStream inputStream) throws SAXException, IOException, ParserConfigurationException  {
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();

		DocumentBuilder db = dbf.newDocumentBuilder();

		Document doc = db.parse(inputStream);
		doc.getDocumentElement().normalize();

		return doc;
	}

	/**
	 * Given an {@code}org.w3c.dom.Document{@code}, writes it to the given {@code}outputStream{@code}
	 * @param document
	 * @param outputStream
	 * @throws TransformerException
	 */
	public static void outputToStream(Document document, OutputStream outputStream) throws TransformerException {
		// Use a Transformer for output
		TransformerFactory tFactory = TransformerFactory.newInstance();
		Transformer transformer = tFactory.newTransformer();
		//    transformer.setOutputProperty(OutputKeys.INDENT, "yes");

		DOMSource source = new DOMSource(document);
		StreamResult result = new StreamResult(outputStream);
		transformer.transform(source, result);
	}

	//static XPath xpath = XPathFactory.newInstance().newXPath();

	/**
	 * Given an element, searches upwards through ancestor Elements till the first Element
	 * matching the requests {@code}parentName{@code} is found.
	 * @param element The starting element
	 * @param parentName The tag name of the requested Element.
	 * @return The found element, or {@code}null{@code} if no matching element is found,
	 */
	public static Element selectParentElement(Element element, String parentName) {
		
	    Node parentNode =  element.getParentNode();
		if (parentNode == null) {
			return null;
		}
		// check that parent is actually an element, else return null
		// this is to prevent ClassCastExceptions if element's parent is not an Element.
		Element parentElement = null;
		if (Node.ELEMENT_NODE == parentNode.getNodeType()){
			parentElement = (Element)parentNode;
		} else {
			return null;
		}
		if (parentElement.getTagName().equals(parentName)) {
			return parentElement;
		}
		return selectParentElement(parentElement, parentName);
	}

	/**
	 * If {@code}xpathExpression{@code} is a plain string with no '/' characterr, this is 
	 * interpreted as a child element name to search for. 
	 * <b/>
	 * If {@code}xpathExpression{@code} is an XPath expression, this is evaluated and is assumed
	 * to identify a single element.
	 * @param element
	 * @param xpathExpression
	 * @return A single element or null if no match or the 1st match if matches more than 1
	 * @throws XPathExpressionException
	 */
	public static Element selectSingleElement(Element element, String xpathExpression) throws XPathExpressionException {
		if (element == null) {
			return null;
		}
		if (xpathExpression.indexOf("/") == -1) {
			NodeList nodeList = element.getChildNodes();
			for (int i = 0; i < nodeList.getLength(); i++) {
				Node node = nodeList.item(i);
				if (node.getNodeType() == Node.ELEMENT_NODE && node.getNodeName().equals(xpathExpression)) {
					return (Element) node;
				}
			}
			//  NodeList nodes = element.getElementsByTagName(xpathExpression);
			//  if (nodes.getLength() > 0) {
			//      return (Element) nodes.item(0);
			//  } else {
			return null;
			//  }
		} else {
			XPath xpath = XPathFactory.newInstance().newXPath();
			Element node = (Element) xpath.evaluate(xpathExpression, element, XPathConstants.NODE);
			return node;
		}
	}

	/**
	 * Gets a list of elements matching {@code}xpathExpression{@code}. If xpathExpression lacks
	 * a '/' character, only immediate children o {@code}element{@code} are searched over.
	 * <br/>
	 * If {@code}xpathExpression{@code} contains an '/' character, a full XPath search is made
	 * @param element
	 * @param xpathExpression
	 * @return A possibly empty but non-null {@code}ArrayList{@code}
	 * @throws XPathExpressionException
	 */
	public static ArrayList<Element> selectElements(Element element, String xpathExpression) throws XPathExpressionException {
		ArrayList<Element> resultVector = new ArrayList<Element>();
		if (element == null) {
			return resultVector;
		}
		if (xpathExpression.indexOf("/") == -1) {
			NodeList nodeList = element.getChildNodes();
			for (int i = 0; i < nodeList.getLength(); i++) {
				Node node = nodeList.item(i);
				if (node.getNodeType() == Node.ELEMENT_NODE && node.getNodeName().equals(xpathExpression)) {
					resultVector.add((Element) node);
				}
			}
		} else {
			XPath xpath = XPathFactory.newInstance().newXPath();
			NodeList nodes = (NodeList) xpath.evaluate(xpathExpression, element, XPathConstants.NODESET);


			for (int i = 0; i < nodes.getLength(); i++) {
				Node node = nodes.item(i);
				resultVector.add((Element) node);
			}
		}
		return resultVector;
	}
}
