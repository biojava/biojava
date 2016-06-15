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
 *
 * @author Scooter
 */
public class XMLHelper {

	static public Element addChildElement(Element parentElement, String elementName) {
		Element childElement = parentElement.getOwnerDocument().createElement(elementName);
		parentElement.appendChild(childElement);
		return childElement;
	}

	static public Document getNewDocument() throws ParserConfigurationException  {

		//Create instance of DocumentBuilderFactory
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		//Get the DocumentBuilder
		DocumentBuilder docBuilder = factory.newDocumentBuilder();
		//Create blank DOM Document
		Document doc = docBuilder.newDocument();
		return doc;
	}

	static public Document loadXML(String fileName) throws SAXException, IOException, ParserConfigurationException  {
		InputStream is = openFile(new File(fileName));
		Document doc = inputStreamToDocument(new BufferedInputStream(is));
		close(is);
		return doc;
	}

	static public Document inputStreamToDocument(InputStream inputStream) throws SAXException, IOException, ParserConfigurationException  {
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();

		DocumentBuilder db = dbf.newDocumentBuilder();

		Document doc = db.parse(inputStream);
		doc.getDocumentElement().normalize();

		return doc;
	}

	static public void outputToStream(Document document, OutputStream outputStream) throws TransformerException {
		// Use a Transformer for output
		TransformerFactory tFactory = TransformerFactory.newInstance();
		Transformer transformer = tFactory.newTransformer();
		//    transformer.setOutputProperty(OutputKeys.INDENT, "yes");

		DOMSource source = new DOMSource(document);
		StreamResult result = new StreamResult(outputStream);
		transformer.transform(source, result);


	}

	static public void outputToStream(Element document, OutputStream outputStream) throws TransformerException  {
		// Use a Transformer for output
		TransformerFactory tFactory = TransformerFactory.newInstance();
		Transformer transformer = tFactory.newTransformer();
		//     transformer.setOutputProperty(OutputKeys.INDENT, "yes");

		DOMSource source = new DOMSource(document);
		StreamResult result = new StreamResult(outputStream);
		transformer.transform(source, result);

	}
	//static XPath xpath = XPathFactory.newInstance().newXPath();

	static public Element selectParentElement(Element element, String parentName) {
		Element parentElement = (Element) element.getParentNode();
		if (parentElement == null) {
			return null;
		}
		if (parentElement.getTagName().equals(parentName)) {
			return parentElement;
		}
		return selectParentElement(parentElement, parentName);
	}

	static public Element selectSingleElement(Element element, String xpathExpression) throws XPathExpressionException {
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

	static public ArrayList<Element> selectElements(Element element, String xpathExpression) throws XPathExpressionException {
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
