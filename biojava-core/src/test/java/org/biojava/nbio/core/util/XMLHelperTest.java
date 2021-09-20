package org.biojava.nbio.core.util;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.xpath.XPathExpressionException;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.w3c.dom.DOMException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

class XMLHelperTest {

    // simple XML used in most of the tests:
    final String TEST_XML = "<root><list><a id=\"1\"/> <a id=\"2\"/> </list></root>";

    @Test
    @DisplayName("Create empty w3dom Document")
    void getNewDocument() throws ParserConfigurationException {
        Document d = XMLHelper.getNewDocument();
        assertNotNull(d);
        assertFalse(d.hasChildNodes());
        assertNull(d.getInputEncoding());
    }

    @Test
    @DisplayName("Create empty w3dom Document")
    void addChildDocument() throws ParserConfigurationException, DOMException {

        Document d = createDocumentWithRootElement();
        Element root = (Element) d.getChildNodes().item(0);

        Element added = XMLHelper.addChildElement(root, "myelement");
        assertNotNull(added);
        assertEquals(root, added.getParentNode());
        assertEquals(added, root.getChildNodes().item(0));
    }

    @Test
    void inputStreamToDocument() throws SAXException, IOException, ParserConfigurationException {
        Document doc = readTestDoc();
        assertParsedDocument(doc);
    }

    Document readTestDoc() throws SAXException, IOException, ParserConfigurationException {
        ByteArrayInputStream bArrayInputStream = new ByteArrayInputStream(TEST_XML.getBytes());
        return XMLHelper.inputStreamToDocument(bArrayInputStream);
    }

    @Test
    void fileToDocument() throws IOException, SAXException, ParserConfigurationException {
        File tmpFile = File.createTempFile("xml", ".xml");
        Files.write(Paths.get(tmpFile.getAbsolutePath()), TEST_XML.getBytes());
        Document doc = XMLHelper.loadXML(tmpFile.getAbsolutePath());
        assertParsedDocument(doc);
    }

    @Test
    void documentToOutputStream() throws SAXException, IOException, ParserConfigurationException, TransformerException {
        ByteArrayOutputStream baos = new ByteArrayOutputStream(100);
        Document doc = readTestDoc();
        XMLHelper.outputToStream(doc, baos);
        assertEquals("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" + TEST_XML,
                new String(baos.toByteArray()));
    }

    @Test
    void selectParentElement() throws SAXException, IOException, ParserConfigurationException {
        Document doc = readTestDoc();
        
        // get a a grandchild element
        NodeList nodes = doc.getElementsByTagName("a");
        
        // can get root node
        Element el = (Element) nodes.item(0);
        Element root = XMLHelper.selectParentElement(el, "root");
        assertNotNull(root);
        
        // non-existing element or if is root node returns null
        assertNull(XMLHelper.selectParentElement(el, "notexisting"));
        assertNull(XMLHelper.selectParentElement(root, "notexisting"));
    }

    @Nested
    class SelectSingleElement {
        @Test
        void selectSingleElement()
                throws SAXException, IOException, ParserConfigurationException, XPathExpressionException {
            Document doc = readTestDoc();
            Element root = (Element) doc.getElementsByTagName("root").item(0);

            // not direct child
            assertNull(XMLHelper.selectSingleElement(root, "a"));

            // direct child
            assertNotNull(XMLHelper.selectSingleElement(root, "list"));

            // xpath match
            Element found = XMLHelper.selectSingleElement(root, "/root/list/a[@id = \"2\"]");
            assertNotNull(found);
            assertEquals("2", found.getAttribute("id"));

            // xpath no match
            Element Notfound = XMLHelper.selectSingleElement(root, "/root/list/a[@id = \"45\"]");
            assertNull(Notfound);

            // xpath returning multiple elements returns 1st element
            Element mult = XMLHelper.selectSingleElement(root, "/root/list/a");
            assertNotNull(mult);
        }

        @Test
        void invalidInput() throws XPathExpressionException {
            assertNull(XMLHelper.selectSingleElement(null, "root"));
        }
    }

    @Nested
    class SelectElements {

       private Document doc = null;
       private Element root = null;

        @BeforeEach
        void before() throws SAXException, IOException, ParserConfigurationException {
             doc = readTestDoc();
             root = (Element) doc.getElementsByTagName("root").item(0);
        }

        @Test
        void selectMultipleElementsWithXPath()
                throws  XPathExpressionException {
            ArrayList<Element> selected = XMLHelper.selectElements(root, "/root/list/a");
            assertEquals(2, selected.size());
        }

        @Test
        void selectMultipleElementsWithXPathSearchesWholeTree()
                throws  XPathExpressionException {
            Element a1  = (Element) doc.getElementsByTagName("a").item(0);
            
            ArrayList<Element> selected = XMLHelper.selectElements(a1, "/root");
            assertEquals(1, selected.size());
            assertEquals("root", selected.get(0).getTagName());
        }

        @Test
        void selectBySimpleTagName() throws XPathExpressionException {
            // search by simple name doesn't search past children
            assertEquals(0, XMLHelper.selectElements(root, "a").size());
            Element list = (Element) doc.getElementsByTagName("list").item(0);

            // 'list' is immediate parent of 'a'
            assertEquals(2, XMLHelper.selectElements(list, "a").size());
        }

        @Test
        void invalidInputtoSelectElements() throws XPathExpressionException {
            assertEquals(0, XMLHelper.selectElements(null, "root").size());
        }
    }

    void assertParsedDocument(Document doc) {
        assertNotNull(doc);
        assertEquals(2, doc.getElementsByTagName("a").getLength());
        assertEquals(1, doc.getElementsByTagName("list").getLength());
    }

    Document createDocumentWithRootElement() throws ParserConfigurationException {
        Document doc = XMLHelper.getNewDocument();
        Element root = doc.createElement("root");
        doc.appendChild(root);
        return doc;
    }
}