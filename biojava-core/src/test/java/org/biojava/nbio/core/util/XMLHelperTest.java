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
import java.nio.file.Path;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerException;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.w3c.dom.DOMException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

class XMLHelperTest {
    
    final String TEST_XML = "<root><list><a id=\"1\"/> <a id=\"2\"/> </list></root>";
    
     @Test
    @DisplayName("Create empty w3dom Document")
    void getNewDocument() throws ParserConfigurationException{
        Document d = XMLHelper.getNewDocument();
        assertNotNull(d);
        assertFalse(d.hasChildNodes());
        assertNull(d.getInputEncoding());
    }


    @Test
    @DisplayName("Create empty w3dom Document")
    void addChildDocument() throws ParserConfigurationException, DOMException{
             
        Document d = createDocumentWithRootElement();
        Element root = (Element)d.getChildNodes().item(0);

        Element added = XMLHelper.addChildElement(root, "myelement");
        assertNotNull(added);
        assertEquals(root, added.getParentNode());
        assertEquals(added, root.getChildNodes().item(0));
    }

    @Test
    void inputStreamToDocument() throws SAXException, IOException, ParserConfigurationException{
        Document doc = readTestDoc();
        assertParsedDocument(doc);
    }

    Document readTestDoc() throws SAXException, IOException, ParserConfigurationException{
        ByteArrayInputStream bArrayInputStream = new ByteArrayInputStream(TEST_XML.getBytes());
        return XMLHelper.inputStreamToDocument(bArrayInputStream);
    }

    @Test 
    void  fileToDocument () throws IOException, SAXException, ParserConfigurationException{
        File tmpFile = File.createTempFile("xml", ".xml");
        Files.write(Path.of(tmpFile.getAbsolutePath()), TEST_XML.getBytes());
        Document doc = XMLHelper.loadXML(tmpFile.getAbsolutePath());
        assertParsedDocument(doc);
    }

    @Test
    void documentToOutputStream() throws SAXException, IOException, ParserConfigurationException, TransformerException{
        ByteArrayOutputStream baos = new ByteArrayOutputStream(100);
        Document doc = readTestDoc();
        XMLHelper.outputToStream(doc, baos);
        assertEquals("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" +
             TEST_XML, new String(baos.toByteArray()));
    }

    @Test
    void selectParentElement() throws SAXException, IOException, ParserConfigurationException{
        Document doc = readTestDoc();
        NodeList nodes = doc.getElementsByTagName("a");
        Element el = (Element)nodes.item(0);
        Element root = XMLHelper.selectParentElement(el, "root");
        assertNotNull(root);
        assertNull(XMLHelper.selectParentElement(el, "notexisting"));
        assertNull(XMLHelper.selectParentElement(root, "notexisting"));

    }
    void assertParsedDocument (Document doc){
        assertNotNull(doc);
        assertEquals(2, doc.getElementsByTagName("a").getLength());
        assertEquals(1, doc.getElementsByTagName("list").getLength());
    }

    Document createDocumentWithRootElement() throws ParserConfigurationException{
        Document d = XMLHelper.getNewDocument();
        Element root = d.createElement("root");
        d.appendChild(root);
        return d;
    }
}