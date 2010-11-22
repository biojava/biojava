package org.biojava.bio.structure.align.xml;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.SortedSet;
import java.util.TreeSet;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;


import org.biojava3.core.util.PrettyXMLWriter;

import org.w3c.dom.Document;

import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

public class RepresentativeXMLConverter {


	public static final String toXML(SortedSet<String> representatives){
		StringWriter sw = new StringWriter();
		PrintWriter writer = new PrintWriter(sw);

		PrettyXMLWriter xml = new PrettyXMLWriter(writer);
		try {
			xml.openTag("representatives");

			for ( String repr : representatives){
				xml.openTag("pdbChain");
				xml.attribute("name", repr);								
				xml.closeTag("pdbChain");
			}
			xml.closeTag("representatives");
		} catch(IOException ex){
			ex.printStackTrace();
		}
		return sw.toString();
	}

	public static final SortedSet<String> fromXML(String xml){
		SortedSet<String> representatives = new TreeSet<String>();
		try {
			DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
			DocumentBuilder db = factory.newDocumentBuilder();
			InputSource inStream = new InputSource();
			inStream.setCharacterStream(new StringReader(xml));
			Document doc = db.parse(inStream);

			// normalize text representation
			doc.getDocumentElement().normalize();


			//Element rootElement = doc.getDocumentElement();

			NodeList listOfPairs = doc.getElementsByTagName("pdbChain");
			//int numArrays = listOfArrays.getLength();

			// go over the blocks
			for(int i=0; i<listOfPairs.getLength() ; i++)
			{
				Node pair       = listOfPairs.item(i);
				//NodeList valList = pair.getChildNodes();
				//int numChildren  = valList.getLength();

				NamedNodeMap map = pair.getAttributes();

				String name =  map.getNamedItem("name").getTextContent();				
				representatives.add(name);
			}

		} catch (Exception e){
			e.printStackTrace();
		}
		
		return representatives;
	}
}
