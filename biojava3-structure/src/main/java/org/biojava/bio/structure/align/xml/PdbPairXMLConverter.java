package org.biojava.bio.structure.align.xml;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.SortedSet;
import java.util.TreeSet;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.biojava.bio.structure.align.client.PdbPair;
import org.biojava.bio.structure.align.fatcat.FatCatRigid;
import org.biojava3.core.util.PrettyXMLWriter;

import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

public class PdbPairXMLConverter {
	
	
	public static final String DEFAULT_METHOD_NAME = FatCatRigid.algorithmName;
	
	public static PdbPairsMessage convertXMLtoPairs(String xml) {
		SortedSet<PdbPair>  pairs = new TreeSet<PdbPair>();
		PdbPairsMessage message = new PdbPairsMessage();
		try
		{
			//Convert string to XML document
			DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
			DocumentBuilder db = factory.newDocumentBuilder();
			InputSource inStream = new InputSource();
			inStream.setCharacterStream(new StringReader(xml));
			Document doc = db.parse(inStream);

			// normalize text representation
			doc.getDocumentElement().normalize();

			NodeList algorithms =  doc.getElementsByTagName("pairs");
			//System.out.println(algorithms.getLength());
			for(int i=0; i<algorithms.getLength() ; i++) {
				
				Node algo       = algorithms.item(i);
				
				NamedNodeMap map = algo.getAttributes();
				
				String name = map.getNamedItem("algorithm").getTextContent();
				
				if ( name != null) {
					
					message.setMethod(name);
				
				}
				
			}
			

			NodeList listOfPairs = doc.getElementsByTagName("pair");
			//System.out.println(listOfPairs.getLength());

			// go over the blocks
			for(int i=0; i<listOfPairs.getLength() ; i++)
			{
				Node pair       = listOfPairs.item(i);
				//NodeList valList = pair.getChildNodes();
				//int numChildren  = valList.getLength();

				NamedNodeMap map = pair.getAttributes();

				String name1 =  map.getNamedItem("name1").getTextContent();
				String name2 =  map.getNamedItem("name2").getTextContent();
				PdbPair pdbPair = new PdbPair(name1, name2);
				pairs.add(pdbPair);
			}

		} catch (Exception e){
			e.printStackTrace();
		}
		message.setPairs(pairs);
		return message;
	}

	public static String convertPairsToXML(SortedSet<PdbPair> pairs, String method){
		StringWriter sw = new StringWriter();
		PrintWriter writer = new PrintWriter(sw);

		if (method == null){
			method = DEFAULT_METHOD_NAME;
		}
		
		
		PrettyXMLWriter xml = new PrettyXMLWriter(writer);
		try {
			
			
			xml.openTag("pairs");
			xml.attribute("algorithm", method);
			for ( PdbPair pair : pairs){
				xml.openTag("pair");
				xml.attribute("name1", pair.getName1());
				xml.attribute("name2", pair.getName2());				
				xml.closeTag("pair");
			}
			xml.closeTag("pairs");
		} catch(IOException ex){
			ex.printStackTrace();
		}

		return sw.toString();
	}


}
