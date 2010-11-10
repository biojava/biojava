package org.biojava3.protmod.io;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;


import org.biojava3.core.util.PrettyXMLWriter;
import org.biojava3.protmod.ModificationCategory;
import org.biojava3.protmod.ModificationCondition;
import org.biojava3.protmod.ModificationOccurrenceType;
import org.biojava3.protmod.ProteinModification;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

public class ProteinModificationXMLConverter {

	public static String toXML(ProteinModification modification) throws IOException{

		StringWriter out = new StringWriter();

		PrettyXMLWriter xml = new PrettyXMLWriter(new PrintWriter(out));
		toXML(modification, xml);

		return out.toString();
	}

	public static void toXML(ProteinModification modification, PrettyXMLWriter xml) throws IOException{
		xml.openTag("proteinModification");
		xml.attribute("id", modification.getId());
		xml.attribute("description", modification.getDescription());
		xml.attribute("category",modification.getCategory().toString());

		xml.attribute("occurenceType",modification.getOccurrenceType().toString());
		if ( modification.getResidId() != null)
			xml.attribute("residID", modification.getResidId());

		if ( modification.getPsimodId() != null)
			xml.attribute("psimodId", modification.getPsimodId());

		if ( modification.getPdbccId() != null) {
			xml.attribute("pdbccID", modification.getPdbccId());
		}
		if ( modification.getSystematicName() != null )
			xml.attribute("systematicName", modification.getSystematicName());

		
		ModificationCondition condition = modification.getCondition();
		if ( condition != null) {

			ModificationConditionXMLConverter.toXML(modification.getCondition(),xml);

		}
		

		xml.closeTag("proteinModification");


	}

	public static ProteinModification fromXML(String xml){

		ProteinModification protMod = null;

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


			//Element rootElement = doc.getDocumentElement();

			NodeList listOfmodifications = doc.getElementsByTagName("proteinModification");
			//int numArrays = listOfArrays.getLength();			
			// go over the blocks
			for(int modPos=0; modPos<listOfmodifications.getLength() ; modPos++)
			{

				Node rootElement       = listOfmodifications.item(modPos);

				protMod = fromXML(rootElement);

			}
		} catch (SAXParseException err) 
		{
			System.out.println ("** Parsing error" + ", line " 
					+ err.getLineNumber () + ", uri " + err.getSystemId ());
			System.out.println(" " + err.getMessage ());
		}
		catch (SAXException e)
		{
			Exception x = e.getException ();
			((x == null) ? e : x).printStackTrace ();
		}
		catch (Throwable t)
		{
			t.printStackTrace ();
		}

		return protMod;
	}

	public static ProteinModification fromXML(Node modificationElement) {
		
		String name = modificationElement.getNodeName();
		if ( ! name.equals("proteinModification"))
			throw new RuntimeException("Provided node is not a proteinModification node, but " + name);
		
		String id = getAttribute(modificationElement,"id");
		
		ProteinModification protMod = ProteinModification.getById(id);
		if ( protMod == null) {

			// we need to register it
			String categoryS = getAttribute(modificationElement,"category");
			String occurenceS = getAttribute(modificationElement, "occurenceType");
			ModificationCategory cat =  ModificationCategory.getByLabel(categoryS);
			ModificationOccurrenceType occType = ModificationOccurrenceType.getByLabel(occurenceS);

			//NodeList listOfConditions = modificationElement.getElementsByTagName("modificationCondition");

			NodeList valList = modificationElement.getChildNodes();
			int numChildren  = valList.getLength();

			
			ModificationCondition condition = null;
			
			for ( int e =0; e< numChildren ; e++){
				Node  listOfConditions = valList.item(e);

				if(!listOfConditions.hasAttributes()) continue;


				if ( listOfConditions.getNodeName().equals("modificationCondition")) {
															
					condition= ModificationConditionXMLConverter.fromXML(listOfConditions);
					
				}
			}
			ProteinModification.register(id, cat, occType, condition);
			protMod = ProteinModification.getById(id);
		}
		return protMod;
	}

	private static String getAttribute(Node node, String attr){
		if( ! node.hasAttributes()) 
			return null;

		NamedNodeMap atts = node.getAttributes();

		if ( atts == null)
			return null;

		Node att = atts.getNamedItem(attr);
		if ( att == null)
			return null;

		String value = att.getTextContent();

		return value;

	}




}
