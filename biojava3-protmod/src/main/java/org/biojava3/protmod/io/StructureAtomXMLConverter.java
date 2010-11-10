package org.biojava3.protmod.io;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

import org.biojava3.core.util.PrettyXMLWriter;


import org.biojava3.protmod.structure.StructureAtom;
import org.biojava3.protmod.structure.StructureGroup;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class StructureAtomXMLConverter {

	public static String toXML(StructureAtom atom) throws IOException{

		StringWriter out = new StringWriter();

		PrettyXMLWriter xml = new PrettyXMLWriter(new PrintWriter(out));
		toXML(atom, xml);

		return out.toString();
	}

	public static void toXML(StructureAtom atom, PrettyXMLWriter xml) throws IOException{
		String name = atom.getAtomName();
		xml.openTag("structureAtom");
		xml.attribute("name", name);
		StructureGroup group = atom.getGroup();
		StructureGroupXMLConverter.toXML(group,xml);
		xml.closeTag("structureAtom");
	}

	public static StructureAtom fromXML(Node structureAtomElement){

		String name = structureAtomElement.getNodeName();
		if ( ! name.equals("structureAtom"))
			throw new RuntimeException("Node is not a structureAtom, but " +name);
		
		String atomName = getAttribute( structureAtomElement,"name");
		StructureGroup group = null;
		
		NodeList valList = structureAtomElement.getChildNodes();
		int numChildren  = valList.getLength();
	
		for ( int e =0; e< numChildren ; e++){
			Node  nodes = valList.item(e);

			if(!nodes.hasAttributes()) continue;


			if ( nodes.getNodeName().equals("structureGroup")) {
				group = StructureGroupXMLConverter.fromXML(nodes);
			}
		}
		StructureAtom atom = new StructureAtom(group, atomName);
		return atom;
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
