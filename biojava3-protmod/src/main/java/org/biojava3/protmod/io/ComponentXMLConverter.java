package org.biojava3.protmod.io;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.HashSet;
import java.util.Set;

import org.biojava3.core.util.PrettyXMLWriter;
import org.biojava3.protmod.Component;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;


public class ComponentXMLConverter {
	public static String toXML(Component component) throws IOException{
		StringWriter out = new StringWriter();

		PrettyXMLWriter xml = new PrettyXMLWriter(new PrintWriter(out));
		toXML(component, xml);

		return out.toString();
	}

	public static void toXML(Component component, PrettyXMLWriter xml) throws IOException{
		xml.openTag("component");
		
		xml.attribute("nTerminal" , component.isNTerminal()+"");
		xml.attribute("cTerminal", component.isCTerminal()+"");
		for (String pdbccId : component.getPdbccIds()){
			xml.openTag("pdbccID");
			xml.attribute("id", pdbccId);
			xml.closeTag("pdbccID");					
		}
		
		xml.closeTag("component");
	}

	public static Component fromXML(String xml){
		return null;
	}

	public static Component fromXML(Node componentN) {
		
		String name = componentN.getNodeName();
		if ( ! name.equals("component"))
			throw new RuntimeException("did not get component element, but " + name);
		
		String type = getAttribute(componentN, "type");
		String nTerminalS = getAttribute(componentN, "nTerminal");
		String cTerminalS = getAttribute(componentN, "cTerminal");
		
		boolean isNTerminal = Boolean.parseBoolean(nTerminalS); 
		boolean isCTerminal = Boolean.parseBoolean(cTerminalS);
                
		Set<String>pdbccIds = new HashSet<String>();

		NodeList valList = componentN.getChildNodes();
		int numChildren  = valList.getLength();


		for ( int e =0; e< numChildren ; e++){
			Node  pdbccN = valList.item(e);

			if(!pdbccN.hasAttributes()) continue;


			if ( pdbccN.getNodeName().equals("pdbccID")) {
				String id = getAttribute(pdbccN, "id");
				pdbccIds.add(id);
			}

		}

		Component c = Component.of(pdbccIds, isNTerminal, isCTerminal);
		return c;

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
