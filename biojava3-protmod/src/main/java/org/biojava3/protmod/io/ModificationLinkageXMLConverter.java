package org.biojava3.protmod.io;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import org.biojava3.core.util.PrettyXMLWriter;

import org.biojava3.protmod.Component;
import org.biojava3.protmod.ModificationLinkage;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class ModificationLinkageXMLConverter {
	public static String toXML(ModificationLinkage linkage) throws IOException{
		StringWriter out = new StringWriter();

		PrettyXMLWriter xml = new PrettyXMLWriter(new PrintWriter(out));
		toXML(linkage, xml);
		
		return out.toString();
	}
	
	public static void toXML(ModificationLinkage linkage, PrettyXMLWriter xml) throws IOException{
		xml.openTag("modificationLinkage");
		
		
		Component c1 = linkage.getComponent1();
		Component c2 = linkage.getComponent2();
		
		xml.attribute("index1", linkage.getIndexOfComponent1()+"");
		xml.attribute("index2", linkage.getIndexOfComponent2()+"");
				
		xml.openTag("component1");		
		ComponentXMLConverter.toXML(c1, xml);
		xml.closeTag("component1");
		xml.openTag("component2");
		
		ComponentXMLConverter.toXML(c2, xml);
		xml.closeTag("component2");
		
		xml.closeTag("modificationLinkage");
	}
	

	
	public static ModificationLinkage fromXML(Node linkageElement) {
		
		Component c1 = getComponent("component1", linkageElement);
		Component c2 = getComponent("component2", linkageElement);
		
		int indexOfComponent1 = Integer.parseInt(getAttribute(linkageElement, "index1"));
		int indexOfComponent2 = Integer.parseInt(getAttribute(linkageElement, "index2"));
		
		List<Component> components = new ArrayList<Component>();
		
		components.add(c1);
		components.add(c2);
		
		ModificationLinkage linkage = new ModificationLinkage(components, indexOfComponent1, indexOfComponent2);
		return linkage;
	}
	
	private static Component getComponent(String elementName, Node linkageElement) {
		NodeList valList = linkageElement.getChildNodes();
		int numChildren  = valList.getLength();
		
		Component c = null;
		for ( int e =0; e< numChildren ; e++){
			Node  componentN = valList.item(e);

			if(!componentN.hasAttributes()) continue;


			if ( componentN.getNodeName().equals(elementName)) {
				c = ComponentXMLConverter.fromXML(componentN)	;
			}
			
		}
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
