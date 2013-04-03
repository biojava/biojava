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


		StringWriter b1 = new StringWriter();
		int pos = 0;
		for (String s : linkage.getPDBNameOfPotentialAtomsOnComponent1()){
			b1.append(s);
			pos++;
			if (pos < linkage.getPDBNameOfPotentialAtomsOnComponent1().size())
				b1.append(", ");

		}

		pos = 0;
		StringWriter b2 = new StringWriter();
		for (String s : linkage.getPDBNameOfPotentialAtomsOnComponent2()){
			b2.append(s);
			pos++;
			if (pos < linkage.getPDBNameOfPotentialAtomsOnComponent2().size())
				b2.append(", ");

		}

		xml.attribute("pdbNameComp1",
				b1.toString());

		xml.attribute("pdbNameComp2",
				b2.toString());

		xml.openTag("component1");		
		ComponentXMLConverter.toXML(c1, xml);
		xml.closeTag("component1");
		xml.openTag("component2");

		ComponentXMLConverter.toXML(c2, xml);
		xml.closeTag("component2");



		xml.closeTag("modificationLinkage");
	}



	public static ModificationLinkage fromXML(Node linkageElement, List<Component> components) {

		String name = linkageElement.getNodeName();
		if ( ! name.equals("modificationLinkage"))
			throw new RuntimeException("did not get modificationLinkage element, but " + name);

		Component c1 = getComponent("component1", linkageElement);
		Component c2 = getComponent("component2", linkageElement);

		int indexOfComponent1 = Integer.parseInt(getAttribute(linkageElement, "index1"));
		int indexOfComponent2 = Integer.parseInt(getAttribute(linkageElement, "index2"));

		String pdbNameComp1 = getAttribute(linkageElement, "pdbNameComp1");
		String pdbNameComp2 = getAttribute(linkageElement, "pdbNameComp2");

		if ( components == null )
			components = new ArrayList<Component>();
			
		if ( components.size() == 0) {

			components.add(c1);
			components.add(c2);
		}
		ModificationLinkage linkage = new ModificationLinkage(components, indexOfComponent1,pdbNameComp1, indexOfComponent2,pdbNameComp2);
		return linkage;
	}

	private static Component getComponent(String elementName, Node linkageElement) {

		String name = linkageElement.getNodeName();
		if ( ! name.equals("modificationLinkage"))
			throw new RuntimeException("did not get modificationLinkage element, but " + name);

		NodeList valList = linkageElement.getChildNodes();
		int numChildren  = valList.getLength();

		Component c = null;
		for ( int e =0; e< numChildren ; e++){
			Node  componentN = valList.item(e);
			if ( componentN.getNodeName().equals(elementName)) {
				NodeList children = componentN.getChildNodes();

				int numComps = children.getLength();

				for ( int i =0; i< numComps ; i++){
					Node compo = children.item(i);
					//System.out.println(compo.getNodeName());
					if ( compo.getNodeName().equals("component")) {
						c = ComponentXMLConverter.fromXML(compo)	;						
					}
				}
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
