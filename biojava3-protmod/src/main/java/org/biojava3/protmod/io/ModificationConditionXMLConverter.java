package org.biojava3.protmod.io;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import org.biojava3.core.util.PrettyXMLWriter;
import org.biojava3.protmod.Component;
import org.biojava3.protmod.ModificationCondition;
import org.biojava3.protmod.ModificationConditionImpl;
import org.biojava3.protmod.ModificationLinkage;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;




public class ModificationConditionXMLConverter {
	public static String toXML(ModificationCondition condition) throws IOException{
		StringWriter out = new StringWriter();

		PrettyXMLWriter xml = new PrettyXMLWriter(new PrintWriter(out));
		toXML(condition, xml);
		
		return out.toString();
	}
	
	public static void toXML(ModificationCondition condition, PrettyXMLWriter xml) throws IOException{
		
		xml.openTag("modificationCondition");
		
		List<Component>  components = condition.getComponents();
		
		for (Component c: components){
			
			ComponentXMLConverter.toXML(c, xml);
		}
		
		List<ModificationLinkage> linkages = condition.getLinkages();
		
		for ( ModificationLinkage linkage: linkages){
			ModificationLinkageXMLConverter.toXML(linkage, xml);
		}
				
		
		xml.closeTag("modificationCondition");
	}
	

	/** assumes that the XML parser is already at a modificationNode
	 * 
	 * @param conditionElement
	 * @return
	 */
	public static ModificationCondition fromXML( Node conditionElement) {
		
		String name = conditionElement.getNodeName();
		if ( ! name.equals("modificationCondition"))
			throw new RuntimeException("Did not get node modificationCondition, but " + name);
		
		NodeList valList = conditionElement.getChildNodes();
		int numChildren  = valList.getLength();
		
		
		
		List<ModificationLinkage> linkages = new ArrayList<ModificationLinkage>();
		List<Component> components = new ArrayList<Component>();
		
		for ( int e =0; e< numChildren ; e++){
			Node  modificationLinkages = valList.item(e);


			if ( modificationLinkages.getNodeName().equals("component")) {
				
				Component component = ComponentXMLConverter.fromXML(modificationLinkages);
				components.add(component);
				
			}

			if ( modificationLinkages.getNodeName().equals("modificationLinkage")) {
				
				ModificationLinkage linkage = ModificationLinkageXMLConverter.fromXML(modificationLinkages, components);
				
				linkages.add(linkage);
				
			}
			
			
		}
		
		return new ModificationConditionImpl(components, linkages);
		
	}
}
