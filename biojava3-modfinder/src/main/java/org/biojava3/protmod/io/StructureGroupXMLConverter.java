package org.biojava3.protmod.io;

import java.io.IOException;

import org.biojava.bio.structure.ResidueNumber;
import org.biojava3.core.util.PrettyXMLWriter;
import org.biojava3.protmod.structure.StructureGroup;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

public class StructureGroupXMLConverter {

	public static void toXML(StructureGroup group, PrettyXMLWriter xml) throws IOException{
		
		xml.openTag("structureGroup");
		xml.attribute("chainID", group.getChainId());
		xml.attribute("pdbName", group.getPDBName());
		if ( group.getInsCode() != null)
			xml.attribute("insCode",group.getInsCode()+"");
		xml.attribute("residueNr", group.getResidueNumber()+"");
		xml.attribute("isAminoAcid", Boolean.toString(group.isAminoAcid()));
		xml.closeTag("structureGroup");
	}

	public static StructureGroup fromXML(Node n) {
		
		
		String chainID = getAttribute(n, "chainID");
		String pdbName = getAttribute(n, "pdbName");
		String insCode = getAttribute(n, "insCode");
		String resN  = getAttribute(n, "residueNr");
		String isAminoAcid = getAttribute(n,"isAminoAcid");
		
		ResidueNumber resNum = new ResidueNumber();
		resNum.setChainId(chainID);
		if ( ( insCode != null) && (! insCode.equals("null")) && insCode.length() == 1)
		resNum.setInsCode(insCode.charAt(0));
		resNum.setSeqNum(Integer.parseInt(resN));
		
		StructureGroup g = new StructureGroup(resNum, pdbName, Boolean.valueOf(isAminoAcid));
		return g;
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
