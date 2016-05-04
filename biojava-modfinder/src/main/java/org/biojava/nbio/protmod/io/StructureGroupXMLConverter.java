/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.protmod.io;

import org.biojava.nbio.protmod.structure.StructureGroup;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.core.util.PrettyXMLWriter;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

import java.io.IOException;

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
		resNum.setChainName(chainID);
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
