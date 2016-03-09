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

import org.biojava.nbio.protmod.structure.StructureAtom;
import org.biojava.nbio.protmod.structure.StructureGroup;
import org.biojava.nbio.core.util.PrettyXMLWriter;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

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
