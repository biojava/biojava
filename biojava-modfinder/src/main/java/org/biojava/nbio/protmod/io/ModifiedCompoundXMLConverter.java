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

import org.biojava.nbio.protmod.ProteinModification;
import org.biojava.nbio.protmod.ProteinModificationRegistry;
import org.biojava.nbio.protmod.structure.*;
import org.biojava.nbio.core.util.PrettyXMLWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.*;

public class ModifiedCompoundXMLConverter {

	private static final Logger logger = LoggerFactory.getLogger(ModifiedCompoundXMLConverter.class);

	public static String toXML(ModifiedCompound mc) throws IOException{

		if ( mc == null) {
			logger.warn("ModifiedCompound == null! ");
			return "<modifiedCompound/>";
		}
		StringWriter out = new StringWriter();

		PrettyXMLWriter xml = new PrettyXMLWriter(new PrintWriter(out));

		ProteinModification protMod = mc.getModification();
		String modificationId = protMod==null?null:protMod.getId();

		xml.openTag("modifiedCompound");
		if ( modificationId != null) {
//			ProteinModificationXMLConverter.toXML(modification, xml);
			xml.openTag("proteinModification");
			xml.attribute("id", modificationId);
			xml.closeTag("proteinModification");
		}


		Set<StructureAtomLinkage > linkages = mc.getAtomLinkages();
		if ( linkages.size() > 0 ) {
			int pos = -1;
			for ( StructureAtomLinkage link: linkages){
				pos ++;
				xml.openTag("linkage");
				xml.attribute("pos", String.valueOf(pos));
				xml.attribute("total", String.valueOf(linkages.size()));
				StructureAtom atom1 = link.getAtom1();
				StructureAtom atom2 = link.getAtom2();
				double distance = link.getDistance();

				xml.attribute("distance", String.valueOf(distance));
				xml.openTag("atom1");
				StructureAtomXMLConverter.toXML(atom1,xml);
				xml.closeTag("atom1");
				xml.openTag("atom2");
				StructureAtomXMLConverter.toXML(atom2,xml);
				xml.closeTag("atom2");
				xml.closeTag("linkage");
			}
		} else {
			// no linkages, need to serialize the residues...
			xml.openTag("linkage");
			xml.closeTag("linkage");
			Set<StructureGroup> groups = mc.getGroups();
			for (StructureGroup group : groups) {
				StructureGroupXMLConverter.toXML(group, xml);
			}
		}




		xml.closeTag("modifiedCompound");
		return out.toString();
	}

	public static ModifiedCompound fromXML(String xml){
		ProteinModification modification = null;
		//Collection<StructureAtomLinkage> linkages = new ArrayList<StructureAtomLinkage>();
		StructureAtomLinkage[] linkages = null;
		List<StructureGroup> structureGroups = new ArrayList<StructureGroup>();
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

			NodeList listOfmodifications = doc.getElementsByTagName("modifiedCompound");
			//int numArrays = listOfArrays.getLength();
			// go over the blocks
			for(int modPos=0; modPos<listOfmodifications.getLength() ; modPos++)
			{

				Node modificationElement       = listOfmodifications.item(modPos);

				NodeList children = modificationElement.getChildNodes();

				int numChildren  = children.getLength();


				for ( int e =0; e< numChildren ; e++){
					Node  listOfConditions = children.item(e);

					if(!listOfConditions.hasAttributes()) continue;


					if ( listOfConditions.getNodeName().equals("proteinModification")) {
						//modification = ProteinModificationXMLConverter.fromXML(listOfConditions);
						String modId = getAttribute(listOfConditions, "id");
						modification = ProteinModificationRegistry.getById(modId);
						if (modification==null) {
							logger.warn("Error: no modification information.");
						}
					} else if ( listOfConditions.getNodeName().equals("linkage")) {
						double dist = Double.parseDouble(getAttribute(listOfConditions, "distance"));
						int pos = Integer.parseInt(getAttribute(listOfConditions,"pos"));
						int total = Integer.parseInt(getAttribute(listOfConditions,"total"));
						if ( linkages == null)
							linkages = new StructureAtomLinkage[total];

						StructureAtom atom1 = getAtom("atom1", listOfConditions);
						StructureAtom atom2 = getAtom("atom2",listOfConditions);
						StructureAtomLinkage linkage = new StructureAtomLinkage(atom1, atom2, dist);
						//linkages.add(linkage);
						linkages[pos] = linkage;
					} else if (listOfConditions.getNodeName().equals("structureGroup")) {
						StructureGroup group = StructureGroupXMLConverter.fromXML(listOfConditions);
						structureGroups.add(group);
//						logger.info("structureGroups size:" + structureGroups.size());
					}
				}
			}
		} catch (SAXParseException err)	{
			logger.error("** Parsing error, line: {}, uri: {}", err.getLineNumber (), err.getSystemId (), err);
		}
		catch (SAXException e) {
			logger.error("Exception: ", e);
		}
		catch (Throwable t) {
			logger.error("Exception: ", t);
		}


		if ( linkages != null) {
			Collection<StructureAtomLinkage> links = Arrays.asList(linkages);
			return new ModifiedCompoundImpl(modification, links);
		} else if ( structureGroups.size() == 1) {
			return new ModifiedCompoundImpl(modification, structureGroups.get(0));
		}
		return null;

	}



	private static StructureAtom getAtom(String elementName, Node n) {

		NodeList children = n.getChildNodes();

		int numChildren  = children.getLength();

		StructureAtom atom = null;
		for ( int e =0; e< numChildren ; e++){
			Node  atoms = children.item(e);

			if ( atoms.getNodeName().equals(elementName)) {
				NodeList child2 = atoms.getChildNodes();
				int numAtoms = child2.getLength();
				//logger.info("got " + numAtoms + " atoms");
				for ( int a=0;a< numAtoms; a++){
					Node atomNode = child2.item(a);
					if(!atomNode.hasAttributes()) continue;
					atom = StructureAtomXMLConverter.fromXML(atomNode);
					return atom;
				}

			}
		}
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
