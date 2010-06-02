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
 * Created on Jun 1, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava.bio.ptm.io;

import java.io.IOException;
import java.io.InputStream;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.biojava.bio.ptm.ModificationCategory;
import org.biojava.bio.ptm.ModificationOccurrenceType;
import org.biojava.bio.ptm.ProteinModification;

import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public final class ProteinModificationXmlReader {
	/**
	 * This is a utility class and thus cannot be instantialized.
	 */
	private ProteinModificationXmlReader() {}
	
	/**
	 * Register common PTMs.
	 */
	private static final String PTM_LIST_XML = "org/biojava/bio/ptm/ptm_list.xml";
	static {
		try {
			InputStream isXml = ProteinModificationXmlReader.class
				.getResourceAsStream(PTM_LIST_XML);
			registerProteinModificationFromXml(isXml);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Read protein modifications from XML file and register them.
	 * @param isXml {@link InputStream} of the XML file.
	 * @throws IOException if failed to read the XML file.
	 * @throws ParserConfigurationException if parse errors occur.
	 * @throws SAXException the {@link DocumentBuilder} cannot be created.
	 */
	public static void registerProteinModificationFromXml(InputStream isXml)
			throws IOException, ParserConfigurationException, SAXException {
		if (isXml==null)
			throw new IllegalArgumentException("Null argument.");
		
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = factory.newDocumentBuilder();
		Document doc = builder.parse(isXml);
		
		NodeList modNodes = doc.getElementsByTagName("Entry");		
		int modSize = modNodes.getLength();
		List<Node> nodes;
		for (int iMod=0; iMod<modSize; iMod++) {
			Node modNode = modNodes.item(iMod);
			Map<String,List<Node>> infoNodes = getChildNodes(modNode);
			
			// ID
			nodes = infoNodes.get("Id");
			if (nodes==null || nodes.size()!=1)
				throw new RuntimeException("Each modification must have exact " +
						"one <Id> field.");
			String id = nodes.get(0).getNodeValue();
			
			// modification category
			nodes = infoNodes.get("Type");
			if (nodes==null || nodes.size()!=1)
				throw new RuntimeException("Each modification must have exact " +
						"one <Type> field.");
			ModificationCategory cat = ModificationCategory.getByLabel(
					nodes.get(0).getNodeValue());
			if (cat==null)
				throw new RuntimeException(nodes.get(0).getTextContent()+
					" is not defined as an modification category.");
			
			// occurrence type
			nodes = infoNodes.get("Occurrence");
			if (nodes==null || nodes.size()!=1)
				throw new RuntimeException("Each modification must have exact " +
						"one <Occurrence> field.");
			ModificationOccurrenceType occType = ModificationOccurrenceType
				.getByLabel(nodes.get(0).getNodeValue());
			if (occType==null)
				throw new RuntimeException(nodes.get(0).getTextContent()+
					" is not defined as an modification occurence type.");
			
			ProteinModification.Builder modBuilder = ProteinModification
				.register(id, cat, occType);
			
			// description
			nodes = infoNodes.get("Description");
			if (nodes!=null && !nodes.isEmpty()) {
				modBuilder.description(nodes.get(0).getNodeValue());
			}
			
			// cross references
			nodes = infoNodes.get("CrossReference");
			if (nodes!=null) {
				for (Node node:nodes) {
					Map<String,List<Node>> xrefInfoNodes = getChildNodes(node);
					
					// source
					List<Node> xrefNode = xrefInfoNodes.get("Source");
					if (xrefNode==null || xrefNode.size()!=1)
						throw new RuntimeException("Error in XML file: " +
							"a cross reference must contain exactly one <Source> field.");
					String xrefDb = xrefNode.get(0).getNodeValue();
					
					// id
					xrefNode = xrefInfoNodes.get("Id");
					if (xrefNode==null || xrefNode.size()!=1)
						throw new RuntimeException("Error in XML file: " +
							"a cross reference must contain exactly one <Id> field.");
					String xrefId = xrefNode.get(0).getNodeValue();
					
					// name
					String xrefName = null;
					xrefNode = xrefInfoNodes.get("Name");
					if (xrefNode!=null && !xrefNode.isEmpty()) {
						xrefName = xrefNode.get(0).getNodeValue();
					}
					
					if (xrefDb.equals("PDBCC")) {
						modBuilder.pdbccId(xrefId).pdbccName(xrefName);
					} else if (xrefDb.equals("RESID")) {
						modBuilder.residId(xrefId).residName(xrefName);
					} else if (xrefDb.equals("PSI-MOD")) {
						modBuilder.psimodId(xrefId).psimodName(xrefName);
					}
				}
			} // end of cross references
			
			// formula
			nodes = infoNodes.get("Formula");
			if (nodes!=null && !nodes.isEmpty()) {
				modBuilder.formula(nodes.get(0).getNodeValue());
			}
			
			// components and atoms
			nodes = infoNodes.get("Components");
			if (nodes!=null) {
				int sizeComp = nodes.size();
				String[] comps = new String[sizeComp];
				for (int iComp=0; iComp<sizeComp; iComp++) {
					Node node = nodes.get(iComp);
					
					// keep track of the labels of components
					Map<String,Integer> mapLabelIndex = new HashMap<String,Integer>();

					Map<String,List<Node>> compInfoNodes = getChildNodes(node);
					
					// components
					List<Node> compNodes = compInfoNodes.get("Component");
					for (Node compNode:compNodes) {
						// comp label
						NamedNodeMap compNodeAttrs = compNode.getAttributes();
						Node labelNode = compNodeAttrs.getNamedItem("label");
						if (labelNode==null)
							throw new RuntimeException("Each component must have a label.");
						String label = labelNode.getNodeValue();
						
						if (mapLabelIndex.containsKey(label))
							throw new RuntimeException("Each component must have a unique label.");
						
						// comp PDBCC ID
						String compId = null;
						List<Node> compIdNodes = getChildNodes(compNode).get("Id");
						if (compIdNodes!=null) {
							for (Node compIdNode : compIdNodes) {
								NamedNodeMap compIdNodeAttr = compIdNode.getAttributes();
								Node compIdSource = compIdNodeAttr.getNamedItem("Source");
								if (compIdSource!=null && compIdSource.getNodeValue().equals("PDBCC")) {
									compId = compIdNode.getTextContent();
									break;
								}
							}
						}
						
						if (compId==null)
							throw new RuntimeException("Each component must have a PDBCC ID");
						
						comps[iComp] = compId;
						
						mapLabelIndex.put(label, iComp);
					}
					
					// bonds
					String[][] bonds = null;
					List<Node> bondNodes = compInfoNodes.get("Bond");
					if (bondNodes!=null) {
						bonds = new String[sizeComp][sizeComp];
						for (Node bondNode:bondNodes) {
							Map<String,List<Node>> bondChildNodes = getChildNodes(bondNode);
							if (bondChildNodes==null)
								throw new RuntimeException("Each bond must contain two atoms");
							
							List<Node> atomNodes = bondChildNodes.get("Atom");
							if (atomNodes==null || atomNodes.size()!=2)
								throw new RuntimeException("Each bond must contain two atoms");
							
							// atom 1
							NamedNodeMap atomNodeAttrs = atomNodes.get(0).getAttributes();
							Node labelNode = atomNodeAttrs.getNamedItem("component");
							if (labelNode==null)
								throw new RuntimeException("Each atom must on a component.");
							String comp1 = labelNode.getNodeValue();
							int iComp1 = mapLabelIndex.get(comp1);
							String atom1 = atomNodes.get(0).getNodeValue();
							
							// atom 2
							atomNodeAttrs = atomNodes.get(1).getAttributes();
							labelNode = atomNodeAttrs.getNamedItem("component");
							if (labelNode==null)
								throw new RuntimeException("Each atom must on a component.");
							String comp2 = labelNode.getNodeValue();
							int iComp2 = mapLabelIndex.get(comp2);
							String atom2 = atomNodes.get(1).getNodeValue();
							
							bonds[iComp1][iComp2] = atom1;
							bonds[iComp2][iComp1] = atom2;							
						}
					}
				}
			} // end of components
		}
	}
	
	/**
	 * Utility method to group child nodes by their names.
	 * @param parent parent node.
	 * @return Map from name to child nodes.
	 */
	private static Map<String,List<Node>> getChildNodes(Node parent) {
		if (parent==null)
			return null;
		
		Map<String,List<Node>> children = new HashMap<String,List<Node>>();
		
		NodeList nodes = parent.getChildNodes();
		int nNodes = nodes.getLength();
		for (int i=0; i<nNodes; i++) {
			Node node = nodes.item(i);
			String name = node.getNodeName();
			List<Node> namesakes = children.get(name);
			if (namesakes==null) {
				namesakes = new ArrayList<Node>();
				children.put(name, namesakes);
			}
			namesakes.add(node);
		}
		
		return children;
	}
}
