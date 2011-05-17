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

package org.biojava3.protmod.io;

import java.io.IOException;
import java.io.InputStream;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.biojava3.protmod.Component;
import org.biojava3.protmod.ModificationCategory;
import org.biojava3.protmod.ModificationCondition;
import org.biojava3.protmod.ModificationConditionImpl;
import org.biojava3.protmod.ModificationOccurrenceType;
import org.biojava3.protmod.ModificationLinkage;
import org.biojava3.protmod.ProteinModificationImpl;
import org.biojava3.protmod.ProteinModificationRegistry;

import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import org.xml.sax.SAXException;

/**
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public final class ProteinModificationXmlReader {
	/**
	 * This is a utility class and thus cannot be instantialized.
	 */
	private ProteinModificationXmlReader() {}
	
	/**
	 * Read protein modifications from XML file and register them.
	 * @param isXml {@link InputStream} of the XML file.
	 * @throws IOException if failed to read the XML file.
	 * @throws ParserConfigurationException if parse errors occur.
	 * @throws SAXException the {@link DocumentBuilder} cannot be created.
	 */
	public static void registerProteinModificationFromXml(InputStream isXml)
			throws IOException, ParserConfigurationException, SAXException {
		if (isXml==null) {
			throw new IllegalArgumentException("Null argument.");
		}
		
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
			if (nodes==null || nodes.size()!=1) {
				throw new RuntimeException("Each modification must have exact " +
						"one <Id> field.");
			}
			String id = nodes.get(0).getTextContent();
			
			// modification category
			nodes = infoNodes.get("Category");
			if (nodes==null || nodes.size()!=1) {
				throw new RuntimeException("Each modification must have exact " +
						"one <Category> field. See Modification "+id+".");
			}
			ModificationCategory cat = ModificationCategory.getByLabel(
					nodes.get(0).getTextContent());
			if (cat==null) {
				throw new RuntimeException(nodes.get(0).getTextContent()+
					" is not defined as an modification category." +
					" See Modification "+id+".");
			}
			
			// occurrence type
			nodes = infoNodes.get("Occurrence");
			if (nodes==null || nodes.size()!=1) {
				throw new RuntimeException("Each modification must have exact " +
						"one <Occurrence> field. See Modification "+id+".");
			}
			ModificationOccurrenceType occType = ModificationOccurrenceType
				.getByLabel(nodes.get(0).getTextContent());
			if (occType==null) {
				throw new RuntimeException(nodes.get(0).getTextContent()+
					" is not defined as an modification occurence type." +
					" See Modification "+id+".");
			}
			
			// condition
			ModificationCondition condition = null;
			{
				nodes = infoNodes.get("Condition");
				if (nodes==null || nodes.size()!=1) {
					throw new RuntimeException("Each modification must have exact " +
							"one <Condition> field. See Modification "+id+".");
				}
				
				Node compsNode = nodes.get(0);
				
				// keep track of the labels of component indices
				Map<String,Integer> mapLabelComp = new HashMap<String,Integer>();

				Map<String,List<Node>> compInfoNodes = getChildNodes(compsNode);
				
				// components
				List<Node> compNodes = compInfoNodes.get("Component");
				int sizeComp = compNodes.size();
				List<Component> comps = new ArrayList<Component>(sizeComp);
				for (int iComp=0; iComp<sizeComp; iComp++) {
					Node compNode = compNodes.get(iComp);
					// comp label
					NamedNodeMap compNodeAttrs = compNode.getAttributes();
					Node labelNode = compNodeAttrs.getNamedItem("component");
					if (labelNode==null) {
						throw new RuntimeException("Each component must have a label." +
								" See Modification "+id+".");
					}
					String label = labelNode.getTextContent();
					
					if (mapLabelComp.containsKey(label)) {
						throw new RuntimeException("Each component must have a unique label." +
								" See Modification "+id+".");
					}
					
					// comp PDBCC ID
					Set<String> compIds = new HashSet<String>();
					List<Node> compIdNodes = getChildNodes(compNode).get("Id");
					if (compIdNodes!=null) {
						for (Node compIdNode : compIdNodes) {
							NamedNodeMap compIdNodeAttr = compIdNode.getAttributes();
							Node compIdSource = compIdNodeAttr.getNamedItem("source");
							if (compIdSource!=null && compIdSource.getTextContent().equals("PDBCC")) {
								String strComps = compIdNode.getTextContent();
								if (strComps.isEmpty()) {
									throw new RuntimeException("Empty component." +
											" See Modification "+id+".");
								}
								compIds.addAll(Arrays.asList(strComps.split(",")));
							}
						}
					}
					
					if (compIds.isEmpty()) {
						throw new RuntimeException("Each component must have a PDBCC ID." +
								" See Modification "+id+".");
					}
					
					// terminal
					boolean nTerminal = false;
					boolean cTerminal = false;
					List<Node> compTermNode = getChildNodes(compNode).get("Terminal");
					if (compTermNode!=null) {
						if (compTermNode.size()!=1) {
							throw new RuntimeException("Only one <Terminal> condition is allowed for " +
									"each component. See Modification "+id+".");
						}
						String nc = compTermNode.get(0).getTextContent();
						if (nc.equals("N")) {
							nTerminal = true;
						} else if (nc.equals("C")) {
							cTerminal = true;
						} else {
							throw new RuntimeException("Only N or C is allowed for <Terminal>." +
									" See Modification "+id+".");
						}
					}

					// register
					Component comp = Component.of(compIds, nTerminal, cTerminal);
					comps.add(comp);						
					mapLabelComp.put(label, comps.size()-1);
				}
				
				// bonds
				List<Node> bondNodes = compInfoNodes.get("Bond");
				List<ModificationLinkage> linkages = null;
				if (bondNodes!=null) {
					int sizeBonds = bondNodes.size();
					linkages = new ArrayList<ModificationLinkage>(sizeBonds);
					for (int iBond=0; iBond<sizeBonds; iBond++) {
						Node bondNode = bondNodes.get(iBond);
						Map<String,List<Node>> bondChildNodes = getChildNodes(bondNode);
						if (bondChildNodes==null) {
							throw new RuntimeException("Each bond must contain two atoms" +
									" See Modification "+id+".");
						}
						
						List<Node> atomNodes = bondChildNodes.get("Atom");
						if (atomNodes==null || atomNodes.size()!=2) {
							throw new RuntimeException("Each bond must contain two atoms" +
									" See Modification "+id+".");
						}
						
						// atom 1
						NamedNodeMap atomNodeAttrs = atomNodes.get(0).getAttributes();
						Node compNode = atomNodeAttrs.getNamedItem("component");
						if (compNode==null) {
							throw new RuntimeException("Each atom must on a component." +
									" See Modification "+id+".");
						}
						String labelComp1 = compNode.getTextContent();
						int iComp1 = mapLabelComp.get(labelComp1);
						
						Node labelNode = atomNodeAttrs.getNamedItem("atom");
						String labelAtom1 = labelNode==null?null:labelNode.getTextContent();
						
						String atom1 = atomNodes.get(0).getTextContent();
						if (atom1.isEmpty()) {
							throw new RuntimeException("Each atom must have a name. Please use wildcard * if unknown." +
									" See Modification "+id+".");
						}
						List<String> potentialAtoms1 = Arrays.asList(atom1.split(","));
						
						// atom 2
						atomNodeAttrs = atomNodes.get(1).getAttributes();
						compNode = atomNodeAttrs.getNamedItem("component");
						if (compNode==null) {
							throw new RuntimeException("Each atom must on a component." +
									" See Modification "+id+".");
						}
						String labelComp2 = compNode.getTextContent();
						int iComp2 = mapLabelComp.get(labelComp2);

						labelNode = atomNodeAttrs.getNamedItem("atom");
						String labelAtom2 = labelNode==null?null:labelNode.getTextContent();
						
						String atom2 = atomNodes.get(1).getTextContent();
						if (atom2.isEmpty()) {
							throw new RuntimeException("Each atom must have a name. Please use wildcard * if unknown." +
									" See Modification "+id+".");
						}
						List<String> potentialAtoms2 = Arrays.asList(atom2.split(","));
						
						// add linkage
						ModificationLinkage linkage = new ModificationLinkage(comps,
								iComp1, potentialAtoms1, labelAtom1,
								iComp2, potentialAtoms2, labelAtom2);
						linkages.add(linkage);
					}
				}
				
				condition = new ModificationConditionImpl(comps, linkages);
			} // end of condition	
			
			ProteinModificationImpl.Builder modBuilder = 
				new ProteinModificationImpl.Builder(id, cat, occType, condition);
			
			// description
			nodes = infoNodes.get("Description");
			if (nodes!=null && !nodes.isEmpty()) {
				modBuilder.setDescription(nodes.get(0).getTextContent());
			}
			
			// cross references
			nodes = infoNodes.get("CrossReference");
			if (nodes!=null) {
				for (Node node:nodes) {
					Map<String,List<Node>> xrefInfoNodes = getChildNodes(node);
					
					// source
					List<Node> xrefNode = xrefInfoNodes.get("Source");
					if (xrefNode==null || xrefNode.size()!=1) {
						throw new RuntimeException("Error in XML file: " +
							"a cross reference must contain exactly one <Source> field." +
							" See Modification "+id+".");
					}
					String xrefDb = xrefNode.get(0).getTextContent();
					
					// id
					xrefNode = xrefInfoNodes.get("Id");
					if (xrefNode==null || xrefNode.size()!=1) {
						throw new RuntimeException("Error in XML file: " +
							"a cross reference must contain exactly one <Id> field." +
							" See Modification "+id+".");
					}
					String xrefId = xrefNode.get(0).getTextContent();
					
					// name
					String xrefName = null;
					xrefNode = xrefInfoNodes.get("Name");
					if (xrefNode!=null && !xrefNode.isEmpty()) {
						xrefName = xrefNode.get(0).getTextContent();
					}
					
					if (xrefDb.equals("PDBCC")) {
						modBuilder.setPdbccId(xrefId).setPdbccName(xrefName);
					} else if (xrefDb.equals("RESID")) {
						modBuilder.setResidId(xrefId).setResidName(xrefName);
					} else if (xrefDb.equals("PSI-MOD")) {
						modBuilder.setPsimodId(xrefId).setPsimodName(xrefName);
					}
				}
			} // end of cross references
			
			// formula
			nodes = infoNodes.get("Formula");
			if (nodes!=null && !nodes.isEmpty()) {
				modBuilder.setFormula(nodes.get(0).getTextContent());
			}
			
			// keywords
			nodes = infoNodes.get("Keyword");
			if (nodes!=null && !nodes.isEmpty()) {
				for (Node node : nodes) {
					modBuilder.addKeyword(node.getTextContent());
				}
			}
			
			ProteinModificationRegistry.register(modBuilder.build());
		}
	}
	
	/**
	 * Utility method to group child nodes by their names.
	 * @param parent parent node.
	 * @return Map from name to child nodes.
	 */
	private static Map<String,List<Node>> getChildNodes(Node parent) {
		if (parent==null)
			return Collections.emptyMap();
		
		Map<String,List<Node>> children = new HashMap<String,List<Node>>();
		
		NodeList nodes = parent.getChildNodes();
		int nNodes = nodes.getLength();
		for (int i=0; i<nNodes; i++) {
			Node node = nodes.item(i);
			if (node.getNodeType()!=Node.ELEMENT_NODE)
				continue;
			
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
