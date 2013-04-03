package org.biojava3.protmod.io;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import java.util.List;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.biojava3.core.util.PrettyXMLWriter;

import org.biojava3.protmod.ProteinModification;
import org.biojava3.protmod.ProteinModificationRegistry;
import org.biojava3.protmod.structure.ModifiedCompound;
import org.biojava3.protmod.structure.ModifiedCompoundImpl;
import org.biojava3.protmod.structure.StructureAtom;
import org.biojava3.protmod.structure.StructureAtomLinkage;
import org.biojava3.protmod.structure.StructureGroup;

import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

public class ModifiedCompoundXMLConverter {


	public static String toXML(ModifiedCompound mc) throws IOException{
		
		if ( mc == null) { 
			System.err.println("ModifiedCompound == null! ");
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
				xml.attribute("pos", pos+"");
				xml.attribute("total", linkages.size()+"");
				StructureAtom atom1 = link.getAtom1();
				StructureAtom atom2 = link.getAtom2();
				double distance = link.getDistance();

				xml.attribute("distance", distance+"");
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
							System.err.println("Error: no modification information.");
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
//						System.out.println("structureGroups size:" + structureGroups.size());
					}
				}
			}


		} catch (SAXParseException err) 
		{
			System.out.println ("** Parsing error" + ", line " 
					+ err.getLineNumber () + ", uri " + err.getSystemId ());
			System.out.println(" " + err.getMessage ());
		}
		catch (SAXException e)
		{
			Exception x = e.getException ();
			((x == null) ? e : x).printStackTrace ();
		}
		catch (Throwable t)
		{
			t.printStackTrace ();
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
				//System.out.println("got " + numAtoms + " atoms");
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
