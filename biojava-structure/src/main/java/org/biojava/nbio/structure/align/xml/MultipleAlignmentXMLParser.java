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
package org.biojava.nbio.structure.align.xml;

import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsemble;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsembleImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.ScoresCache;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

/**
 * Parse an XML file representing a {@link MultipleAlignmentEnsemble}, so
 * that the original alignment can be recovered.
 * <p>
 * Atoms need to be downloaded, either manually or using the method
 * getAtomArrays() in MultipleAlignmentEnsemble.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class MultipleAlignmentXMLParser {

	/**
	 * Creates a list of MultipleAlignment ensembles from an XML file.
	 * This recovers only the information that was previously stored.
	 * If the Atoms are needed, the method getAtomArrays() will automatically
	 * download the structures from the stored structure identifiers.
	 * 
	 * @param xml String XML file containing any number of ensembles
	 * @return List of ensembles in the file
	 * @throws ParserConfigurationException
	 * @throws SAXException
	 * @throws IOException
	 */
	public static List<MultipleAlignmentEnsemble> parseXMLfile(String xml)
			throws ParserConfigurationException, SAXException, IOException {

		List<MultipleAlignmentEnsemble> ensembles = 
				new ArrayList<MultipleAlignmentEnsemble>();

		//Convert string to XML document
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = factory.newDocumentBuilder();
		InputSource inStream = new InputSource();
		inStream.setCharacterStream(new StringReader(xml));
		Document doc = db.parse(inStream);
		doc.getDocumentElement().normalize();

		//In case there are more than one ensemble in the document (generalize)
		NodeList listOfEnsembles = 
				doc.getElementsByTagName("MultipleAlignmentEnsemble");

		//Explore all the ensembles, if multiple ones
		for (int e=0; e<listOfEnsembles.getLength(); e++) {

			Node root = listOfEnsembles.item(e);
			MultipleAlignmentEnsemble ensemble = parseEnsemble(root);
			ensembles.add(ensemble);
		}
		return ensembles;
	}

	public static MultipleAlignmentEnsemble parseEnsemble(Node root){

		MultipleAlignmentEnsemble ensemble = 
				new MultipleAlignmentEnsembleImpl();
		
		parseHeader(root, ensemble);
		
		NodeList children = root.getChildNodes();

		for (int i=0; i<children.getLength(); i++) {

			Node child = children.item(i);
			if (child.getNodeName().equals("MultipleAlignment")){
				parseMultipleAlignment(child, ensemble);
			}
			else if (child.getNodeName().equals("Structures")){
				parseStructures(child, ensemble);
			}
			else if (child.getNodeName().equals("ScoresCache")){
				parseScoresCache(child, ensemble);
			}
		}
		
		return ensemble;
	}

	public static MultipleAlignment parseMultipleAlignment(Node root,
			MultipleAlignmentEnsemble ensemble) {

		MultipleAlignment msa = new MultipleAlignmentImpl(ensemble);
		NodeList children = root.getChildNodes();
		
		for (int i=0; i<children.getLength(); i++) {

			Node child = children.item(i);
			
			if (child.getNodeName().equals("BlockSet")){
				parseBlockSet(child, msa);
			}
			else if (child.getNodeName().equals("ScoresCache")){
				parseScoresCache(child, msa);
			}
		}
		return msa;
	}

	public static BlockSet parseBlockSet(Node root, MultipleAlignment msa) {

		BlockSet bs = new BlockSetImpl(msa);
		List<Matrix4d> transforms = new ArrayList<Matrix4d>();
		NodeList children = root.getChildNodes();
		
		for (int i=0; i<children.getLength(); i++) {

			Node child = children.item(i);
			
			if (child.getNodeName().equals("Block")){
				parseBlock(child, bs);
			}
			else if (child.getNodeName().equals("Matrix4d")){
				Matrix4d t = parseMatrix4d(child);
				transforms.add(t);
			}
			else if (child.getNodeName().equals("ScoresCache")){
				parseScoresCache(child, bs);
			}
		}
		//Because if it is 0 means that there were no transformations
		if (transforms.size() != 0){
			bs.setTransformations(transforms);
		}
		return bs;
	}

	public static Block parseBlock(Node root, BlockSet blockSet) {

		Block b = new BlockImpl(blockSet);
		List<List<Integer>> alignRes = new ArrayList<List<Integer>>();
		b.setAlignRes(alignRes);
		NodeList children = root.getChildNodes();

		for(int i=0; i<children.getLength(); i++) {

			Node child = children.item(i);
			if (child.getNodeName().contains("eqr")){
				
				NamedNodeMap atts = child.getAttributes();
	
				int str = 1;
				Node node = atts.getNamedItem("str"+str);
	
				while (node!=null){
	
					if (alignRes.size() < str) {
						alignRes.add(new ArrayList<Integer>());
					}
					
					String residue = node.getTextContent();
					if (residue.equals("null")){
						alignRes.get(str-1).add(null);
					} else {
						alignRes.get(str-1).add(new Integer(residue));
					}
					
					str++;
					node = atts.getNamedItem("str"+str);
				}
			}
			else if (child.getNodeName().equals("ScoresCache")){
				parseScoresCache(child, b);
			}
		}
		return b;
	}

	public static Matrix4d parseMatrix4d(Node node) {

		Matrix4d m = new Matrix4d();
		NamedNodeMap atts = node.getAttributes();

		for (int x=0; x<4; x++){
			for (int y=0; y<4; y++){
				String key = "mat"+(x+1)+(y+1);
				String value = atts.getNamedItem(key).getTextContent();
				m.setElement(x, y, new Double(value));
			}
		}
		return m;
	}

	public static void parseScoresCache(Node root, ScoresCache cache) {

		NodeList children = root.getChildNodes();

		for (int i=0; i<children.getLength(); i++) {

			Node child = children.item(i);
			NamedNodeMap atts = child.getAttributes();
			if (atts != null) {
				Node score = atts.getNamedItem("value");
				Double value = new Double(score.getTextContent());
				cache.putScore(child.getNodeName(), value);
			}
		}
	}
	
	public static void parseHeader(Node node, 
			MultipleAlignmentEnsemble ensemble) {
		
		NamedNodeMap atts = node.getAttributes();

		String algo = atts.getNamedItem("Algorithm").getTextContent();
		if (!algo.equals("null")){
			ensemble.setAlgorithmName(algo);
		}

		String version = atts.getNamedItem("Version").getTextContent();
		if (!version.equals("null")){
			ensemble.setVersion(version);
		}

		String ioTime = atts.getNamedItem("IOTime").getTextContent();
		if (!ioTime.equals("null")){
			ensemble.setIoTime(new Long(ioTime));
		}

		String time = atts.getNamedItem("CalculationTime").getTextContent();
		if (!time.equals("null")){
			ensemble.setCalculationTime(new Long(time));
		}
	}
	
	public static void parseStructures(Node root, 
			MultipleAlignmentEnsemble ensemble) {
		
		List<String> names = new ArrayList<String>();
		ensemble.setStructureNames(names);

		NamedNodeMap atts = root.getAttributes();

		int str = 1;
		Node node = atts.getNamedItem("name"+str);

		while (node!=null){

			String name = node.getTextContent();
			names.add(name);

			str++;
			node = atts.getNamedItem("name"+str);
		}
	}

}