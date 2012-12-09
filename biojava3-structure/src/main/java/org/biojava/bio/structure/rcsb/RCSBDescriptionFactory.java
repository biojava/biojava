/**
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
 * Created on 2012-11-20
 * Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.bio.structure.rcsb;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

/**
 * Fetches information from <a href="http://www.pdb.org/pdb/software/rest.do#descPDB">RCSB's RESTful Web Service Interface</a>.
 * A factory for {@link RCSBDescription RCSBDescriptions} from {@code describeMol} XML files. The factory methods will
 * return null if the data was not found (rather than throwing an exception); client code should test for this. This is
 * for consistency: if the factory could not read some part (corresponding to a field in a class in
 * {@code rcsb.descriptions}) of the XML file, either because it was blank or contained an error that could not be
 * safely ignored, that field will simply be null. This holds even for numerical values. On some parse errors, the error
 * will additionally be printed to standard error.
 * 
 * Example usage:
 * <pre>
 * RCSBDescription description = RCSBDescriptionFactory.get("1w0p");
 * System.out.println(description.getPdbId()); // prints "1w0p"
 * </pre>
 * 
 * @see <a href="http://www.pdb.org/pdb/software/rest.do#descPDB">RCSB RESTful</a>
 * 
 * @author dmyerstu
 * @since 3.0.6
 */
public class RCSBDescriptionFactory {

	// this IS needed
	private static boolean documentBuilderFactorySet = false;

	private static final String URL_STUB = "http://www.rcsb.org/pdb/rest/describeMol?structureId=";

	/**
	 * @return An {@link RCSBDescription} from the XML file loaded as {@code stream}. Prefer calling {@link #get(String)} if
	 *         you want data directly from RCSB's RESTful service.
	 * @see RCSBDescriptionFactory#get(String)
	 */
	public static RCSBDescription get(InputStream stream) {

		NodeList data;
		try {
			data = getNodes(stream);
		} catch (IOException e) {
			printError(e);
			return null;
		}
		
		// first get the main info
		RCSBDescription description = new RCSBDescription();
		Element structureIdE = null;
		for (int i = 0; i < data.getLength(); i++) {
			if (data.item(i).getNodeType() != 1) continue;
			structureIdE = (Element) data.item(i);
			if (structureIdE.getNodeName().equals("structureId")) {
				description.setPdbId(structureIdE.getAttribute("id"));
			}
		}

		// now get polymers
		data = structureIdE.getChildNodes();
		Element polymerE = null;
		for (int i = 0; i < data.getLength(); i++) {
			if (data.item(i).getNodeType() != 1) continue;
			polymerE = (Element) data.item(i);
			if (polymerE.getNodeName().equals("polymer")) {
				RCSBPolymer polymer = makePolymer(polymerE);
				description.addPolymer(polymer);
			}
		}

		return description;

	}

	/**
	 * @return An {@link RCSBDescription} from the XML file at
	 *         {@code "http://www.pdb.org/pdb/rest/describeMol?structureId=pdbId"}. This is the preferred factory
	 *         method, unless a different URL or input source is required.
	 * @see RCSBDescriptionFactory#get(InputStream)
	 */
	public static RCSBDescription get(String pdbId) {

		InputStream is;
		try {
			URL url = new URL(URL_STUB + pdbId);
			is = url.openConnection().getInputStream();
		} catch (IOException e) {
			printError(e);
			return null;
		}
		return get(is);
	}

	/**
	 * @param stream
	 * @return A {@link NodeList} of top-level {@link Node Nodes} in {@code stream}.
	 * @throws IOException
	 */
	private static NodeList getNodes(InputStream stream) throws IOException {

		if (!documentBuilderFactorySet) { // it's really stupid, but we have to do this
			System.setProperty("javax.xml.parsers.DocumentBuilderFactory",
					"com.sun.org.apache.xerces.internal.jaxp.DocumentBuilderFactoryImpl");
			documentBuilderFactorySet = true;
		}
		DocumentBuilderFactory builderFactory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = null;
		Document document = null;
		try {
			builder = builderFactory.newDocumentBuilder();
		} catch (ParserConfigurationException e) {
			printError(e);
			stream.close();
			throw new IOException(e);
		}
		try {
			document = builder.parse(stream);
		} catch (SAXException e) {
			System.out.println(e.getMessage());
			printError(e);
			stream.close();
			throw new IOException(e);
		}
		Node root = document.getDocumentElement();
		return root.getChildNodes();
	}

	private static RCSBMacromolecule makeMolecule(Element moleculeE) {
		RCSBMacromolecule molecule = new RCSBMacromolecule();
		molecule.setName(moleculeE.getAttribute("name"));
		Element element = null;
		NodeList data = moleculeE.getChildNodes();
		for (int i = 0; i < data.getLength(); i++) {
			if (data.item(i).getNodeType() != 1) continue;
			element = (Element) data.item(i);
			if (element.getNodeName().equals("accession")) {
				molecule.addAccession(element.getAttribute("id"));
			}
		}
		return molecule;
	}

	private static RCSBPolymer makePolymer(Element polymerE) {

		RCSBPolymer polymer = new RCSBPolymer();
		polymer.setIndex(toInt(polymerE.getAttribute("entityNr")));
		polymer.setLength(toInt(polymerE.getAttribute("length")));
		polymer.setWeight(toDouble(polymerE.getAttribute("weight")));
		polymer.setType(toStr(polymerE.getAttribute("type")));

		Element element = null;
		NodeList data = polymerE.getChildNodes();
		for (int i = 0; i < data.getLength(); i++) {
			if (data.item(i).getNodeType() != 1) continue;
			element = (Element) data.item(i);
			if (element.getNodeName().equals("chain")) {
				parseChains(polymer, element.getAttribute("id"));
			} else if (element.getNodeName().equals("Taxonomy")) {
				String name = element.getAttribute("name");
				int id = toInt(element.getAttribute("id"));
				RCSBTaxonomy taxonomy = new RCSBTaxonomy(name, id);
				polymer.setTaxonomy(taxonomy);
			} else if (element.getNodeName().equals("macroMolecule")) {
				RCSBMacromolecule molecule = makeMolecule(element);
				polymer.setMolecule(molecule);
			} else if (element.getNodeName().equals("polymerDescription")) {
				polymer.setDescription(element.getAttribute("description"));
			} else if (element.getNodeName().equals("enzClass")) {
				polymer.setEnzClass(element.getAttribute("ec"));
			} else if (element.getNodeName().equals("synonym")) {
				parseSynonyms(polymer, element.getAttribute("name"));
			}
		}
		return polymer;
	}

	private static void parseChains(RCSBPolymer polymer, String string) {
		String[] parts = string.split("\\s*,\\s*");
		for (String part : parts) {
			if (part.length() == 1) {
				polymer.addChain(part.charAt(0));
			} else {
				printError(new Exception("Chain id contained more than one character."));
			}
		}
	}

	/**
	 * Prints an error message for {@code e} that shows causes and suppressed messages recursively.
	 * Just a little more useful than {@code e.printStackTrace()}.
	 * 
	 * @param e
	 */
	public static void printError(Exception e) {
		System.err.println(printError(e, ""));
	}

	/**
	 * @see #printError(Exception)
	 */
	private static String printError(Exception e, String tabs) {
		StringBuilder sb = new StringBuilder();
		Throwable prime = e;
		while (prime != null) {
			if (tabs.length() > 0) sb.append(tabs + "Cause:" + "\n");
			sb.append(tabs + prime.getClass().getSimpleName());
			if (prime.getMessage() != null) sb.append(": " + prime.getMessage());
			sb.append("\n");
			if (prime instanceof Exception) {
				StackTraceElement[] trace = ((Exception) prime).getStackTrace();
				for (StackTraceElement element : trace) {
					sb.append(tabs + element.toString() + "\n");
				}
			}
			prime = prime.getCause();
			tabs += "\t";
		}
		sb.append("\n");
		return sb.toString();
	}

	private static void parseSynonyms(RCSBPolymer polymer, String string) {
		String[] parts = string.split("\\s*,\\s*");
		for (String part : parts) {
			polymer.addSynonym(part);
		}
	}

	private static Double toDouble(String s) {
		if (s == "") return null;
		try {
			return Double.parseDouble(s);
		} catch (NumberFormatException e) {
			printError(e);
		}
		return null;
	}

	private static Integer toInt(String s) {
		if (s == "") return null;
		try {
			return Integer.parseInt(s);
		} catch (NumberFormatException e) {
			printError(e);
		}
		return null;
	}

	/**
	 * @param s
	 * @return {@code s}, or null if {@code s} is the empty string
	 */
	private static String toStr(String s) {
		if (s == "") return null;
		return s;
	}

}
