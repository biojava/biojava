/**
 * BioJava development code
 *
 * This code may be freely distributed and modified under the terms of the GNU Lesser General Public Licence. This
 * should be distributed with the code. If you do not have a copy, see:
 *
 * http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual authors. These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims, or to join the biojava-l mailing list, visit the home page
 * at:
 *
 * http://www.biojava.org/
 *
 * Created on 2012-11-20 Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.nbio.structure.rcsb;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

/**
 * Fetches information from <a href="http://www.pdb.org/pdb/software/rest.do#descPDB">RCSB's RESTful Web Service
 * Interface</a>. A factory for {@link RCSBDescription RCSBDescriptions} from {@code describeMol} XML files. The factory
 * methods will return null if the data was not found (rather than throwing an exception); client code should test for
 * this. This is for consistency: if the factory could not read some part (corresponding to a field in a class in
 * {@code rcsb.descriptions}) of the XML file, either because it was blank or contained an error that could not be
 * safely ignored, that field will simply be null. This holds even for numerical values. On some parse errors, the error
 * will additionally be printed to standard error.
 *
 * Example usage:
 *
 * <pre>
 * RCSBDescription description = RCSBDescriptionFactory.get(&quot;1w0p&quot;);
 * RCSBLigand firstLigand = ligands.getLigands().get(0);
 * System.out.println(description.getPdbId()); // prints &quot;1w0p&quot;
 * </pre>
 *
 * @see <a href="http://www.pdb.org/pdb/software/rest.do#descPDB">RCSB RESTful</a>
 *
 *      TODO: Handle queries with more than 1 PDB Id.
 *
 * @author dmyerstu
 * @since 3.0.6
 */
public class RCSBDescriptionFactory {

	private static final Logger logger = LoggerFactory.getLogger(RCSBDescriptionFactory.class);

	private static final String URL_STUB = "http://www.rcsb.org/pdb/rest/describeMol?structureId=";

	/**
	 * @return An {@link RCSBDescription} from the XML file loaded as {@code stream}. Prefer calling
	 *         {@link #get(String)} if you want data directly from RCSB's RESTful service.
	 * @see RCSBDescriptionFactory#get(String)
	 */
	public static RCSBDescription get(InputStream stream) {

		NodeList data;
		try {
			data = ReadUtils.getNodes(stream);
		} catch (IOException e) {
			logger.warn("Couldn't parse XML", e);
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
			logger.warn("Couldn't open connection", e);
			return null;
		}
		return get(is);
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
		polymer.setIndex(ReadUtils.toInt(polymerE.getAttribute("entityNr")));
		polymer.setLength(ReadUtils.toInt(polymerE.getAttribute("length")));
		polymer.setWeight(ReadUtils.toDouble(polymerE.getAttribute("weight")));
		polymer.setType(ReadUtils.toStr(polymerE.getAttribute("type")));

		Element element = null;
		NodeList data = polymerE.getChildNodes();
		for (int i = 0; i < data.getLength(); i++) {
			if (data.item(i).getNodeType() != 1) continue;
			element = (Element) data.item(i);
			if (element.getNodeName().equals("chain")) {
				parseChains(polymer, element.getAttribute("id"));
			} else if (element.getNodeName().equals("Taxonomy")) {
				String name = element.getAttribute("name");
				int id = ReadUtils.toInt(element.getAttribute("id"));
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
				logger.warn("Chain id contained more than one character");
			}
		}
	}

	private static void parseSynonyms(RCSBPolymer polymer, String string) {
		String[] parts = string.split("\\s*,\\s*");
		for (String part : parts) {
			polymer.addSynonym(part);
		}
	}

}
