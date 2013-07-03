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
 * Created on 2013-06-13
 * Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.bio.structure.rcsb;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * Fetches information from <a href="http://www.pdb.org/pdb/software/rest.do#descPDB">RCSB's RESTful Web Service Interface</a>.
 * A factory for {@link RCSBLigands RCSBLigands} from {@code ligandInfo} XML files. The factory methods will
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

public class RCSBLigandsFactory {

	private static final String URL_STUB = "http://www.pdb.org/pdb/rest/ligandInfo?structureId=";

	/**
	 * @return An {@link RCSBLigands} from the XML file loaded as {@code stream}. Prefer calling {@link #get(String)} if
	 *         you want data directly from RCSB's RESTful service.
	 * @see RCSBDescriptionFactory#get(String)
	 */
	public static RCSBLigands get(InputStream stream) {

		NodeList data;
		try {
			data = ReadUtils.getNodes(stream);
		} catch (IOException e) {
			ReadUtils.printError(e);
			return null;
		}
		
		// first get the ligandInfo
		RCSBLigands ligands = new RCSBLigands();
		Element structureIdE = null;
		for (int i = 0; i < data.getLength(); i++) {
			if (data.item(i).getNodeType() != 1) continue;
			structureIdE = (Element) data.item(i);
			if (structureIdE.getNodeName().equals("ligandInfo")) {
				break;
			}
		}

		// now get individual ligands
		data = structureIdE.getChildNodes();
		Element ligandE = null;
		for (int i = 0; i < data.getLength(); i++) {
			if (data.item(i).getNodeType() != 1) continue;
			ligandE = (Element) data.item(i);
			if (ligandE.getNodeName().equals("ligand")) {
				if (ligands.getPdbId() == null) {
					ligands.setPdbId(ligandE.getAttribute("structureId"));
				}
				RCSBLigand ligand = makeLigand(ligandE);
				ligands.addLigand(ligand);
			}
		}

		return ligands;

	}

	private static RCSBLigand makeLigand(Element ligandE) {
		RCSBLigand ligand = new RCSBLigand();
		ligand.setId(ligandE.getAttribute("chemicalID"));
		ligand.setType(ligandE.getAttribute("type"));
		ligand.setWeight(ReadUtils.toDouble(ligandE.getAttribute("molecularWeight")));
		Element element = null;
		NodeList data = ligandE.getChildNodes();
		for (int i = 0; i < data.getLength(); i++) {
			if (data.item(i).getNodeType() != 1) continue;
			element = (Element) data.item(i);
			if (element.getNodeName().equals("chemicalName")) {
				ligand.setName(element.getTextContent());
			} else if (element.getNodeName().equals("formula")) {
				ligand.setFormula(element.getTextContent());
			} else if (element.getNodeName().equals("InChIKey")) {
				ligand.setInChIKey(element.getTextContent());
			} else if (element.getNodeName().equals("InChI")) {
				ligand.setInChI(element.getTextContent());
			} else if (element.getNodeName().equals("smiles")) {
				ligand.setSmiles(element.getTextContent());
			}
		}
		return ligand;
	}

	/**
	 * @return An {@link RCSBLigands} from the XML file at
	 *         {@code "http://www.pdb.org/pdb/rest/describeMol?structureId=pdbId"}. This is the preferred factory
	 *         method, unless a different URL or input source is required.
	 * @see RCSBDescriptionFactory#get(InputStream)
	 */
	public static RCSBLigands get(String pdbId) {

		InputStream is;
		try {
			URL url = new URL(URL_STUB + pdbId);
			is = url.openConnection().getInputStream();
		} catch (IOException e) {
			ReadUtils.printError(e);
			return null;
		}
		return get(is);
	}

}
