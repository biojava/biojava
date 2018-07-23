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
 * Created on 2013-06-13 Created by Douglas Myers-Turnbull
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
import java.util.ArrayList;
import java.util.List;

/**
 * Fetches information from <a href="http://www.pdb.org/pdb/software/rest.do#descPDB">RCSB's RESTful Web Service
 * Interface</a>. A factory for {@link RCSBLigands RCSBLigands} from {@code ligandInfo} XML files. The factory methods
 * will return null if the data was not found (rather than throwing an exception); client code should test for this.
 * This is for consistency: if the factory could not read some part (corresponding to a field in a class in
 * {@code rcsb.descriptions}) of the XML file, either because it was blank or contained an error that could not be
 * safely ignored, that field will simply be null. This holds even for numerical values. On some parse errors, the error
 * will additionally be printed to standard error.
 *
 * Example usage:
 *
 * <pre>
 * RCSBLigands ligands = RCSBLigandsFactory.getFromPdbIds(&quot;1w0p&quot;);
 * List&lt;RCSBLigand&gt; list = ligands.getLigands();
 * System.out.println(list.get(0).getFormula()); // prints &quot;CA 2&quot;
 * System.out.println(list.get(1).getFormula()); // prints &quot;C11 H19 N O9&quot;
 * </pre>
 *
 * @see <a href="http://www.pdb.org/pdb/software/rest.do#descPDB">RCSB RESTful</a>
 *
 * @author dmyerstu
 * @since 3.0.6
 */

public class RCSBLigandsFactory {

	private static final String HET_URL_STUB = "http://www.rcsb.org/pdb/rest/describeHet?chemicalID=";

	private static final Logger logger = LoggerFactory.getLogger(RCSBLigandsFactory.class);

	private static final String PDB_URL_STUB = "http://www.rcsb.org/pdb/rest/ligandInfo?structureId=";

	/**
	 * @return A list of {@link RCSBLigand RCSBLigands} from the XML file loaded as {@code stream}. Prefer calling
	 *         {@link #getFromHeteroAtomId(String)} if you want data directly from RCSB's RESTful service.
	 * @see RCSBDescriptionFactory#get(String)
	 */
	public static RCSBLigand getFromHeteroAtomId(InputStream stream) {
		return getFromHeteroAtomIds(stream).get(0);
	}

	/**
	 * @return An {@link RCSBLigands} from the XML file at
	 *         {@code "http://www.pdb.org/pdb/rest/describeHet?chemicalID=hetid"}. This is the preferred factory method,
	 *         unless a different URL or input source is required.
	 * @see RCSBDescriptionFactory#get(InputStream)
	 */
	public static RCSBLigand getFromHeteroAtomId(String heteroAtomId) {
		return getFromHeteroAtomIds(heteroAtomId).get(0);
	}

	/**
	 * @return A list of {@link RCSBLigand RCSBLigands} from the XML file loaded as {@code stream}. Prefer calling
	 *         {@link #getFromHeteroAtomId(String)} if you want data directly from RCSB's RESTful service.
	 * @see RCSBDescriptionFactory#get(String)
	 */
	public static List<RCSBLigand> getFromHeteroAtomIds(InputStream stream) {

		NodeList data;
		try {
			data = ReadUtils.getNodes(stream);
		} catch (IOException e) {
			logger.warn("Couldn't parse XML", e);
			return null;
		}

		List<RCSBLigand> ligands = new ArrayList<RCSBLigand>();

		// first get the ligandInfo
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
				RCSBLigand ligand = makeLigand(ligandE);
				ligands.add(ligand);
			}
		}

		return ligands;

	}

	/**
	 * @return An {@link RCSBLigands} from the XML file at
	 *         {@code "http://www.pdb.org/pdb/rest/describeHet?chemicalID=hetid"}. This is the preferred factory method,
	 *         unless a different URL or input source is required.
	 * @see RCSBDescriptionFactory#get(InputStream)
	 */
	public static List<RCSBLigand> getFromHeteroAtomIds(List<String> heteroAtomIds) {
		String[] x = new String[heteroAtomIds.size()];
		heteroAtomIds.toArray(x);
		return getFromHeteroAtomIds(x); // somewhat cheating here
	}

	/**
	 * @return An {@link RCSBLigands} from the XML file at
	 *         {@code "http://www.pdb.org/pdb/rest/describeHet?chemicalID=hetid"}. This is the preferred factory method,
	 *         unless a different URL or input source is required.
	 * @see RCSBDescriptionFactory#get(InputStream)
	 */
	public static List<RCSBLigand> getFromHeteroAtomIds(String... heteroAtomIds) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < heteroAtomIds.length; i++) {
			if (i > 0) sb.append(",");
			sb.append(heteroAtomIds[i]);
		}
		InputStream is;
		try {
			URL url = new URL(HET_URL_STUB + sb.toString());
			is = url.openConnection().getInputStream();
		} catch (IOException e) {
			logger.warn("Couldn't open connection", e);
			return null;
		}
		return getFromHeteroAtomIds(is);
	}

	/**
	 * @return An {@link RCSBLigands} from the XML file loaded as {@code stream}. Prefer calling
	 *         {@link #getFromPdbId(String)} if you want data directly from RCSB's RESTful service.
	 * @see RCSBDescriptionFactory#get(String)
	 */
	public static RCSBLigands getFromPdbId(InputStream stream) {

		NodeList data;
		try {
			data = ReadUtils.getNodes(stream);
		} catch (IOException e) {
			logger.warn("Couldn't parse XML", e);
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

	/**
	 * @return An {@link RCSBLigands} from the XML file at
	 *         {@code "http://www.pdb.org/pdb/rest/describeMol?structureId=pdbId"}. This is the preferred factory
	 *         method, unless a different URL or input source is required.
	 * @see RCSBDescriptionFactory#get(InputStream)
	 */
	public static RCSBLigands getFromPdbId(String pdbId) {
		InputStream is;
		try {
			URL url = new URL(PDB_URL_STUB + pdbId);
			is = url.openConnection().getInputStream();
		} catch (IOException e) {
			logger.warn("Couldn't open connection", e);
			return null;
		}
		return getFromPdbId(is);
	}

	/**
	 * @return An {@link RCSBLigands} from the XML file loaded as {@code stream}. Prefer calling
	 *         {@link #getFromPdbId(String)} if you want data directly from RCSB's RESTful service.
	 * @see RCSBDescriptionFactory#get(String)
	 */
	public static List<RCSBLigands> getFromPdbIds(InputStream stream) {

		NodeList dataaa;
		try {
			dataaa = ReadUtils.getNodes(stream);
		} catch (IOException e) {
			logger.warn("Couldn't parse XML", e);
			return null;
		}

		// first we have to handle the element "ligandsInEntry", which is not present if we have only 1 structure

		List<RCSBLigands> ligandsList = new ArrayList<RCSBLigands>();

		Element structureIdE = null;

		for (int k = 0; k < dataaa.getLength(); k++) {

			if (dataaa.item(k).getNodeType() != 1) continue;
			structureIdE = (Element) dataaa.item(k);
			if (structureIdE.getNodeName().equals("structureId")) {

				// now get the ligandInfo
				NodeList data = structureIdE.getChildNodes();
				RCSBLigands ligands = new RCSBLigands();
				Element ligandIdE = null;
				for (int i = 0; i < data.getLength(); i++) {
					if (data.item(i).getNodeType() != 1) continue;
					ligandIdE = (Element) data.item(i);
					if (ligandIdE.getNodeName().equals("ligandInfo")) {
						break;
					}
				}

				// now get individual ligands
				data = ligandIdE.getChildNodes();
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

				ligandsList.add(ligands);

			}
		}

		return ligandsList;

	}

	/**
	 * @return An {@link RCSBLigands} from the XML file at
	 *         {@code "http://www.pdb.org/pdb/rest/describeMol?structureId=pdbId"}. This is the preferred factory
	 *         method, unless a different URL or input source is required.
	 * @see RCSBDescriptionFactory#get(InputStream)
	 */
	public static List<RCSBLigands> getFromPdbIds(List<String> pdbIds) {
		String[] x = new String[pdbIds.size()];
		pdbIds.toArray(x);
		return getFromPdbIds(x);
	}

	/**
	 * @return An {@link RCSBLigands} from the XML file at
	 *         {@code "http://www.pdb.org/pdb/rest/describeMol?structureId=pdbId"}. This is the preferred factory
	 *         method, unless a different URL or input source is required.
	 * @see RCSBDescriptionFactory#get(InputStream)
	 */
	public static RCSBLigands getFromPdbIds(String pdbId) {
		InputStream is;
		try {
			URL url = new URL(PDB_URL_STUB + pdbId);
			is = url.openConnection().getInputStream();
		} catch (IOException e) {
			logger.warn("Couldn't open connection", e);
			return null;
		}
		return getFromPdbId(is);
	}

	/**
	 * @return An {@link RCSBLigands} from the XML file at
	 *         {@code "http://www.pdb.org/pdb/rest/describeMol?structureId=pdbId"}. This is the preferred factory
	 *         method, unless a different URL or input source is required.
	 * @see RCSBDescriptionFactory#get(InputStream)
	 */
	public static List<RCSBLigands> getFromPdbIds(String... pdbIds) {
		InputStream is;
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < pdbIds.length; i++) {
			if (i > 0) sb.append(",");
			sb.append(pdbIds[i]);
		}
		try {
			URL url = new URL(PDB_URL_STUB + sb.toString());
			is = url.openConnection().getInputStream();
		} catch (IOException e) {
			logger.warn("Couldn't open connection", e);
			return null;
		}
		return getFromPdbIds(is);
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

}
