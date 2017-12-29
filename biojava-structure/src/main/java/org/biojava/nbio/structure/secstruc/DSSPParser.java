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
package org.biojava.nbio.structure.secstruc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Class to parse a DSSP file (output of the DSSP program),
 * that contains the secondary structure assignment of a structure.
 * <p>
 * This class has been ported from the OWL Java library for
 * Structural Bioinformatics (https://github.com/eppic-team/owl).
 * <p>
 * As of September 2015, the DSSP source code and executables can
 * be downloaded from http://swift.cmbi.ru.nl/gv/dssp/.
 *
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class DSSPParser {

	private static final Logger logger =
			LoggerFactory.getLogger(DSSPParser.class);

	/**
	 * Parse a DSSP output file and return the secondary structure
	 * annotation as a List of {@link SecStrucState} objects.
	 *
	 * @param dsspIs an InputStream to a DSSP file
	 * @param structure Structure object associated to the dssp
	 * @param assign assigns the SS to the structure if true
	 * @return a List of SS annotation objects
	 * @throws StructureException
	 * @throws IOException
	 */
	public static List<SecStrucState> parseInputStream(InputStream dsspIs,
			Structure structure, boolean assign)
					throws IOException, StructureException {
			
		BufferedReader reader = new BufferedReader(new InputStreamReader(dsspIs));
		return generalParse(reader, structure, assign);
	}

	/**
	 * Parse a DSSP output file and return the secondary structure
	 * annotation as a List of {@link SecStrucState} objects.
	 *
	 * @param dsspPath path to the DSSP file to parse
	 * @param structure Structure object associated to the dssp
	 * @param assign assigns the SS to the structure if true
	 * @return a List of SS annotation objects
	 * @throws StructureException
	 * @throws IOException
	 */
	public static List<SecStrucState> parseFile(String dsspPath,
			Structure structure, boolean assign)
					throws IOException, StructureException {

		File file = new File(dsspPath);
		Reader read = new FileReader(file);
		BufferedReader reader = new BufferedReader(read);
		return generalParse(reader, structure, assign);
	}

	/**
	 * Fetch and parse the DSSP file of the specified pdb code
	 * from the PDB web server and return the secondary structure
	 * annotation as a List of {@link SecStrucState} objects.
	 *
	 * @param pdb path to the DSSP file to parse
	 * @param structure Structure object associated to the dssp
	 * @param assign assigns the SS to the structure if true
	 * @return a List of SS annotation objects
	 * @throws StructureException
	 * @throws IOException
	 */
	public static List<SecStrucState> fetch(String pdb,
			Structure structure, boolean assign)
					throws IOException, StructureException {

		URL url = new URL("http://files.rcsb.org/dssp/" + 
				pdb.toLowerCase().substring(1, 3) + "/" + 
				pdb.toLowerCase() + "/" +
				pdb.toLowerCase() + ".dssp.gz");
		InputStream in = new GZIPInputStream(url.openStream());
		Reader read = new InputStreamReader(in);
		BufferedReader reader = new BufferedReader(read);
		return generalParse(reader, structure, assign);
	}

	/**
	 * Parse a DSSP format String and return the secondary structure
	 * annotation as a List of {@link SecStrucState} objects.
	 *
	 * @param dsspOut String with the DSSP output to parse
	 * @param structure Structure object associated to the dssp
	 * @param assign assigns the SS to the structure if true
	 * @return a List of SS annotation objects
	 * @throws StructureException
	 * @throws IOException
	 */
	public static List<SecStrucState> parseString(String dsspOut,
			Structure structure, boolean assign)
					throws IOException, StructureException {

		Reader read = new StringReader(dsspOut);
		BufferedReader reader = new BufferedReader(read);
		return generalParse(reader, structure, assign);
	}

	private static List<SecStrucState> generalParse(BufferedReader reader,
			Structure structure, boolean assign)
					throws IOException, StructureException {

		String startLine = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC";
		String line;

		List<SecStrucState> secstruc = new ArrayList<SecStrucState>();

		//Find the first line of the DSSP output
		while((line = reader.readLine()) != null) {
			if(line.startsWith(startLine)) break;
		}

		while((line = reader.readLine()) != null) {

			String indexStr = line.substring(0,5).trim();
			String resNumStr = line.substring(5,10).trim();

			//Only happens if dssp inserts a line indicating a chain break
			if(!resNumStr.equals("")) {

				int index = Integer.parseInt(indexStr);
				//Get the group of the structure corresponding to the residue
				int resNum = Integer.parseInt(resNumStr);
				char insCode = line.charAt(10);
				String chainId = line.substring(11,13).trim();
				ResidueNumber r = new ResidueNumber(chainId, resNum, insCode);
				Group parent = structure.getPolyChainByPDB(chainId)
						.getGroupByPDB(r);
				SecStrucType ssType =
						SecStrucType.fromCharacter(line.charAt(16));

				SecStrucState ss = new SecStrucState(parent,
						SecStrucInfo.DSSP_ASSIGNMENT, ssType);

				//Parse the Bridge partners - TODO parallel or antiparallel?
				String bp = line.substring(25,29).trim();
				if (bp != "") {
					BetaBridge bb = new BetaBridge(
							index, Integer.valueOf(bp), BridgeType.parallel);
					ss.addBridge(bb);
				} else logger.warn("Unable to parse beta Bridge for resn "+index);

				bp = line.substring(29,33).trim();
				if (bp != "") {
					BetaBridge bb = new BetaBridge(
							index, Integer.valueOf(bp), BridgeType.parallel);
					ss.addBridge(bb);
				} else logger.warn("Unable to parse beta Bridge for resn "+index);

				//Parse the energy terms of donor and acceptor
				for (int i=0; i<4; i++){

					int a = 42 + i*11;
					int b = a + 8;

					String val = line.substring(a,b).trim();
					if (val == "") {
						logger.warn("Unable to parse energy for resn "+index);
						continue;
					}

					String[] p = val.split(",");

					int partner = Integer.parseInt(p[0]);
					if (partner != 0) partner += index;
					double energy = Double.valueOf(p[1]) * 1000.0;

					switch(i){
					case 0:
						ss.getAccept1().setPartner(partner);
						ss.getAccept1().setEnergy(energy);
						break;
					case 1:
						ss.getDonor1().setPartner(partner);
						ss.getDonor1().setEnergy(energy);
						break;
					case 2:
						ss.getAccept2().setPartner(partner);
						ss.getAccept2().setEnergy(energy);
						break;
					case 3:
						ss.getDonor2().setPartner(partner);
						ss.getDonor1().setEnergy(energy);
						break;
					}
				}

				//Angle properties
				String val = line.substring(91,97).trim();
				if (val != "") ss.setKappa(Float.valueOf(val));
				else logger.warn("Unable to parse kappa for resn "+index);

				val = line.substring(103,109).trim();
				if (val != "") ss.setPhi(Float.valueOf(val));
				else logger.warn("Unable to parse phi for resn "+index);

				val = line.substring(109,116).trim();
				if (val != "") ss.setPsi(Float.valueOf(val));
				else logger.warn("Unable to parse psi for resn "+index);

				if (assign) parent.setProperty(Group.SEC_STRUC, ss);
				secstruc.add(ss);
			}
		}

		reader.close();
		return secstruc;
	}

}
