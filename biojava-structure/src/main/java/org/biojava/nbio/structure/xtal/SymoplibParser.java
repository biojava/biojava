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
package org.biojava.nbio.structure.xtal;

import org.biojava.nbio.structure.xtal.io.SpaceGroupMapRoot;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.xml.bind.JAXBException;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 * A class containing static methods to parse the symop.lib file from the
 * CCP4 package. The file contains the transformations belonging to all
 * protein crystallography space groups.
 *
 * See http://structure.usc.edu/ccp4/symlib.html for documentation
 *
 * @author duarte_j
 *
 */
public class SymoplibParser {

	private static final Logger logger = LoggerFactory.getLogger(SymoplibParser.class);

	private static final String newline = System.getProperty("line.separator");

	private static final String SPACE_GROUPS_FILE = "org/biojava/nbio/structure/xtal/spacegroups.xml";

	private static final Pattern namePat = Pattern.compile(".*\\s([A-Z]+)(\\s'.+')?\\s+'(.+)'.*");

	private static  TreeMap<Integer, SpaceGroup> sgs = parseSpaceGroupsXML();


	private static HashMap<String, SpaceGroup> name2sgs; // map for lookups based on short names

	/**
	 * Gets the space group for the given standard identifier.
	 * See for example http://en.wikipedia.org/wiki/Space_group
	 * @param id
	 * @return
	 */
	public static SpaceGroup getSpaceGroup(int id) {
		return sgs.get(id);
	}


	/**
	 * Load all SpaceGroup information from the file spacegroups.xml
	 *
	 * @return a map providing information for all spacegroups
	 */
	private static TreeMap<Integer, SpaceGroup> parseSpaceGroupsXML() {

		// NOTE: if the space group file is requested by some part of the code (i.e. this method is called) and
		//       there is a problem in reading it, then that's truly a FATAL problem, since this is not a user file
		//       but a file that's part of the distribution: it MUST be there and MUST have the right format. A failure
		//       to read it is more of a "compilation" error than a runtime error. That's the reason that System.exit
		//       is called (which otherwise usually is not a good idea).
		//
		//       The rest of the application will simply not work: there are 3 options to handle it
		//	     a) returning null and then a NullPointer will happen down the line and thus a not very clear
		//          error message will be printed
		//       b) throw the exception forward and catch it in the final main but that would also be bad because
		//          this is a file that the user didn't input but that should be part of the distribution
		//		 c) call System.exit(1) and "crash" the application with a human-understandable error message

		InputStream spaceGroupIS = SymoplibParser.class.getClassLoader().getResourceAsStream(SPACE_GROUPS_FILE);

		if ( spaceGroupIS == null) {
			logger.error("Fatal error! Could not find resource: " + SPACE_GROUPS_FILE + ". This probably means that your biojava jar file is corrupt or incorrectly built.");
			System.exit(1);
		}

		TreeMap<Integer, SpaceGroup> map = new TreeMap<Integer, SpaceGroup>();

		try {
			map = parseSpaceGroupsXML(spaceGroupIS);
		} catch (IOException e) {
			logger.error("Fatal error! Could not parse resource: "+SPACE_GROUPS_FILE+". Error: "+e.getMessage());
			System.exit(1);
		} catch (JAXBException e) {
			logger.error("Fatal error! Could not parse resource: "+SPACE_GROUPS_FILE+". Problem in xml formatting: "+e.getMessage());
			System.exit(1);
		}

		name2sgs = new HashMap<String, SpaceGroup>();

		for (SpaceGroup sg:map.values()) {

			sg.initializeCellTranslations();
			name2sgs.put(sg.getShortSymbol(), sg);
			if (sg.getAltShortSymbol()!=null) {
				// we add also alternative name to map so we can look it up
				name2sgs.put(sg.getAltShortSymbol(), sg);
			}
		}

		return map;

	}


	/**
	 * Load all SpaceGroup information from the file spacegroups.xml
	 *
	 * @return a map providing information for all spacegroups
	 */
	public static TreeMap<Integer, SpaceGroup> parseSpaceGroupsXML(
			InputStream spaceGroupIS) throws IOException, JAXBException {

		String xml = convertStreamToString(spaceGroupIS);

		SpaceGroupMapRoot spaceGroups = SpaceGroupMapRoot.fromXML(xml);
		return spaceGroups.getMapProperty();

	}


	private static String convertStreamToString(InputStream stream) throws IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
		StringBuilder sb = new StringBuilder();

		String line = null;

		while ((line = reader.readLine()) != null) {
			sb.append(line).append(newline);
		}

		return sb.toString();
	}

	/**
	 * Get the space group for the given international short name, using
	 * the PDB format, e.g. 'P 21 21 21' or 'C 1 c 1'
	 * @param shortName
	 * @return the SpaceGroup or null if the shortName is not valid
	 */
	public static SpaceGroup getSpaceGroup(String shortName) {
		if (shortName==null || shortName.length()<=2) return null;

		// PDB uses group "P 1-" for 13 racemic mixture entries (as of Sep2011), e.g. 3e7r
		// they call the space group "P 1-" unusually (symop.lib and everyone else call it "P -1")
		if (shortName.equals("P 1-")) shortName="P -1";

		// enantiomorphic space groups contain sometime letters indicating glide planes which should always be lower case
		// in some PDB entries like 4gwv they are in upper case, we fix that here: convert any non-first letter to lower case
		shortName = shortName.substring(0, 1)+shortName.substring(1).toLowerCase();

		return name2sgs.get(shortName);
	}

	public static TreeMap<Integer,SpaceGroup> getAllSpaceGroups() {
		return sgs;
	}


	/**
	 * A parser for the symop.lib file provided by CCP4. Note: this file is not getting re-distributed by BioJava.
	 * It can be downloaded from:
	 *
	 *  http://www.ccp4.ac.uk/cvs/viewvc.cgi/libccp4/data/symop.lib?revision=1.10&view=markup
	 *
	 * Note: this file is not needed by BioJava. BioJava loads equivalent information from the file spacegroups.xml
	 *
	 * @param symoplibIS
	 * @return
	 */
	public static TreeMap<Integer,SpaceGroup> parseSymopLib(InputStream symoplibIS) {
		TreeMap<Integer, SpaceGroup> map = new TreeMap<Integer, SpaceGroup>();
		name2sgs = new HashMap<String, SpaceGroup>();
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(symoplibIS));
			String line;
			SpaceGroup currentSG = null;
			while ((line=br.readLine())!=null) {
				if (!line.startsWith(" ")) {
					if (currentSG!=null) {
						map.put(currentSG.getId(),currentSG);
						name2sgs.put(currentSG.getShortSymbol(), currentSG);
						if (currentSG.getAltShortSymbol()!=null) {
							// we add also alternative name to map so we can look it up
							name2sgs.put(currentSG.getAltShortSymbol(), currentSG);
						}
					}

					int idxFirstSpace = line.indexOf(' ');
					int idxSecondSpace = line.indexOf(' ',idxFirstSpace+1);
					int idxThirdSpace = line.indexOf(' ',idxSecondSpace+1);
					int id = Integer.parseInt(line.substring(0, idxFirstSpace));
					int multiplicity = Integer.parseInt(line.substring(idxFirstSpace+1, idxSecondSpace));
					int primitiveMultiplicity = Integer.parseInt(line.substring(idxSecondSpace+1, idxThirdSpace));
					Matcher m = namePat.matcher(line);
					String shortSymbol = null;
					String altShortSymbol = null;
					String brav = null;
					if (m.matches()) {
						brav = m.group(1);
						altShortSymbol = m.group(2); // null if there is no match
						if (altShortSymbol!=null) altShortSymbol = altShortSymbol.trim().replaceAll("'", "");
						shortSymbol = m.group(3);
					}
					currentSG = new SpaceGroup(id, multiplicity, primitiveMultiplicity, shortSymbol, altShortSymbol, BravaisLattice.getByName(brav));
				} else {
					currentSG.addTransformation(line.trim());
				}
			}
			br.close();
			// and we add the last SG
			map.put(currentSG.getId(), currentSG);
			name2sgs.put(currentSG.getShortSymbol(), currentSG);
			if (currentSG.getAltShortSymbol()!=null) {
				// we add also alternative name to map so we can look it up
				name2sgs.put(currentSG.getAltShortSymbol(), currentSG);
			}

		} catch (IOException e) {
			logger.error("Fatal error! Can't read symop.lib file. Error: "+e.getMessage()+". ");
			System.exit(1);
		}

		for (SpaceGroup sg:map.values()) {
			sg.initializeCellTranslations();
		}
		return map;
	}


}
