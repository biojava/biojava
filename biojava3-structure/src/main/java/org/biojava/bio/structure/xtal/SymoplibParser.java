package org.biojava.bio.structure.xtal;

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
	
	// the symop library from the CCP4 package
	private static final String SYMOPFILE = "/org/biojava/bio/structure/xtal/symop.lib";
	
	private static final Pattern namePat = Pattern.compile(".*\\s([A-Z]+)(\\s'.+')?\\s+'(.+)'.*");
	
	private static final InputStream symoplibIS = SymoplibParser.class.getResourceAsStream(SYMOPFILE);
	private static final TreeMap<Integer, SpaceGroup> sgs = parseSymopLib();
	
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
	 * Gets the space group for the given international short name, using
	 * the PDB format, e.g. 'P 21 21 21' or 'C 1 c 1'
	 * @param shortName
	 * @return the SpaceGroup or null if the shortName is not valid
	 */
	public static SpaceGroup getSpaceGroup(String shortName) {
		if (shortName==null) return null;
		
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
	
	private static TreeMap<Integer,SpaceGroup> parseSymopLib() {
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
			System.err.println("Fatal error! Can't read resource file "+SYMOPFILE+". Error: "+e.getMessage()+". Exiting.");
			System.exit(1);
		}

		for (SpaceGroup sg:map.values()) {
			sg.initializeCellTranslations();
		}
		return map;
	}
	
	

}
