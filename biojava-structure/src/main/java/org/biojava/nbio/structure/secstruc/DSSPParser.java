package org.biojava.nbio.structure.secstruc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;

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
	 * @param dsspPath path to the DSSP file to parse
	 * @param structure Structure object associated to the dssp
	 * @param assign assigns the SS to the structure if true
	 * @return a List of SS annotation objects
	 * @throws StructureException 
	 * 
	 */
	public static List<SecStrucState> parseDSSP(String dsspPath, 
			Structure structure, boolean assign) 
					throws IOException, StructureException {

		String startLine = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC";
		String line;
		
		File file = new File(dsspPath);
		Reader reader = new FileReader(file);
		BufferedReader dsspOutput = new BufferedReader(reader);
		
		List<SecStrucState> secstruc = new ArrayList<SecStrucState>();
		
		//Find the first line of the DSSP output
		while((line = dsspOutput.readLine()) != null) {
			if(line.startsWith(startLine)) break;
		}
		
		while((line = dsspOutput.readLine()) != null) {
			
			String indexStr = line.substring(0,5).trim();
			String resNumStr = line.substring(5,10).trim();
			
			//Only happens if dssp inserts a line indicating a chain break
			if(!resNumStr.equals("")) {	
				
				int index = Integer.valueOf(indexStr);
				//Get the group of the structure corresponding to the residue
				int resNum = Integer.valueOf(resNumStr);
				char insCode = line.charAt(10);
				String chainId = line.charAt(11)+"";
				ResidueNumber r = new ResidueNumber(chainId, resNum, insCode);
				Group parent = structure.getChainByPDB(chainId)
						.getGroupByPDB(r);
				SecStrucType ssType = 
						SecStrucType.fromCharacter(line.charAt(16));
				
				SecStrucState ss = new SecStrucState(parent, 
						SecStrucInfo.DSSP_FILE_ASSIGNMENT, ssType);
				
				//TODO ignore for now
				int BP1 = Integer.valueOf(line.charAt(28));
				int BP2 = Integer.valueOf(line.charAt(32));
				
				//Parse the energy terms of donor and acceptor
				for (int i=0; i<4; i++){
					
					int a = 42 + i*11;
					int b = a + 8;
					
					String val = line.substring(a,b).trim();
					if (val == "") continue;
					
					String[] p = val.split(",");
					
					int partner = Integer.valueOf(p[0]);
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
				
				val = line.substring(103,109).trim();
				if (val != "") ss.setPhi(Float.valueOf(val));
				
				val = line.substring(109,116).trim();
				if (val != "") ss.setPsi(Float.valueOf(val));
				
				if (assign) parent.setProperty(Group.SEC_STRUC, ss);
				secstruc.add(ss);
			}
		}
		
		reader.close();
		dsspOutput.close();
		return secstruc;
	}

}