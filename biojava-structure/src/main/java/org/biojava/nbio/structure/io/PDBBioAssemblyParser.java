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
package org.biojava.nbio.structure.io;

import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;

/** 
 * Parses REMARK 350 records in a PDB file and creates transformations to 
 * construct the quaternary structure of a protein from an asymmetric unit
 * 
 * @author Peter Rose
 * @author Andreas Prlic
 *
 */
public class PDBBioAssemblyParser {	
	
	private static final Logger logger = LoggerFactory.getLogger(PDBBioAssemblyParser.class);
	
	private Integer currentBioMolecule = null;
	private List<String> currentChainIDs = new ArrayList<String>();
	private Matrix currentMatrix = null;
	private double[] shift = null;
	private Map<Integer,BioAssemblyInfo> transformationMap = new HashMap<Integer, BioAssemblyInfo>();
	private int modelNumber = 1;
	private int currentMmSize;
	
	private List<BiologicalAssemblyTransformation> transformations;
	
	/**
	 * Parses REMARK 350 line. See format description:
	 * http://www.wwpdb.org/documentation/format33/remarks2.html
	 * 
	 * @param line
	 */
	public void pdb_REMARK_350_Handler(String line) {
				
		if (line.startsWith("REMARK 350 BIOMOLECULE:")) {
		    initialize();
			currentBioMolecule = Integer.parseInt(line.substring(24).trim());
			
		} else if ( line.matches("REMARK 350 \\w+ DETERMINED BIOLOGICAL UNIT:.*" ) ||
					line.matches("REMARK 350 \\w+ DETERMINED QUATERNARY STRUCTURE:.*" )) {
			// text can be : 
			// author determined biological unit
			// software determined quaternary structure
			currentMmSize = getMmSize(line);
		} else if ( line.startsWith("REMARK 350 APPLY THE FOLLOWING TO CHAINS:")) {
			currentChainIDs.clear();
			addToCurrentChainList(line);	
			
		} else if ( line.startsWith("REMARK 350 IN ADDITION APPLY THE FOLLOWING TO CHAINS:")) {
			currentChainIDs.clear();
			addToCurrentChainList(line);
	
		} else if ( line.startsWith("REMARK 350") && line.contains("AND CHAINS:")) {
			addToCurrentChainList(line);
			
		} else if ( line.startsWith("REMARK 350   BIOMT")) {
	        if (readMatrix(line)) {
	        	saveMatrix();
	        	modelNumber++;
	        }
		}	
	}
	
	/**
	 * Returns a map of bioassembly transformations
	 * @return
	 */
	public Map<Integer, BioAssemblyInfo> getTransformationMap() {
		return transformationMap;
	}
	
	/**
	 * Parses a row of a BIOMT matrix in a REMARK 350 record.
	 * Example: REMARK 350   BIOMT1   2  1.000000  0.000000  0.000000        0.00000 
	 * @param line
	 * @return true if 3rd line of matrix has been parsed (matrix is complete)
	 */
	private boolean readMatrix(String line) {
		// split by one or more spaces
		String[] items = line.split("[ ]+");
			
		// parse BIOMTx, where x is the position in the matrix
		String pos = items[2].substring(5);
		int row = Integer.parseInt(pos);
		if (row == 1) {
			currentMatrix = Matrix.identity(3,3);		
			shift = new double[3];
		}
		
		currentMatrix.set((row-1), 0,Float.parseFloat(items[4]));
		currentMatrix.set((row-1), 1,Float.parseFloat(items[5]));
		currentMatrix.set((row-1), 2,Float.parseFloat(items[6]));
		shift[row-1] = Float.parseFloat(items[7]);

		// return true if 3rd row of matrix has been processed
		return row == 3;
	}
	
	/**
	 * Saves transformation matrix for the list of current chains
	 */
	private void saveMatrix() {

		for (String chainId : currentChainIDs) {
			BiologicalAssemblyTransformation transformation = new BiologicalAssemblyTransformation();
			transformation.setRotationMatrix(currentMatrix.getArray());
			transformation.setTranslation(shift);
			transformation.setId(String.valueOf(modelNumber));
			transformation.setChainId(chainId);
			transformations.add(transformation);
		}
		
		if (!transformationMap.containsKey(currentBioMolecule)) {
			BioAssemblyInfo bioAssembly = new BioAssemblyInfo();
			bioAssembly.setId(currentBioMolecule);
			if (currentMmSize==0) {
				logger.warn("No macromolecular size could be parsed for biological assembly {}",currentBioMolecule);
			}
			bioAssembly.setMacromolecularSize(currentMmSize); 
			bioAssembly.setTransforms(transformations);
			transformationMap.put(currentBioMolecule,bioAssembly);
		}
	}
	
	private int getMmSize(String line) {
		int index = line.indexOf(':');
		String mmString = line.substring(index+1,line.length()-1).trim().toLowerCase();
		return getSizefromString(mmString);
	}

	private static int getSizefromString(String oligomer){
		int size=0;
		
		oligomer = oligomer.toLowerCase();
		
		if (oligomer.equals("monomeric")) {
		    size = 1;
		} else if (oligomer.equals("dimeric")) {
		    size = 2;
		} else if (oligomer.equals("trimeric")) {
		    size = 3;
		} else if (oligomer.equals("tetrameric")) {
		    size = 4;
		} else if (oligomer.equals("pentameric")) {
		    size = 5;
		} else if (oligomer.equals("hexameric")) {
		    size = 6;
		} else if (oligomer.equals("heptameric")) {
		    size = 7;
		} else if (oligomer.equals("octameric")) {
		    size = 8;
		} else if (oligomer.equals("nonameric")) {
		    size = 9;
		} else if (oligomer.equals("decameric")) {
		    size = 10;
		} else if (oligomer.equals("undecameric")) {
		    size = 11;
		} else if (oligomer.equals("dodecameric")) {
		    size = 12;
		} else if (oligomer.equals("tridecameric")) {
		    size = 13;
		} else if (oligomer.equals("tetradecameric")) {
		    size = 14;
		} else if (oligomer.equals("pentadecameric")) {
		    size = 15;
		} else if (oligomer.equals("hexadecameric")) {
		    size = 16;
		} else if (oligomer.equals("heptadecameric")) {
		    size = 17;
		} else if (oligomer.equals("octadecameric")) {
		    size = 18;
		} else if (oligomer.equals("nonadecameric")) {
		    size = 19;
		} else if (oligomer.equals("eicosameric")) {
		    size = 20;
		} else if( oligomer.matches("(\\d+).*")) {
		    size = Integer.parseInt((oligomer.replaceAll("(\\d+).*", "$1")));
		}
		return size;
	}
	
	/**
	 * Parses list of chain ids (A, B, C, etc.)
	 */
	private void addToCurrentChainList(String line) {
		int index = line.indexOf(":");
		String chainList = line.substring(index+1).trim();
        // split by spaces or commas
		String[] chainIds = chainList.split("[ ,]+");
		currentChainIDs.addAll(Arrays.asList(chainIds));
	}
	
	private void initialize() {
		transformations = new ArrayList<BiologicalAssemblyTransformation>();	
		currentMatrix = Matrix.identity(3,3);
		currentBioMolecule = null;
		shift = new double[3];
		modelNumber = 1;
		currentMmSize = 0;
	}
}
