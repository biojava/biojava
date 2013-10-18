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
package org.biojava.bio.structure.io;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.jama.Matrix;
import org.biojava.bio.structure.quaternary.ModelTransformationMatrix;

/** 
 * Parses REMARK 350 records in a PDB file and creates transformations to 
 * construct the quaternary structure of a protein from an asymmetric unit
 * 
 * @author Peter Rose
 * @author Andreas Prlic
 *
 */
public class PDBBioAssemblyParser {	
	private Integer currentBioMolecule = null;
	private List<String> currentChainIDs = new ArrayList<String>();
	private Matrix currentMatrix = null;
	private double[] shift = null;
	private Map<Integer,List<ModelTransformationMatrix>> transformationMap = new HashMap<Integer, List<ModelTransformationMatrix>>();
	private int modelNumber = 1;
	
	private List<ModelTransformationMatrix> transformations;
	
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
	public Map<Integer, List<ModelTransformationMatrix>> getTransformationMap() {
		return transformationMap;
	}
	
	/**
	 * Returns the number of bioassemblies
	 * @return number of bioassemblies
	 */
	public int getNrBioAssemblies(){
		return transformationMap.size();
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
		
		// note, BioJava uses a transposed form of the rotation matrix
		currentMatrix.set(0,(row-1),Float.parseFloat(items[4]));
		currentMatrix.set(1,(row-1),Float.parseFloat(items[5]));
		currentMatrix.set(2,(row-1),Float.parseFloat(items[6]));	
		shift[row-1] = Float.parseFloat(items[7]);

		// return true if 3rd row of matrix has been processed
		return row == 3;
	}
	
	/**
	 * Saves transformation matrix for the list of current chains
	 */
	private void saveMatrix() {
		for (String chainId : currentChainIDs) {
			ModelTransformationMatrix transformation = new ModelTransformationMatrix();
			transformation.setMatrix(currentMatrix);
			transformation.setVector(shift);
			transformation.setId(String.valueOf(modelNumber));
			transformation.setChainId(chainId);
			transformations.add(transformation);
		}
			
		transformationMap.put(currentBioMolecule,transformations);
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
		System.out.println("chainids" + currentChainIDs);
	}
	
	private void initialize() {
		transformations = new ArrayList<ModelTransformationMatrix>();	
		currentMatrix = Matrix.identity(3,3);
		currentBioMolecule = null;
		shift = new double[3];
		modelNumber = 1;
	}
}
