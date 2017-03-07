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
//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;

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

	//private static final Logger logger = LoggerFactory.getLogger(PDBBioAssemblyParser.class);

	private Integer currentBioMolecule = null;
	private List<String> currentChainIDs = new ArrayList<String>();
	private Matrix currentMatrix = null;
	private double[] shift = null;
	private Map<Integer,BioAssemblyInfo> transformationMap = new HashMap<Integer, BioAssemblyInfo>();
	private int modelNumber = 1;

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

		} 
		// not parsing anymore the size (from biojava 5.0), thus this is not needed anymore
		// eventually if needed this could be used to
		// infer if bioassembly is author or software determined
		//else if ( line.matches("REMARK 350 \\w+ DETERMINED BIOLOGICAL UNIT:.*" ) ||
		//			line.matches("REMARK 350 \\w+ DETERMINED QUATERNARY STRUCTURE:.*" )) {
			// text can be :
			// author determined biological unit
			// software determined quaternary structure
		//} 
		else if ( line.startsWith("REMARK 350 APPLY THE FOLLOWING TO CHAINS:")) {
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
			bioAssembly.setTransforms(transformations);
			transformationMap.put(currentBioMolecule,bioAssembly);
		}
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
	}
	
	/**
	 * Set the macromolecularSize fields of the parsed bioassemblies.
	 * This can only be called after the full PDB file has been read so that
	 * all the info for all bioassemblies has been gathered.
	 * Note that an explicit method to set the field is necessary here because 
	 * in PDB files the transformations contain only the author chain ids, corresponding 
	 * to polymeric chains, whilst in mmCIF files the transformations 
	 * contain all asym ids of both polymers and non-polymers.
	 */
	public void setMacromolecularSizes() {
		for (BioAssemblyInfo bioAssembly : transformationMap.values()) {
			bioAssembly.setMacromolecularSize(bioAssembly.getTransforms().size());
		}
	}
}
