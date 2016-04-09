package org.biojava.nbio.structure.io.mmtf;

import java.io.FileNotFoundException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ExperimentalTechnique;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.PDBCrystallographicInfo;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.biojava.nbio.structure.secstruc.DSSPParser;
import org.biojava.nbio.structure.secstruc.SecStrucCalc;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.SpaceGroup;
import org.rcsb.mmtf.utils.CodecUtils;

/**
 * A utils class of functions needed for Biojava to read and write to mmtf.
 * @author Anthony Bradley
 *
 */
public class MmtfUtils {
	/**
	 * Set up the configuration parameters for BioJava.
	 */
	public static AtomCache setUpBioJava() {
		// Set up the atom cache etc
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		FileParsingParameters params = cache.getFileParsingParams();
		params.setCreateAtomBonds(true);
		params.setAlignSeqRes(true);
		params.setParseBioAssembly(true);
		params.setUseInternalChainId(true);
// MOVE INTO BIOJAVA IF NEED BE
//		CustomChemCompProvider cc = new CustomChemCompProvider();
//		ChemCompGroupFactory.setChemCompProvider(cc);
//		cc.checkDoFirstInstall();
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
		return cache;
	}

	/**
	 * Set up the configuration parameters for BioJava. - with an extra URL
	 */
	public static AtomCache setUpBioJava(String extraUrl) {
		// Set up the atom cache etc
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		FileParsingParameters params = cache.getFileParsingParams();
		params.setCreateAtomBonds(true);
		params.setAlignSeqRes(true);
		params.setParseBioAssembly(true);
		params.setUseInternalChainId(true);
		// MOVE INTO BIOJAVA IF NEED BE
//		CustomChemCompProvider cc = new CustomChemCompProvider();
//		ChemCompGroupFactory.setChemCompProvider(cc);
//		cc.checkDoFirstInstall();
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
		return cache;
	}


	/**
	 * This sets all microheterogeneous groups (previously alternate location groups) as separate groups.
	 * @param bioJavaStruct
	 */
	public static void fixMicroheterogenity(Structure bioJavaStruct) {
		// Loop through the models
		for (int i=0; i<bioJavaStruct.nrModels(); i++){
			// Then the chains
			List<Chain> chains = bioJavaStruct.getModel(i);
			for (Chain c : chains) {
				// Build a new list of groups
				List<Group> outGroups = new ArrayList<>();
				for (Group g : c.getAtomGroups()) {
					List<Group> removeList = new ArrayList<>();
					for (Group altLoc : g.getAltLocs()) {	  
						// Check if they are not equal -> microheterogenity
						if(! altLoc.getPDBName().equals(g.getPDBName())) {
							// Now add this group to the main list
							removeList.add(altLoc);
						}
					}
					// Add this group
					outGroups.add(g);
					// Remove any microhet alt locs
					g.getAltLocs().removeAll(removeList);
					// Add these microhet alt locs
					outGroups.addAll(removeList);
				}
				c.setAtomGroups(outGroups);
			}
		}
	}

	/**
	 * Function to get all the atoms in the strucutre as a list.
	 *
	 * @param bioJavaStruct the bio java struct
	 * @return the all atoms
	 */
	public static List<Atom> getAllAtoms(Structure bioJavaStruct) {
		// Get all the atoms
		List<Atom> theseAtoms = new ArrayList<Atom>();
		for (int i=0; i<bioJavaStruct.nrModels(); i++){
			List<Chain> chains = bioJavaStruct.getModel(i);
			for (Chain c : chains) {
				for (Group g : c.getAtomGroups()) {
					for(Atom a: getAtomsForGroup(g)){
						theseAtoms.add(a);					
					}
				}
			}
		}
		return theseAtoms;
	}

	/**
	 * Function to get a list of atoms for a group.
	 *
	 * @param inputGroup the Biojava Group to consider
	 * @return the atoms for the input Biojava Group
	 */
	public static List<Atom> getAtomsForGroup(Group inputGroup) {
		Set<Atom> uniqueAtoms = new HashSet<Atom>();
		List<Atom> theseAtoms = new ArrayList<Atom>();
		for(Atom a: inputGroup.getAtoms()){
			theseAtoms.add(a);
			uniqueAtoms.add(a);
		}
		List<Group> altLocs = inputGroup.getAltLocs();
		for(Group thisG: altLocs){
			for(Atom a: thisG.getAtoms()){
				if(uniqueAtoms.contains(a)){ 
					continue;
				}
				theseAtoms.add(a);
			}
		}
		return theseAtoms;
	}
	
	
	/**
	 * Function to generate the secondary structure for a Biojava structure object.
	 * @param bioJavaStruct the Biojava structure for which it is to be calculate.
	 */
	public static void calculateDsspSecondaryStructure(Structure bioJavaStruct) {
		SecStrucCalc ssp = new SecStrucCalc();
		try{
			ssp.calculate(bioJavaStruct, true);
		}
		catch(StructureException e) {
			try{
				DSSPParser.fetch(bioJavaStruct.getPDBCode(), bioJavaStruct, true); //download from PDB the DSSP result
			}
			catch(FileNotFoundException enew){
			}
			catch(Exception bige){
				System.out.println(bige);
			}
		}

	}
	
	/**
	 * Get the string representation of a space group.
	 * @param spaceGroup the input SpaceGroup object
	 * @return the space group as a string.
	 */
	public static String getSpaceGroupAsString(SpaceGroup spaceGroup) {
		if(spaceGroup==null){
			return "NA";
		}
		else{
			return spaceGroup.getShortSymbol();
		}
	}

	/**
	 * Get the length six array of the unit cell information.
	 * @param xtalInfo the input PDBCrystallographicInfo object
	 * @return the length six float array
	 */
	public static float[] getUnitCellAsArray(PDBCrystallographicInfo xtalInfo) {
		CrystalCell xtalCell = xtalInfo.getCrystalCell();
		if(xtalCell==null){
			return null;
		}else{
			float[] inputUnitCell = new float[6];
			inputUnitCell[0] = (float) xtalCell.getA();
			inputUnitCell[1] = (float) xtalCell.getB();
			inputUnitCell[2] = (float) xtalCell.getC();
			inputUnitCell[3] = (float) xtalCell.getAlpha();
			inputUnitCell[4] = (float) xtalCell.getBeta();
			inputUnitCell[5] = (float) xtalCell.getGamma();
			return inputUnitCell;
		}
	}

	/**
	 * Converts the set of experimental techniques to an array of strings.
	 * @param experimentalTechniques the input set of experimental techniques
	 * @return the array of strings describing the methods used.
	 */
	public static String[] techniquesToStringArray(Set<ExperimentalTechnique> experimentalTechniques) {
		String[] outArray = new String[experimentalTechniques.size()];
		int index = 0;
		for (ExperimentalTechnique experimentalTechnique : experimentalTechniques) {
			outArray[index] = experimentalTechnique.getName();
			index++;
		}
		return outArray;
	}

	/**
	 * Covert a Date object to ISO time format.
	 * @param inputDate The input date object
	 * @return the time in ISO time format
	 */
	public static String dateToIsoString(Date inputDate) {
		//TODO THIS COULD BE A STATIC METHOD IN A UTILS CLASS
		DateFormat dateStringFormat = new SimpleDateFormat("yyyy-MM-dd");
		return dateStringFormat.format(inputDate);
	}
	
	
	/**
	 * Derive a map mapping chain ids to the chain index.
	 * @param bioJavaStruct the input structure for the map
	 * @return the output map
	 */
	public static Map<String, Integer> getChainIdToIndexMap(Structure bioJavaStruct) {
		// First build a map of asymid -> chain index
		Map<String,Integer> chainIdToIndexMapOne = new HashMap<>();
		int chainCounter = 0;
		for (int i=0; i<bioJavaStruct.nrModels(); i++) {
			for (Chain chain : bioJavaStruct.getChains(i)) {
				String idOne = chain.getChainID();
				if (!chainIdToIndexMapOne.containsKey(idOne)) { 
					chainIdToIndexMapOne.put(idOne, chainCounter);
				}
				chainCounter++;
			}
		}
		return chainIdToIndexMapOne;
	}

	/**
	 * Convert a bioassembly information into a map of transform -> chainindices it relates to
	 * @param bioassemblyInfo
	 * @param chainIdToIndexMap 
	 * @return
	 */
	public static Map<double[], int[]> getTransformMap(BioAssemblyInfo bioassemblyInfo, Map<String, Integer> chainIdToIndexMap) {
		Map<Matrix4d, List<Integer>> matMap = new HashMap<>();
		List<BiologicalAssemblyTransformation> transforms = bioassemblyInfo.getTransforms();
		for (BiologicalAssemblyTransformation transformation : transforms) {
			//			double[] tranMatrix = convertToDoubleArray();
			Matrix4d transMatrix = transformation.getTransformationMatrix();
			int chainIndex = chainIdToIndexMap.get(transformation.getChainId());
			if(matMap.containsKey(transMatrix)){
				matMap.get(transMatrix).add(chainIndex);
			}
			else{
				List<Integer> chainIdList = new ArrayList<>();
				chainIdList.add(chainIndex);
				matMap.put(transMatrix, chainIdList);
			}
		}
		Map<double[], int[]> outMap = new HashMap<>();
		for (Entry<Matrix4d, List<Integer>> entry : matMap.entrySet()) {
			outMap.put(convertToDoubleArray(entry.getKey()), CodecUtils.convertToIntArray(entry.getValue()));
		}
		return outMap;
	}

	/**
	 * Convert a four-d matrix to a double array. Row-packed.
	 * @param transformationMatrix the input matrix4d object
	 * @return the double array (16 long).
	 */
	private static double[] convertToDoubleArray(Matrix4d transformationMatrix) {
		// Initialise the output array
		double[] outArray = new double[16];
		// Iterate over the matrix
		for(int i=0; i<4; i++){
			for(int j=0; j<4; j++){
				// Now set this element
				outArray[i*4+j] = transformationMatrix.getElement(i,j);
			}
		}
		return outArray;
	}
	
	/**
	 * Get a list of all the chains in a structure.
	 * @param structure the input structure
	 * @return the list of chains
	 */
	public static List<Chain> getAllChains(Structure structure) {
		List<Chain> chainList = new ArrayList<>();
		for (int i=0; i<structure.nrModels(); i++) {
			chainList.addAll(structure.getChains(i));
		}
		return chainList;
	}

	/**
	 * Count the total number of groups in the structure
	 * @param structure the input structure
	 * @return the total number of groups
	 */
	public static int getNumGroups(Structure structure) {
		int count = 0;
		for(int i=0; i<structure.nrModels(); i++) {
			for(Chain chain : structure.getChains(i)){
				count+= chain.getAtomGroups().size();
			}
		}
		return count;
	}	
}
