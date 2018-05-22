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
package org.biojava.nbio.structure.io.mmtf;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Bond;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ExperimentalTechnique;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.NucleotideImpl;
import org.biojava.nbio.structure.PDBCrystallographicInfo;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.biojava.nbio.structure.secstruc.DSSPParser;
import org.biojava.nbio.structure.secstruc.SecStrucCalc;
import org.biojava.nbio.structure.secstruc.SecStrucState;
import org.biojava.nbio.structure.secstruc.SecStrucType;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.SpaceGroup;
import org.rcsb.mmtf.dataholders.DsspType;
import org.rcsb.mmtf.utils.CodecUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A utils class of functions needed for Biojava to read and write to mmtf.
 * @author Anthony Bradley
 *
 */
public class MmtfUtils {
	
	private static final Logger LOGGER = LoggerFactory.getLogger(MmtfUtils.class);
	
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
		DownloadChemCompProvider cc = new DownloadChemCompProvider();
		ChemCompGroupFactory.setChemCompProvider(cc);
		cc.checkDoFirstInstall();
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
		return cache;
	}

	/**
	 * Set up the configuration parameters for BioJava.
	 * @param extraUrl the string describing the URL (or file path) from which
	 * to get missing CCD entries.
	 */
	public static AtomCache setUpBioJava(String extraUrl) {
		// Set up the atom cache etc
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		FileParsingParameters params = cache.getFileParsingParams();
		params.setCreateAtomBonds(true);
		params.setAlignSeqRes(true);
		params.setParseBioAssembly(true);
		DownloadChemCompProvider.serverBaseUrl = extraUrl;
		DownloadChemCompProvider.useDefaultUrlLayout = false;
		DownloadChemCompProvider cc = new DownloadChemCompProvider();
		ChemCompGroupFactory.setChemCompProvider(cc);
		cc.checkDoFirstInstall();
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
		return cache;
	}


	/**
	 * This sets all microheterogeneous groups 
	 * (previously alternate location groups) as separate groups.
	 * This is required because mmtf groups cannot have multiple HET codes.
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
	 * Generate the secondary structure for a Biojava structure object.
	 * @param bioJavaStruct the Biojava structure for which it is to be calculate.
	 */
	public static void calculateDsspSecondaryStructure(Structure bioJavaStruct) {
		SecStrucCalc ssp = new SecStrucCalc();

		try{
			ssp.calculate(bioJavaStruct, true);
		}
		catch(StructureException e) {
			LOGGER.warn("Could not calculate secondary structure (error {}). Will try to get a DSSP file from the RCSB web server instead.", e.getMessage());
			
			try {
				DSSPParser.fetch(bioJavaStruct.getPDBCode(), bioJavaStruct, true); //download from PDB the DSSP result
			} catch(Exception bige){
				LOGGER.warn("Could not get a DSSP file from RCSB web server. There will not be secondary structure assignment for this structure ({}). Error: {}", bioJavaStruct.getPDBCode(), bige.getMessage());
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
		if(experimentalTechniques==null){
			return new String[0];
		}
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
		DateFormat dateStringFormat = new SimpleDateFormat("yyyy-MM-dd");
		return dateStringFormat.format(inputDate);
	}

	/**
	 * Convert a bioassembly information into a map of transform, chainindices it relates to.
	 * @param bioassemblyInfo  the bioassembly info object for this structure
	 * @param chainIdToIndexMap the map of chain ids to the index that chain corresponds to.
	 * @return the bioassembly information (as primitive types).
	 */
	public static Map<double[], int[]> getTransformMap(BioAssemblyInfo bioassemblyInfo, Map<String, Integer> chainIdToIndexMap) {
	    Map<Matrix4d, List<Integer>> matMap = new LinkedHashMap<>();
		List<BiologicalAssemblyTransformation> transforms = bioassemblyInfo.getTransforms();
		for (BiologicalAssemblyTransformation transformation : transforms) {
			Matrix4d transMatrix = transformation.getTransformationMatrix();
			String transChainId = transformation.getChainId();
			if (!chainIdToIndexMap.containsKey(transChainId)){
				continue;
			}
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

	    Map<double[], int[]> outMap = new LinkedHashMap<>();
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
	public static double[] convertToDoubleArray(Matrix4d transformationMatrix) {
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


	/**
	 * Get a list of atoms for a group. Only add each atom once.
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
	 * Find the number of bonds in a group
	 * @param atomsInGroup the list of atoms in the group
	 * @return the number of atoms in the group
	 */
	public static int getNumBondsInGroup(List<Atom> atomsInGroup) {
		int bondCounter = 0;
		for(Atom atom : atomsInGroup) { 
			if(atom.getBonds()==null){
				continue;
			}
			for(Bond bond : atom.getBonds()) {
				// Now set the bonding information.
				Atom other = bond.getOther(atom);
				// If both atoms are in the group
				if (atomsInGroup.indexOf(other)!=-1){
					Integer firstBondIndex = atomsInGroup.indexOf(atom);
					Integer secondBondIndex = atomsInGroup.indexOf(other);
					// Don't add the same bond twice
					if (firstBondIndex<secondBondIndex){
						bondCounter++;
					}
				}
			}
		}
		return bondCounter;
	}

	/**
	 * Get the secondary structure as defined by DSSP.
	 * @param group the input group to be calculated
	 * @return the integer index of the group type.
	 */
	public static int getSecStructType(Group group) {
		SecStrucState props = (SecStrucState) group.getProperty("secstruc");
		if(props==null){
			return DsspType.NULL_ENTRY.getDsspIndex();
		}
		return DsspType.dsspTypeFromString(props.getType().name).getDsspIndex();
	}

	/**
	 * Get the secondary structure as defined by DSSP.
	 * @param group the input group to be calculated
	 * @param the integer index of the group type.
	 */
	public static void setSecStructType(Group group, int dsspIndex) {
		SecStrucType secStrucType = getSecStructTypeFromDsspIndex(dsspIndex);
		SecStrucState secStrucState = new SecStrucState(group, "MMTF_ASSIGNED", secStrucType);
		if(secStrucType!=null){
			group.setProperty("secstruc", secStrucState);
		}
		else{
		}
	}


	/**
	 * Set the DSSP type based on a numerical index.
	 * @param dsspIndex the integer index of the type to set
	 * @return the instance of the SecStrucType object holding this secondary
	 * structure type.
	 */
	public static SecStrucType getSecStructTypeFromDsspIndex(int dsspIndex) {
		String dsspType = DsspType.dsspTypeFromInt(dsspIndex).getDsspType();
		for(SecStrucType secStrucType : SecStrucType.values())
		{
			if(dsspType==secStrucType.name)
			{
				return secStrucType;
			}
		}
		// Return a null entry.
		return null;
	}

	/**
	 * Get summary information for the structure.
	 * @param structure the structure for which to get the information.
	 */
	public static MmtfSummaryDataBean getStructureInfo(Structure structure) {
		MmtfSummaryDataBean mmtfSummaryDataBean = new MmtfSummaryDataBean();
		// Get all the atoms
		List<Atom> theseAtoms = new ArrayList<>();
		List<Chain> allChains = new ArrayList<>();
		Map<String, Integer> chainIdToIndexMap = new LinkedHashMap<>();
		int chainCounter = 0;
		int bondCount = 0;
		mmtfSummaryDataBean.setAllAtoms(theseAtoms);
		mmtfSummaryDataBean.setAllChains(allChains);
		mmtfSummaryDataBean.setChainIdToIndexMap(chainIdToIndexMap);
		for (int i=0; i<structure.nrModels(); i++){
			List<Chain> chains = structure.getModel(i);
			allChains.addAll(chains);
			for (Chain chain : chains) {
				String idOne = chain.getId();
				if (!chainIdToIndexMap.containsKey(idOne)) { 
					chainIdToIndexMap.put(idOne, chainCounter);
				}
				chainCounter++;
				for (Group g : chain.getAtomGroups()) {
					for(Atom atom: getAtomsForGroup(g)){
						theseAtoms.add(atom);		
						// If both atoms are in the group
						if (atom.getBonds()!=null){
							bondCount+=atom.getBonds().size();
						}
					}
				}
			}
		}
		// Assumes all bonds are referenced twice
		mmtfSummaryDataBean.setNumBonds(bondCount/2);
		return mmtfSummaryDataBean;

	}

	/**
	 * Get a list of N 4*4 matrices from a single list of doubles of length 16*N.
	 * @param ncsOperMatrixList the input list of doubles
	 * @return the list of 4*4 matrics 
	 */
	public static Matrix4d[] getNcsAsMatrix4d(double[][] ncsOperMatrixList) {
		if(ncsOperMatrixList==null){
			return null;
		}
		int numMats = ncsOperMatrixList.length;
		if(numMats==0){
			return null;
		}
		if(numMats==1 && ncsOperMatrixList[0].length==0){
			return null;
		}
		Matrix4d[] outList = new Matrix4d[numMats];
		for(int i=0; i<numMats; i++){
			outList[i] = new Matrix4d(ncsOperMatrixList[i]);
		}
		return outList;
	}

	/**
	 * Get a list of length N*16 of a list of Matrix4d*N.
	 * @param ncsOperators the {@link Matrix4d} list 
	 * @return the list of length N*16 of the list of matrices
	 */
	public static double[][] getNcsAsArray(Matrix4d[] ncsOperators) {
		if(ncsOperators==null){
			return new double[0][0];
		}
		double[][] outList = new double[ncsOperators.length][16];
		for(int i=0; i<ncsOperators.length;i++){
			outList[i] = convertToDoubleArray(ncsOperators[i]);
		}
		return outList;
	}

	/**
	 * Insert the group in the given position in the sequence.
	 * @param chain the chain to add the seq res group to
	 * @param group the group to add
	 * @param sequenceIndexId the index to add it in
	 */
	public static void insertSeqResGroup(Chain chain, Group group, int sequenceIndexId) {
		List<Group> seqResGroups = chain.getSeqResGroups();
		addGroupAtId(seqResGroups, group, sequenceIndexId);
	}

	/**
	 * Add the missing groups to the SeqResGroups.
	 * @param modelChain the chain to add the information for
	 * @param sequence the sequence of the construct
	 */
	public static void addSeqRes(Chain modelChain, String sequence) {
		List<Group> seqResGroups = modelChain.getSeqResGroups();
		GroupType chainType = getChainType(modelChain.getAtomGroups());
		for(int i=0; i<sequence.length(); i++){
			char singleLetterCode = sequence.charAt(i);
			Group group = null;
			if(seqResGroups.size()<=i){
			}
			else{
				group=seqResGroups.get(i);
			}
			if(group!=null){
				continue;
			}
			group = getSeqResGroup(modelChain, singleLetterCode, chainType);
			addGroupAtId(seqResGroups, group, i);
			seqResGroups.set(i, group);
		}
	}

	private static GroupType getChainType(List<Group> groups) {
		for(Group group : groups) {
			if(group==null){
				continue;
			}
			else if(group.getType()!=GroupType.HETATM){
				return group.getType();
			}
		}
		return GroupType.HETATM;
	}

	private static <T> void addGroupAtId(List<T> seqResGroups, T group, int sequenceIndexId) {
		while(seqResGroups.size()<=sequenceIndexId){
			seqResGroups.add(null);
		}
		if(sequenceIndexId>=0){
			seqResGroups.set(sequenceIndexId, group);
		}		
	}
	
	private static Group getSeqResGroup(Chain modelChain, char singleLetterCode, GroupType type) {
		if(type==GroupType.AMINOACID){
			AminoAcidImpl a = new AminoAcidImpl();
			a.setRecordType(AminoAcid.SEQRESRECORD);
			a.setAminoType(singleLetterCode);
			ChemComp chemComp = new ChemComp();
			chemComp.setOne_letter_code(""+singleLetterCode);
			a.setChemComp(chemComp);
			return a;

		} else if (type==GroupType.NUCLEOTIDE) {
			NucleotideImpl n = new NucleotideImpl();
			ChemComp chemComp = new ChemComp();
			chemComp.setOne_letter_code(""+singleLetterCode);
			n.setChemComp(chemComp);
			return n;
		}
		else{
			return null;
		}
	}
}