package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;

/**
 * A MultipleAlignment is a Data Structure to store the core of a multiple structure alignment, as a return type.
 * 
 * It is described as a collection of {@link BlockSet} that define the aligned positions, 
 *                    a collection of structure identifiers (i,e. Atom arrays),
 *                    information about the alignment (Score, RMSD, etc),
 *                    and creation properties (algorithm, version, etc).
 *
 * @author Aleix Lafita
 * 
 */
public class MultipleAlignment implements Serializable, Cloneable{

	private static final long serialVersionUID = 3432043794125805139L;
	public static final String DEFAULT_ALGORITHM_NAME = "CE-MC";  //there is no algorithm yet... TODO

	//Creation Properties
	String algorithmName;
	String version;
	long ioTime;
	long calculationTime;
	long id;
	
	//Structure Identifiers
	List<String> structureNames;  			//names of the structures in PDB or SCOP format
	List<Atom[]> atomArrays;      			//arrays of atoms for every structure in the alignment
	int size;								//number of structures

	//Aligned Positions
	List<BlockSet> blockSets;				//aligned positions of the alignment
	List<String> alnSequences; 				//sequence alignment for every structure as a String with gaps (cash)
	
	//Alignment Information
	double rmsd;							//RMSD of the multiple alignment
	double tmScore;							//TM-score of the multiple alignment
	int length;								//total length of the alignment, including gaps = number of aligned positions
	int coreLength;							//length of the ungapped aligned positions = length - gap length
	double coverage;						//fraction of the total structure length aligned
	double similarity;						//similarity measure for all the structures (global value)
	
	//Algorithm Specific Information
	double algScore;						//algorithm score for the multiple alignment
	double probability;
	
	/**
	 * Constructor
	 */
	public MultipleAlignment() {
		
		algorithmName = DEFAULT_ALGORITHM_NAME;
		version = "1.0";
		
		structureNames = new ArrayList<String>();
		atomArrays = new ArrayList<Atom[]>();
		
		blockSets = new ArrayList<BlockSet>();
		alnSequences = new ArrayList<String>();
		
	}
	
	/**
	 * Copy constructor.
	 */
	public MultipleAlignment(MultipleAlignment ma) {
		
		this.setAlgorithmName(ma.getAlgorithmName());
		this.setAlgScore(ma.getAlgScore());
		this.setAlnSequences(new ArrayList<String>(ma.getAlnSequences()));
		this.setAtomArrays(new ArrayList<Atom[]>(ma.getAtomArrays()));
		
		//Ensure a proper cloning of all the BlockSet objects
		List<BlockSet> bSets = new ArrayList<BlockSet>();
		for (BlockSet bs:ma.getBlockSets()){
			bSets.add(bs.clone());
		}
		this.setBlockSets(bSets);
		
		this.setCalculationTime(ma.getCalculationTime());
		this.setCoreLength(ma.getCoreLength());
		this.setCoverage(ma.getCoverage());
		this.setId(ma.getId());
		this.setIoTime(ma.getIoTime());
		this.setLength(ma.getLength());
		
		this.setProbability(ma.getProbability());
		this.setRmsd(ma.getRmsd());
		this.setSimilarity(ma.getSimilarity());
		this.setStructureNames(ma.getStructureNames());
		this.setSize(ma.getSize());
		
		this.setTmScore(ma.getTmScore());
		this.setVersion(ma.getVersion());
		
	}
	
	/**
	 * Creates and returns a copy of this object. Uses the copy constructor.
	 */
	@Override
	public MultipleAlignment clone() {
		
		return new MultipleAlignment(this);
	}

	@Override
	public String toString() {
		return "MultipleAlignment [algorithmName=" + algorithmName
				+ ", version=" + version + ", ioTime=" + ioTime
				+ ", calculationTime=" + calculationTime + ", id=" + id
				+ ", structureNames=" + structureNames + ", atomArrays="
				+ atomArrays + ", size=" + size + ", blockSets=" + blockSets
				+ ", alnSequences=" + alnSequences + ", rmsd=" + rmsd
				+ ", tmScore=" + tmScore + ", length=" + length
				+ ", coreLength=" + coreLength + ", coverage=" + coverage
				+ ", similarity=" + similarity + ", algScore=" + algScore
				+ ", probability=" + probability + "]";
	}

	//Getters and Setters **************************************************************************************
	
	public String getAlgorithmName() {
		return algorithmName;
	}

	public void setAlgorithmName(String algorithmName) {
		this.algorithmName = algorithmName;
	}

	public String getVersion() {
		return version;
	}

	public void setVersion(String version) {
		this.version = version;
	}

	public long getIoTime() {
		return ioTime;
	}

	public void setIoTime(long ioTime) {
		this.ioTime = ioTime;
	}

	public long getCalculationTime() {
		return calculationTime;
	}

	public void setCalculationTime(long calculationTime) {
		this.calculationTime = calculationTime;
	}

	public long getId() {
		return id;
	}

	public void setId(long id) {
		this.id = id;
	}

	public List<String> getStructureNames() {
		return structureNames;
	}

	public void setStructureNames(List<String> structureNames) {
		this.structureNames = structureNames;
	}

	public List<Atom[]> getAtomArrays() {
		return atomArrays;
	}

	public void setAtomArrays(List<Atom[]> atomArrays) {
		this.atomArrays = atomArrays;
	}

	public List<BlockSet> getBlockSets() {
		return blockSets;
	}

	public void setBlockSets(List<BlockSet> blockSets) {
		this.blockSets = blockSets;
	}

	public List<String> getAlnSequences() {
		return alnSequences;
	}

	public void setAlnSequences(List<String> alnSequences) {
		this.alnSequences = alnSequences;
	}

	public double getRmsd() {
		return rmsd;
	}

	public void setRmsd(double rmsd) {
		this.rmsd = rmsd;
	}

	public double getTmScore() {
		return tmScore;
	}

	public void setTmScore(double tmScore) {
		this.tmScore = tmScore;
	}

	/**
	 * Total length of the multiple alignment, including gaps.
	 */
	public int getLength() {
		return length;
	}

	public void setLength(int length) {
		this.length = length;
	}

	public double getCoverage() {
		return coverage;
	}

	public void setCoverage(double coverage) {
		this.coverage = coverage;
	}

	public double getSimilarity() {
		return similarity;
	}

	public void setSimilarity(double similarity) {
		this.similarity = similarity;
	}

	public double getAlgScore() {
		return algScore;
	}

	public void setAlgScore(double algScore) {
		this.algScore = algScore;
	}

	public double getProbability() {
		return probability;
	}

	public void setProbability(double probability) {
		this.probability = probability;
	}
	
	/**
	 * Number of structures in the multiple alignment.
	 */
	public int getSize() {
		return size;
	}

	public void setSize(int size) {
		this.size = size;
	}

	/**
	 * Length of the ungapped aligned positions, which means that all structures have an aligned residue in those positions.
	 */
	public int getCoreLength() {
		return coreLength;
	}

	public void setCoreLength(int coreLength) {
		this.coreLength = coreLength;
	}

	
}
