package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.jama.Matrix;

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

	//Creation Properties
	String algorithmName;
	String version;
	long ioTime;
	long id;
	
	//Structure Identifiers
	List<String> structureNames;  			//names of the structures in PDB or SCOP format
	List<Atom[]> atomArrays;      			//arrays of atoms for every structure in the alignment
	List<Matrix> distanceMatrix; 			//A list of n (l*l)-matrices that store the distance between every pair of residues for every protein
											//n=nr. structures; l=length of the protein

	//Aligned Positions
	List<BlockSet> blockSets;				//aligned positions. It is a list because it can store more than one alternative MSTA. Index 0 is the optimal alignment.
	List<String> alnSequences; 				//sequence alignment for every structure as a String with gaps (cash)
	
	//Alignment Information - TODO This information seems more reasonable to be in Pose or BlockSet, because it is BlockSet specific.
	double rmsd;							//RMSD of the optimal multiple alignment
	double tmScore;							//TM-score of the optimal multiple alignment
	int length;								//total length of the optimal alignment, including gaps = number of aligned positions
	double coverage;						//fraction of the total structure length aligned
	double similarity;						//similarity measure for all the structures (global value)
	
	//Algorithm Specific Information
	double algScore;						//algorithm score for the multiple alignment
	double probability;
	long calculationTime;					//running time of the algorithm
	
	/**
	 * Constructor.
	 * @return MultipleAlignment a MultipleAlignment instance.
	 */
	public MultipleAlignment() {
		
		algorithmName = null;
		version = "1.0";
		ioTime = 0;
		id = 0;
		
		structureNames = null;
		atomArrays = null;
		distanceMatrix = null;
		
		blockSets = null;
		alnSequences = null;
		
		rmsd = 0;
		tmScore = 0;
		length = 0;
		coverage =0;
		similarity = 0;
		
		algScore = 0;
		probability = 0;
		calculationTime = 0;
	}
	
	/**
	 * Copy constructor.
	 * @return MultipleAlignment identical copy of the input MultipleAlignment.
	 */
	public MultipleAlignment(MultipleAlignment ma) {
		
		algorithmName = ma.getAlgorithmName();
		algScore = ma.getAlgScore();
		alnSequences = new ArrayList<String>(ma.getAlnSequences());
		atomArrays = new ArrayList<Atom[]>(ma.getAtomArrays());
		distanceMatrix = new ArrayList<Matrix>(ma.getDistanceMatrix());
		
		blockSets = null;
		if (ma.getBlockSets()!=null){
			//Make a deep copy of everything
			this.blockSets = new ArrayList<BlockSet>();
			for (BlockSet bs:ma.getBlockSets()){
				blockSets.add((BlockSet) bs.clone());
			}
		}
		
		structureNames = null;
				
		rmsd = ma.getRmsd();
		tmScore = ma.getTmScore();
		length = ma.length();
		coverage = ma.getCoverage();
		similarity = ma.getSimilarity();
		
		algScore = ma.getAlgScore();
		probability = ma.getProbability();
		calculationTime = ma.getCalculationTime();
		
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
				+ ", version=" + version + ", ioTime=" + ioTime + ", id=" + id
				+ ", structureNames=" + structureNames + ", atomArrays="
				+ atomArrays + ", distanceMatrix=" + distanceMatrix
				+ ", blockSets=" + blockSets + ", alnSequences=" + alnSequences
				+ ", rmsd=" + rmsd + ", tmScore=" + tmScore + ", length="
				+ length + ", coverage=" + coverage + ", similarity="
				+ similarity + ", algScore=" + algScore + ", probability="
				+ probability + ", calculationTime=" + calculationTime + "]";
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
	 * Returns the number of aligned positions in the optimal alignment.
	 * @return int number of aligned positions
	 * @see #size()
	 * @see #getBlockSetNum()
	 */
	public int length() {
		return length;
	}
	
	/**
	 * Returns the number of aligned structures in the MultipleAlignment.
	 * @return int number of aligned structures
	 * @see #length()
	 * @see #getBlockSetNum()
	 */
	public int size() {
		if (atomArrays == null) return 0;
		else return atomArrays.size();
	}
	
	/**
	 * Returns the number of BlockSets in the MultipleAlignment.
	 * @return int number of BlockSets
	 * @see #length()
	 * @see #size()
	 */
	public int getBlockSetNum(){
		if (blockSets == null) return 0;
		else return blockSets.size();
	}

	public List<Matrix> getDistanceMatrix() {
		return distanceMatrix;
	}

	public void setDistanceMatrix(List<Matrix> distanceMatrix) {
		this.distanceMatrix = distanceMatrix;
	}

}
