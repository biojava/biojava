package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.Pose.PoseMethod;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * A general implementation of a MultipleAlignment.
 *
 * @author Aleix Lafita
 * 
 */
public class MultipleAlignmentImpl implements Serializable, MultipleAlignment{

	private static final long serialVersionUID = 3432043794125805139L;

	MultipleAlignmentEnsemble parent;
	
	//Multiple Alignment Positions
	List<BlockSet> blockSets;				//aligned positions. It is a list because it can store more than one alternative MSTA. Index 0 is the optimal alignment.
	
	//Cache variables (can be updated)
	List<String> alnSequences; 				//sequence alignment for every structure as a String with gaps (cache)
	int length;								//total length of the alignment, including gaps: sum of BlockSet lengths (cache)
	int coreLength;							//length of the alignment without gaps, only columns without gaps: sum of BlockSet core lengths (cache)
	Pose pose;
	
	//Algorithm Specific Information
	double algScore;						//algorithm score for the multiple alignment
	double probability;
	
	/**
	 * Constructor.
	 * @param ensemble parent EnsembleMSTA.
	 * @return MultipleAlignment a MultipleAlignment instance part of an EnsembleMSTA.
	 */
	public MultipleAlignmentImpl(MultipleAlignmentEnsemble ensemble) {
		
		parent = ensemble;
		
		blockSets = null;
		alnSequences = null;
		
		algScore = 0;
		probability = 0;
		
		length = -1;						//Value -1 reserved to indicate that has to be calculated
		coreLength = -1;
	}
	
	/**
	 * Copy constructor.
	 * @param ma MultipleAlignmentImpl to copy.
	 * @return MultipleAlignmentImpl identical copy of the input MultipleAlignmentImpl.
	 */
	public MultipleAlignmentImpl(MultipleAlignmentImpl ma) {
		
		parent = ma.parent;
		
		alnSequences = new ArrayList<String>(ma.getAlnSequences());
		
		blockSets = null;
		if (ma.blockSets!=null){
			//Make a deep copy of everything
			this.blockSets = new ArrayList<BlockSet>();
			for (BlockSet bs:ma.blockSets){
				blockSets.add((BlockSet) bs.clone());
			}
		}
		
		algScore = ma.getAlgScore();
		probability = ma.getProbability();		
	}
	
	@Override
	public Object clone() {
		return new MultipleAlignmentImpl(this);
	}

	@Override
	public String toString() {
		return "MultipleAlignmentImpl [parent=" + parent + ", blockSets="
				+ blockSets + ", alnSequences=" + alnSequences + ", pose="
				+ pose + ", algScore=" + algScore + ", probability="
				+ probability + "]";
	}

	@Override
	public String getAlgorithmName() {
		return parent.getAlgorithmName();
	}

	@Override
	public String getVersion() {
		return parent.getVersion();
	}

	@Override
	public long getIoTime() {
		return parent.getIoTime();
	}

	@Override
	public long getCalculationTime() {
		return parent.getCalculationTime();
	}

	@Override
	public long getId() {
		return parent.getId();
	}

	@Override
	public List<String> getStructureNames() {
		return parent.getStructureNames();
	}

	@Override
	public List<Atom[]> getAtomArrays() {
		return parent.getAtomArrays();
	}

	@Override
	public List<BlockSet> getBlockSets() {
		if (blockSets == null) blockSets = new ArrayList<BlockSet>();
		return blockSets;
	}

	@Override
	public void setBlockSets(List<BlockSet> blockSets) {
		this.blockSets = blockSets;
	}

	@Override
	public List<String> getAlnSequences() {
		if (alnSequences == null) updateAlnSequences();
		return alnSequences;
	}

	@Override
	public void updateAlnSequences() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Pose getPose() throws StructureAlignmentException {
		if (pose == null) pose = new PoseMSA(this);
		return pose;
	}

	@Override
	public void updatePose(PoseMethod method) throws StructureException, StructureAlignmentException {
		if (pose == null) pose = new PoseMSA(this);
		pose.updatePose(method);
	}

	@Override
	public double getAlgScore() {
		return algScore;
	}

	@Override
	public void setAlgScore(double algScore) {
		this.algScore = algScore;
	}

	@Override
	public double getProbability() {
		return probability;
	}

	@Override
	public void setProbability(double probability) {
		this.probability = probability;
	}

	@Override
	public int size() throws StructureAlignmentException {
		return parent.size();
	}
	
	@Override
	public List<Matrix> getDistanceMatrix() {
		return parent.getDistanceMatrix();
	}

	@Override
	public int length() throws StructureAlignmentException {
		if (length == -1) updateLength();
		return length;
	}

	@Override
	public int getBlockSetNum() throws StructureAlignmentException {
		if (blockSets == null) throw new StructureAlignmentException("Empty BlockSet: getBlocks() == null.");
		else return blockSets.size();
	}

	@Override
	public int getCoreLength() throws StructureAlignmentException {
		if (coreLength == -1) updateCoreLength();
		return coreLength;
	}

	@Override
	public void updateLength() throws StructureAlignmentException {
		if(getBlockSetNum()==0) throw new StructureAlignmentException("Empty MultipleAlignment: getBlockSetNum() == 0.");
		//Try to calculate it from the BlockSet information
		else {
			length = 0;
			for (BlockSet blockSet:blockSets) length += blockSet.length();
		}
	}

	@Override
	public void updateCoreLength() throws StructureAlignmentException {
		if(getBlockSetNum()==0) throw new StructureAlignmentException("Empty MultipleAlignment: getBlockSetNum() == 0.");
		//Try to calculate it from the BlockSet information
		else {
			coreLength = 0;
			for (BlockSet blockSet:blockSets) coreLength += blockSet.getCoreLength();
		}
	}

	@Override
	public void updateCache(PoseMethod method) throws StructureAlignmentException, StructureException {
		updatePose(method);
		updateAlnSequences();
		updateCoreLength();
		updateLength();
	}
	
}
