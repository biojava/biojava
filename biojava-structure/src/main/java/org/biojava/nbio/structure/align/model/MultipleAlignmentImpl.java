package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.Pose.PoseMethod;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * A general implementation of a {@link MultipleAlignment}.
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
	double algScore;						//algorithm score for this multiple alignment
	double probability;
	
	/**
	 * Default Constructor. Empty alignment. No structures assigned.
	 * @return MultipleAlignment an empty MultipleAlignment instance.
	 */
	public MultipleAlignmentImpl() {
		this(new MultipleAlignmentEnsembleImpl());  //assign an empty ensemble.
	}
	
	/**
	 * Constructor that allows creating a MultipleAlignment instance without caring about the ensemble.
	 * Ideal for dealing with one MultipleAlignment alternative only.
	 * @param atomArrays List of Atom arrays of the structures.
	 * @return MultipleAlignment a MultipleAlignment instance part of an MultipleAlignmentEnsemble.
	 */
	public MultipleAlignmentImpl(List<Atom[]> atomArrays) {
		this(new MultipleAlignmentEnsembleImpl(atomArrays));
	}
	
	/**
	 * Constructor linking to an existing ensemble.
	 * @param ensemble parent MultipleAlignmentEnsemble.
	 * @return MultipleAlignment a MultipleAlignment instance part of an ensemble.
	 */
	public MultipleAlignmentImpl(MultipleAlignmentEnsemble ensemble) {
		
		parent = ensemble;
		//TODO is this a good idea? Cross-link the two objects automatically when instanciated.
		if (parent!=null) parent.getMultipleAlignments().add(this);
		
		blockSets = null;
		alnSequences = null;
		pose = null;
		
		algScore = 0;
		probability = 0;
		
		length = -1;						//Value -1 reserved to indicate that has to be calculated
		coreLength = -1;
	}
	
	/**
	 * Constructor from an AFPChain instance. Creates an equivalent pairwise alignment.
	 * @param ensemble parent MultipleAlignmentEnsemble.
	 * @return MultipleAlignment a MultipleAlignment instance part of an MultipleAlignmentEnsemble.
	 * @throws StructureAlignmentException 
	 * @throws StructureException 
	 */
	public MultipleAlignmentImpl(AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureAlignmentException, StructureException {
		//Copy all the creation and algorithm information
		this(new MultipleAlignmentEnsembleImpl(Arrays.asList(ca1,ca2)));
		parent.setAlgorithmName(afpChain.getAlgorithmName());
		parent.setVersion(afpChain.getVersion());
		parent.setId(afpChain.getId());
		parent.setCalculationTime(afpChain.getCalculationTime());
		algScore = afpChain.getAlignScore();
		probability = afpChain.getProbability();
		
		//Create a BlockSet for every block in AFPChain TODO we dont't know if blocks correspond to flexible or CP.
		for (int bs=0; bs<afpChain.getBlockNum(); bs++){
			BlockSet blockSet = new BlockSetImpl(this);
			Block block = new BlockImpl(blockSet);
			block.getAlignRes().add(new ArrayList<Integer>()); //add the two chains
			block.getAlignRes().add(new ArrayList<Integer>());
			
			for (int i=0; i<afpChain.getOptLen()[bs]; i++){
				block.getAlignRes().get(0).add(afpChain.getOptAln()[bs][0][i]);
				block.getAlignRes().get(1).add(afpChain.getOptAln()[bs][1][i]);
			}
		}
		updateCache(PoseMethod.REFERENCE);
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
				+ blockSets + ", alnSequences=" + alnSequences + ", length="
				+ length + ", coreLength=" + coreLength + ", pose=" + pose
				+ ", algScore=" + algScore + ", probability=" + probability
				+ "]";
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
	public List<Matrix> getDistanceMatrix() throws StructureAlignmentException {
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

	@Override
	public MultipleAlignmentEnsemble getParent() {
		return parent;
	}

	@Override
	public void setParent(MultipleAlignmentEnsemble parent) {
		this.parent = parent;
		//TODO is this a good idea? Cross-link the two objects automatically when parent is changed.
		if (parent!=null) parent.getMultipleAlignments().add(this);
	}
	
}
