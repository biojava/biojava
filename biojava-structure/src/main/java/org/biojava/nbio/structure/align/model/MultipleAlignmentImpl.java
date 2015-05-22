package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.Pose.PoseMethod;

/**
 * A general implementation of a {@link MultipleAlignment}.
 *
 * @author Aleix Lafita
 * 
 */
public class MultipleAlignmentImpl extends AbstractScoresCache implements Serializable, MultipleAlignment, Cloneable{

	private static final long serialVersionUID = 3432043794125805139L;

	MultipleAlignmentEnsemble parent;

	//Multiple Alignment Positions
	List<BlockSet> blockSets;				//aligned positions. It is a list because it can store more than one alternative MSTA. Index 0 is the optimal alignment.

	//Cache variables (can be updated)
	List<String> alnSequences; 				//sequence alignment for every structure as a String with gaps (cache)
	int length;								//total length of the alignment, including gaps: sum of BlockSet lengths (cache)
	int coreLength;							//length of the alignment without gaps, only columns without gaps: sum of BlockSet core lengths (cache)
	List<Matrix4d> pose;

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
	 * @param structureNames List of Structure identifiers.
	 * @return MultipleAlignment a MultipleAlignment instance part of its own MultipleAlignmentEnsemble.
	 */
	public MultipleAlignmentImpl(List<String> structureNames) {
		this(new MultipleAlignmentEnsembleImpl(structureNames));
	}
	
	/**
	 * Constructor linking to an existing ensemble.
	 * 
	 * Automatically adds this alignment to the parent ensemble
	 * @param ensemble parent MultipleAlignmentEnsemble.
	 * @return MultipleAlignment a MultipleAlignment instance part of an ensemble.
	 */
	public MultipleAlignmentImpl(MultipleAlignmentEnsemble ensemble) {
		super();

		parent = ensemble;
		if (parent!=null) parent.getMultipleAlignments().add(this);

		blockSets = null;
		alnSequences = null;
		pose = null;

		length = -1;						//Value -1 reserved to indicate that has to be calculated
		coreLength = -1;
	}

	/**
	 * Copy constructor.
	 * @param ma MultipleAlignmentImpl to copy.
	 * @return MultipleAlignmentImpl identical copy of the input MultipleAlignmentImpl.
	 */
	public MultipleAlignmentImpl(MultipleAlignmentImpl ma) {
		super(ma);
		
		parent = ma.parent;
		
		alnSequences = null;
		if (ma.alnSequences!=null) alnSequences = new ArrayList<String>(ma.alnSequences);
		
		pose = null;  //Because the pose is a cache variable it has to be updated/calculated again.
		
		blockSets = null;
		if (ma.blockSets!=null){
			//Make a deep copy of everything
			this.blockSets = new ArrayList<BlockSet>();
			for (BlockSet bs:ma.blockSets){
				BlockSet newBS = bs.clone();
				newBS.setMultipleAlignment(this); //This automatically adds the newBS to the blockSets list
			}
		}
		
		
	}
	
	@Override
	public MultipleAlignmentImpl clone() {
		return new MultipleAlignmentImpl(this);
	}

	@Override
	public String toString() {
		return "MultipleAlignmentImpl [parent=" + parent + ", blockSets="
				+ blockSets + ", alnSequences=" + alnSequences + ", length="
				+ length + ", coreLength=" + coreLength + ", pose=" + pose
				+ "]";
	}

	@Override
	public List<BlockSet> getBlockSets() {
		if (blockSets == null) blockSets = new ArrayList<BlockSet>();
		return blockSets;
	}
	
	/**
	 * Convenience method to get a list of all blocks from all blocksets
	 * @return
	 * @throws StructureAlignmentException
	 */
	@Override
	public List<Block> getBlocks() {
		List<Block> blocks = new ArrayList<Block>();
		for(BlockSet bs : getBlockSets()) {
			blocks.addAll(bs.getBlocks());
		}
		return blocks;
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

	public void updateAlnSequences() {
		// TODO Auto-generated method stub
		
	}

	/**
	 * Returns a transformation matrix for each structure giving the
	 * 3D superimposition information of the multiple structure alignment.
	 * @return the 3D superimposition information of the alignment
	 * @throws StructureAlignmentException 
	 */
	@Override
	public List<Matrix4d> getTransformations()
			throws StructureAlignmentException {
		return pose;
	}

	/**
	 * Set a new superposition for the structures.
	 * 
	 * This may trigger other properties to update which depend on the superposition.
	 * @param matrices
	 */
	@Override
	public void setTransformations(List<Matrix4d> matrices) throws StructureAlignmentException {
		if(size() != matrices.size()) {
			throw new StructureAlignmentException("Wrong number of structures for this alignment");
		}
		// set properties that depend on the pose to null
		//TODO clearPose();
		pose = matrices;
	}

	@Override
	public int size() throws StructureAlignmentException {
		return parent.size();
	}

	@Override
	public int length() throws StructureAlignmentException {
		if (length == -1) updateLength();
		return length;
	}

	@Override
	public int getCoreLength() throws StructureAlignmentException {
		if (coreLength == -1) updateCoreLength();
		return coreLength;
	}

	protected void updateLength() throws StructureAlignmentException {
		if(getBlockSets().size()==0) throw new StructureAlignmentException("Empty MultipleAlignment: getBlockSetNum() == 0.");
		//Try to calculate it from the BlockSet information
		else {
			length = 0;
			for (BlockSet blockSet:blockSets) length += blockSet.length();
		}
	}

	protected void updateCoreLength() throws StructureAlignmentException {
		if(getBlockSets().size()==0) throw new StructureAlignmentException("Empty MultipleAlignment: getBlockSetNum() == 0.");
		//Try to calculate it from the BlockSet information
		else {
			coreLength = 0;
			for (BlockSet blockSet:blockSets) coreLength += blockSet.getCoreLength();
		}
	}

	protected void updateCache(PoseMethod method) throws StructureAlignmentException, StructureException {
		//TODO updatePose(method);
		updateAlnSequences();
		updateCoreLength();
		updateLength();
	}

	@Override
	public MultipleAlignmentEnsemble getEnsemble() {
		return parent;
	}

	@Override
	public void setEnsemble(MultipleAlignmentEnsemble parent) {
		//Delete the alignment instance from the parent list
		if (parent!=null) parent.getMultipleAlignments().remove(this);
		this.parent = parent;
		//Cross-link the two objects automatically when parent is changed.
		if (parent!=null) parent.getMultipleAlignments().add(this);
	}
	
}
