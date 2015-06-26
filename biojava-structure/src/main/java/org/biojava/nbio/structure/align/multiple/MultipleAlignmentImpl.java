package org.biojava.nbio.structure.align.multiple;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.StructureException;

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
		pose = null;
		if (ma.pose != null){
			//Make a deep copy of everything
			this.pose = new ArrayList<Matrix4d>();
			for (Matrix4d trans:ma.pose){
				Matrix4d newTrans = (Matrix4d) trans.clone();
				pose.add(newTrans);
			}
		}
		
		length = -1;
		coreLength = -1;
		
		blockSets = null;
		if (ma.blockSets!=null){
			//Make a deep copy of everything
			this.blockSets = new ArrayList<BlockSet>();
			for (BlockSet bs:ma.blockSets){
				BlockSet newBS = bs.clone();
				newBS.setMultipleAlignment(this);
				this.blockSets.add(newBS);
			}
		}
	}
	
	@Override
	public void clear() {
		super.clear();
		length = -1;
		coreLength = -1;
		pose = null;
		for(BlockSet a : getBlockSets()) {
			a.clear();
		}
	}
	
	@Override
	public MultipleAlignmentImpl clone() {
		return new MultipleAlignmentImpl(this);
	}
	
	@Override
	public String toString() {
		String resume = "Structures:" + parent.getStructureNames() + 
				" \nAlgorithm:" + parent.getAlgorithmName() + "_" + parent.getVersion() + 
				" \nBlockSets: "+ getBlockSets().size() + 
				" \nBlocks: " + getBlocks().size() +
				" \nLength: " + length() +
				" \nCore Length: "+ getCoreLength();
		for (String score:getScores()){
			resume += " \n"+score+": "+ String.format("%.2f", getScore(score));
		}
		return resume;
	}

	@Override
	public List<BlockSet> getBlockSets() {
		if (blockSets == null) blockSets = new ArrayList<BlockSet>();
		return blockSets;
	}
	
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
	public List<Matrix4d> getTransformations(){
		return pose;
	}

	@Override
	public void setTransformations(List<Matrix4d> matrices) {
		if(size() != matrices.size()) throw new IllegalArgumentException("Wrong number of structures for this alignment");
		clear();
		pose = matrices;
	}

	@Override
	public int size() {
		return parent.size();
	}

	@Override
	public int length() {
		if (length < 0 ) updateLength();
		return length;
	}

	@Override
	public int getCoreLength() {
		if (coreLength < 0) updateCoreLength();
		return coreLength;
	}


	/**
	 * Force recalculation of the length (aligned columns) based on the block lengths
	 */
	protected void updateLength() {
		if(getBlockSets().size()==0) throw new IndexOutOfBoundsException("Empty MultipleAlignment: getBlockSetNum() == 0.");
		//Try to calculate it from the BlockSet information
		else {
			length = 0;
			for (BlockSet blockSet:blockSets) length += blockSet.length();
		}
	}

	/**
	 * Force recalculation of the core length (ungapped columns) based on the block core lengths
	 */
	protected void updateCoreLength() {
		if(getBlockSets().size()==0) throw new IndexOutOfBoundsException("Empty MultipleAlignment: getBlockSetNum() == 0.");
		//Try to calculate it from the BlockSet information
		else {
			coreLength = 0;
			for (BlockSet blockSet:blockSets) coreLength += blockSet.getCoreLength();
		}
	}

	/**
	 * Updates all cached properties
	 * @throws StructureException
	 */
	protected void updateCache() {
		updateCoreLength();
		updateLength();
	}

	@Override
	public MultipleAlignmentEnsemble getEnsemble() {
		return parent;
	}
	
	@Override
	public void setEnsemble(MultipleAlignmentEnsemble parent) {
		this.parent = parent;
	}
}
