package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;

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
	
	
	/**
	 * Clear scores and cached properties. Recursively clears member blocksets.
	 */
	@Override
	public void clear() {
		super.clear();
		alnSequences = null;
		length = -1;
		coreLength = -1;
		for(BlockSet a : getBlockSets()) {
			a.clear();
		}
	}
	
	@Override
	public MultipleAlignmentImpl clone() {
		return new MultipleAlignmentImpl(this);
	}


	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "MultipleAlignmentImpl [blockSets=" + blockSets + ", length="
				+ length + ", coreLength=" + coreLength + "]";
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
		clear();
		pose = matrices;
	}

	@Override
	public int size() throws StructureAlignmentException {
		return parent.size();
	}

	@Override
	public int length() throws StructureAlignmentException {
		if (length < 0 ) updateLength();
		return length;
	}

	@Override
	public int getCoreLength() throws StructureAlignmentException {
		if (coreLength < 0) updateCoreLength();
		return coreLength;
	}

	/**
	 * Forces recalculation of the AlnSequences based on the alignment blocks.
	 */
	public void updateAlnSequences() {
		
		//This method is only to try the aligPanel. It does not calculate the correct sequence alignment
		alnSequences = new ArrayList<String>();
		try {
			List<Atom[]> atoms = getEnsemble().getAtomArrays();
			for (int row=0; row<size(); row++){
				String seq = "";
				for (BlockSet bs :getBlockSets()){
					for (Block b: bs.getBlocks()){
						for (int res=0; res<b.length(); res++){
							if (b.getAlignRes().get(row).get(res) != null)
								seq += StructureTools.get1LetterCode(atoms.get(row)[res].getGroup().getPDBName());
							else seq += "-";
						}
					}
				}
				alnSequences.add(seq);
			}
		} catch (StructureAlignmentException e) {
			e.printStackTrace();
		}
	}


	/**
	 * Force recalculation of the length (aligned columns) based on the block lengths
	 * @throws StructureAlignmentException
	 */
	protected void updateLength() throws StructureAlignmentException {
		if(getBlockSets().size()==0) throw new StructureAlignmentException("Empty MultipleAlignment: getBlockSetNum() == 0.");
		//Try to calculate it from the BlockSet information
		else {
			length = 0;
			for (BlockSet blockSet:blockSets) length += blockSet.length();
		}
	}

	/**
	 * Force recalculation of the core length (ungapped columns) based on the block core lengths
	 * @throws StructureAlignmentException
	 */
	protected void updateCoreLength() throws StructureAlignmentException {
		if(getBlockSets().size()==0) throw new StructureAlignmentException("Empty MultipleAlignment: getBlockSetNum() == 0.");
		//Try to calculate it from the BlockSet information
		else {
			coreLength = 0;
			for (BlockSet blockSet:blockSets) coreLength += blockSet.getCoreLength();
		}
	}

	/**
	 * Updates all cached properties
	 * @throws StructureAlignmentException
	 * @throws StructureException
	 */
	protected void updateCache() throws StructureAlignmentException, StructureException {
		updateAlnSequences();
		updateCoreLength();
		updateLength();
	}

	@Override
	public MultipleAlignmentEnsemble getEnsemble() {
		return parent;
	}

	/** 
     * Set the back-reference to its parent Ensemble.
     * 
     * Neither removes this alignment from its previous ensemble, if any, nor
     * adds it to the new parent. Calling code should assure that links to
     * and from the ensemble are consistent and free of memory leaks.
     * @param parent the parent MultipleAlignmentEnsemble.
     * @see #getEnsemble()
     */
	@Override
	public void setEnsemble(MultipleAlignmentEnsemble parent) {
		this.parent = parent;
	}

}
