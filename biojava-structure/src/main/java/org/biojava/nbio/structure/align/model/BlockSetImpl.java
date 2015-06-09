package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.StructureException;

/**
 * A general implementation of a BlockSet to store multiple alignments.
 *
 * @author Aleix Lafita
 * 
 */
public class BlockSetImpl extends AbstractScoresCache implements Serializable, BlockSet, Cloneable{

	private static final long serialVersionUID = -1015791986000076089L;
	
	private MultipleAlignment parent;		//link to the MultipleAlignment parent instance
	private List<Block> blocks;				//aligned positions as a list of Blocks
	
	//Cache variables (can be updated)
	private List<Matrix4d> pose;			//3D superimposition information pose
	private int length;						//total number of aligned positions, including gaps = sum of blocks lengths (cache)
	private int coreLength;					//number of aligned positions without gaps (cache)
	
	/**
	 * Constructor.
	 * @param multipleAlignment the parent MultipleAlignment of the BlockImpl instance.
	 * @return BlockSetImpl a BlockSetImpl instance linked to its parent MultipleAlignment.
	 * @throws StructureAlignmentException 
	 */
	public BlockSetImpl(MultipleAlignment multipleAlignment) {
		
		parent = multipleAlignment;
		if (parent!=null) parent.getBlockSets().add(this);
		blocks = null;
		
		//Cache variables (can be updated)
		pose = null;
		length = -1;						//Value -1 reserved to indicate that has to be calculated
		coreLength = -1;
	}
	
	/**
	 * Copy constructor.
	 * @param bs BlockSetImpl object to be copied.
	 * @return BlockSetImpl an identical copy of the input BlockSetImpl object.
	 */
	public BlockSetImpl(BlockSetImpl bs){
		
		this.parent = bs.parent;
		this.length = bs.length;
		this.coreLength = bs.coreLength;
		
		this.pose = null;  pose = null;  //Because the pose is a cache variable it has to be updated/calculated again.
		
		blocks = null;
		if (bs.blocks!=null){
			//Make a deep copy of everything
			this.blocks = new ArrayList<Block>();
			for (Block b:bs.blocks){
				Block newB = b.clone();
				newB.setBlockSet(this); //This automatically adds the newB to the blocks list
			}
		}
	}
	
	/**
	 * Clear scores and cached properties. Recursively clears member blocks.
	 */
	@Override
	public void clear() {
		super.clear();
		length = -1;
		coreLength = -1;
		for(Block a : getBlocks()) {
			a.clear();
		}
	}
	
	@Override
	public BlockSetImpl clone(){
		return new BlockSetImpl(this);
	}
	
	@Override
	public String toString() {
		return "BlockSetImpl [blocks=" + blocks
				+ ", pose=" + pose + ", length=" + length + ", coreLength="
				+ coreLength + "]";
	}

	@Override
	public MultipleAlignment getMultipleAlignment() {
		return parent;
	}

	@Override
	public void setMultipleAlignment(MultipleAlignment parent) {
		this.parent = parent;
	}

	@Override
	public List<Block> getBlocks() {
		if (blocks==null) blocks = new ArrayList<Block>();
		return blocks;
	}

	@Override
	public void setBlocks(List<Block> blocks) {
		this.blocks = blocks;
		for(Block b:blocks) {
			b.setBlockSet(this);
		}
	}
	
	/**
	 * Returns a transformation matrix for each structure giving the
	 * 3D superimposition information of the multiple structure alignment.
	 * @return the 3D superimposition information of the alignment
	 * @throws StructureAlignmentException 
	 */
	@Override
	public List<Matrix4d> getTransformations() throws StructureAlignmentException {
		return pose;
	}
	
	/**
	 * Set a new superposition for the structures.
	 * 
	 * This may trigger other properties to update which depend on the superposition.
	 * @param matrices
	 */
	@Override
	public void setTransformations(List<Matrix4d> transformations) throws StructureAlignmentException {
		pose = transformations;
	}

	@Override
	public int length() throws StructureAlignmentException {
		if (length == -1) updateLength();
		return length;
	}

	@Override
	public int size() throws StructureAlignmentException {
		//Get the size from the variables that can contain the information
		if (parent != null) return parent.size();
		else if (getBlocks().size()==0) throw new StructureAlignmentException("Empty BlockSet: parent==null && getBlockNum() == 0.");
		else return blocks.get(0).size();
	}

	@Override
	public int getCoreLength() throws StructureAlignmentException {
		if (coreLength == -1) updateCoreLength();
		return coreLength;
	}

	protected void updateLength() throws StructureAlignmentException {
		if(getBlocks().size()==0) throw new StructureAlignmentException("Empty BlockSet: getBlockNum() == 0.");
		//Try to calculate it from the Block information
		else {
			length = 0;
			for (Block block:blocks) length += block.length();
		}
	}

	protected void updateCoreLength() throws StructureAlignmentException {
		if(getBlocks().size()==0) throw new StructureAlignmentException("Empty BlockSet: getBlockNum() == 0.");
		//Try to calculate it from the Block information
		else {
			coreLength = 0;
			for (Block block:blocks) coreLength += block.getCoreLength();
		}
	}

	protected void updateCache() throws StructureAlignmentException, StructureException {
		updateCoreLength();
		updateLength();
	}
	
}
