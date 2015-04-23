package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * A BlockSet is a Data Structure to store the aligned positions of a multiple alignment as a collection of {@link Block}.
 * It allows flexible alignments, non-sequential alignments and circular permutations, thanks to the multiple block format.
 * It is part of a {@link MultipleAlignment} instance, named as parent.
 *
 * @author Aleix Lafita
 * 
 */
public class BlockSet implements Serializable, Cloneable{

	private static final long serialVersionUID = -1015791986000076089L;
	
	MultipleAlignment parent;		//Link to the MultipleAlignment Instance as parent
	List<Block> blocks;				//Aligned positions as a list of Blocks
	Pose pose;						//Superimposition pose (cash)
	int length;						//total number of aligned positions, including gaps = sum of columns in the blocks
	
	/**
	 * Constructor.
	 * @param multipleAlignment: the parent MultipleAlignment instance.
	 */
	public BlockSet(MultipleAlignment multipleAlignment){
		
		parent = multipleAlignment;
		blocks = new ArrayList<Block>();
		pose = new Pose(this);
		
	}
	
	/**
	 * Copy constructor.
	 * @param bs: the BlockSet from which the information is copied.
	 */
	public BlockSet(BlockSet bs){
		
		this.parent = bs.parent;
		this.pose = bs.pose.clone();
		this.length = bs.length;
		
		//Ensure a proper cloning of all the Block objects
		List<Block> blocks = new ArrayList<Block>();
		for (Block b:bs.blocks){
			blocks.add(b.clone());
		}
		this.blocks = blocks;
		
	}
	
	/**
	 * Creates and returns a copy of this object. Uses the copy constructor.
	 */
	public BlockSet clone(){
		
		return new BlockSet(this);
	}
	
	@Override
	public String toString() {
		return "BlockSet [parent=" + parent + ", blocks=" + blocks + ", pose="
				+ pose + ", length=" + length + "]";
	}

	//Getters and Setters **************************************************************************************
	
	public MultipleAlignment getMultipleAlignment() {
		return parent;
	}
	
	public void setMultipleAlignment(MultipleAlignment parent) {
		this.parent = parent;
	}
	
	public List<Block> getBlocks() {
		return blocks;
	}
	
	public void setBlocks(List<Block> blocks) {
		this.blocks = blocks;
	}
	
	public Pose getPose() {
		return pose;
	}
	
	public void setPose(Pose pose) {
		this.pose = pose;
	}

	/**
	 * Total number of aligned positions, including gaps, in the blocks (sum of columns in all blocks)
	 */
	public int getLength() {
		return length;
	}

	public void setLength(int length) {
		this.length = length;
	}
	
	/**
	 * Return the number of Blocks in the BlockSet.
	 */
	public int getBlockNum() {
		return blocks.size();
	}
}
