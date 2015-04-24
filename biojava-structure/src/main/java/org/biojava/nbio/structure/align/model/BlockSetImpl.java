package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * A general implementation of a BlockSet to store multiple alignments.
 *
 * @author Aleix Lafita
 * 
 */
public class BlockSetImpl implements Serializable, BlockSet{

	private static final long serialVersionUID = -1015791986000076089L;
	
	private MultipleAlignment parent;		//link to the MultipleAlignment parent instance
	private List<Block> blocks;				//aligned positions as a list of Blocks
	private Pose pose;						//3D superimposition information pose (cache)
	private int length;						//total number of aligned positions, including gaps = sum of blocks lengths
	
	/**
	 * Constructor.
	 * @param multipleAlignment the parent MultipleAlignment of the BlockImpl instance.
	 * @return BlockSetImpl a BlockSetImpl instance linked to its parent MultipleAlignment.
	 */
	public BlockSetImpl(MultipleAlignment multipleAlignment){
		
		parent = multipleAlignment;
		blocks = null;
		pose = new PoseImpl(this);
		length = -1;						//Value -1 reserved to indicate that has to be calculated TODO correct?
	}
	
	/**
	 * Copy constructor.
	 * @param bs BlockSetImpl object to be copied.
	 * @return BlockSetImpl an identical copy of the input BlockSetImpl object.
	 */
	public BlockSetImpl(BlockSetImpl bs){
		
		this.parent = bs.getMultipleAlignment();
		this.pose = (Pose) bs.getPose().clone();
		this.length = bs.length();
		
		//Ensure a proper cloning of all the Block objects
		List<Block> blocks = new ArrayList<Block>();
		for (Block b:bs.getBlocks()){
			blocks.add((Block) b.clone());
		}
		this.blocks = blocks;
	}
	
	@Override
	public Object clone(){
		return new BlockSetImpl(this);
	}
	
	@Override
	public String toString() {
		return "BlockSetImpl [parent=" + parent + ", blocks=" + blocks
				+ ", pose=" + pose + ", length=" + length + "]";
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
		return blocks;
	}

	@Override
	public void setBlocks(List<Block> blocks) {
		this.blocks = blocks;
	}

	@Override
	public Pose getPose() {
		return pose;
	}

	@Override
	public void setPose(Pose pose) {
		this.pose = pose;
	}

	@Override
	public int length() {
		if (length == -1){
			if (blocks == null || getBlockNum()==0) return 0;
			//Try to calculate it from the block information
			else {
				length = 0;
				for (Block block:blocks) length += block.length();
			}
		}
		return length;
	}

	@Override
	public int size() {
		//Get the size from the variables that can contain the information
		if (parent != null) return parent.size();
		else if (blocks != null && getBlockNum()>0) return blocks.get(0).size();
		else return 0;
	}

	@Override
	public int getBlockNum() {
		if (blocks == null) return 0;
		else return blocks.size();
	}

}
