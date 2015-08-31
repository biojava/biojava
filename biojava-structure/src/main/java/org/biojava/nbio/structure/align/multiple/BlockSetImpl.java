package org.biojava.nbio.structure.align.multiple;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

/**
 * A general implementation of a BlockSet to store a flexible part of a 
 * multiple alignment.
 *
 * @author Aleix Lafita
 * @since 4.1.0
 * 
 */
public class BlockSetImpl extends AbstractScoresCache 
implements Serializable, BlockSet, Cloneable {

	private static final long serialVersionUID = -1015791986000076089L;

	//Member variables
	private MultipleAlignment parent;
	private List<Block> blocks;

	//Cache variables (can be updated)
	private List<Matrix4d> pose; //Transformation matrices
	private int length;
	private int coreLength;

	/**
	 * Constructor. Links also the parent to this instance by adding the
	 * BlockSet to the parent's List.
	 * 
	 * @param alignment MultipleAlignment parent of the BlockSet.
	 * @return BlockSet an instance linked to the parent alignment.
	 */
	public BlockSetImpl(MultipleAlignment alignment) {

		parent = alignment;
		if (parent!=null) parent.getBlockSets().add(this);
		blocks = null;

		pose = null;
		length = -1; //Value -1 reserved to indicate that has to be calculated
		coreLength = -1;
	}

	/**
	 * Copy constructor. Makes also a deep copy of all constituent 
	 * {@link Block}s.
	 * 
	 * @param bs BlockSet object to be copied.
	 * @return BlockSet an identical copy of the input object.
	 */
	public BlockSetImpl(BlockSetImpl bs){

		super(bs); //This copies the cached scores
		this.parent = bs.parent;
		this.length = bs.length;
		this.coreLength = bs.coreLength;

		this.pose = null;
		if (bs.pose != null){
			//Make a deep copy of everything
			this.pose = new ArrayList<Matrix4d>();
			for (Matrix4d trans:bs.pose){
				Matrix4d newTrans = (Matrix4d) trans.clone();
				pose.add(newTrans);
			}
		}

		blocks = null;
		if (bs.blocks!=null){
			//Make a deep copy of everything
			this.blocks = new ArrayList<Block>();
			for (Block b:bs.blocks){
				Block newB = b.clone();
				newB.setBlockSet(this);
				this.blocks.add(newB);
			}
		}
	}

	@Override
	public void clear() {
		super.clear();
		length = -1;
		coreLength = -1;
		pose = null;
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

	@Override
	public List<Matrix4d> getTransformations() {
		return pose;
	}

	@Override
	public void setTransformations(List<Matrix4d> transformations) {
		if(size() != transformations.size()) {
			throw new IllegalArgumentException(
					"Wrong number of structures for this alignment");
		}
		pose = transformations;
	}

	@Override
	public int length() {
		if (length == -1) updateLength();
		return length;
	}

	@Override
	public int size() {
		//Get the size from the variables that can contain the information
		if (parent != null) return parent.size();
		else if (getBlocks().size()==0) {
			throw new IndexOutOfBoundsException(
					"Empty BlockSet: number of Blocks == 0.");
		}
		else return blocks.get(0).size();
	}

	@Override
	public int getCoreLength() {
		if (coreLength == -1) updateCoreLength();
		return coreLength;
	}

	protected void updateLength() {
		if(getBlocks().size()==0) {
			throw new IndexOutOfBoundsException(
					"Empty BlockSet: number of Blocks == 0.");
		}
		//Try to calculate it from the Block information
		else {
			length = 0;
			for (Block block:blocks) length += block.length();
		}
	}

	protected void updateCoreLength() {
		if(getBlocks().size()==0) {
			throw new IndexOutOfBoundsException(
					"Empty BlockSet: number of Blocks == 0.");
		}
		//Try to calculate it from the Block information
		else {
			coreLength = 0;
			for (Block block:blocks) coreLength += block.getCoreLength();
		}
	}

	protected void updateCache() {
		updateCoreLength();
		updateLength();
	}
}
