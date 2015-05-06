package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.Pose.PoseMethod;

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
	
	//Cache variables (can be updated)
	private Pose pose;						//3D superimposition information Pose
	private int length;						//total number of aligned positions, including gaps = sum of blocks lengths (cache)
	private int coreLength;					//number of aligned positions without gaps (cache)
	
	/**
	 * Constructor.
	 * @param multipleAlignment the parent MultipleAlignment of the BlockImpl instance.
	 * @return BlockSetImpl a BlockSetImpl instance linked to its parent MultipleAlignment.
	 * @throws StructureAlignmentException 
	 */
	public BlockSetImpl(MultipleAlignment multipleAlignment) throws StructureAlignmentException{
		
		parent = multipleAlignment;
		if (parent!=null) parent.getBlockSets().add(this);
		blocks = null;
		
		//Cache variables (can be updated)
		pose = new PoseBS(this);
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
				Block newB = (Block) b.clone();
				newB.setBlockSet(this); //This automatically adds the newB to the blocks list
			}
		}
	}
	
	@Override
	public Object clone(){
		return new BlockSetImpl(this);
	}
	
	@Override
	public String toString() {
		return "BlockSetImpl [parent=" + parent + ", blocks=" + blocks
				+ ", pose=" + pose + ", length=" + length + ", coreLength="
				+ coreLength + "]";
	}

	@Override
	public MultipleAlignment getMultipleAlignment() {
		return parent;
	}

	@Override
	public void setMultipleAlignment(MultipleAlignment parent) {
		//Delete the alignment instance from the parent list
		if (parent!=null) parent.getBlockSets().remove(this);
		this.parent = parent;
		//Cross-link parent and this instance
		if (parent!=null) parent.getBlockSets().add(this);
	}

	@Override
	public List<Block> getBlocks() {
		if (blocks==null) blocks = new ArrayList<Block>();
		return blocks;
	}

	@Override
	public void setBlocks(List<Block> blocks) {
		this.blocks = blocks;
	}

	@Override
	public Pose getPose() throws StructureAlignmentException {
		if (pose == null) pose = new PoseBS(this);
		return pose;
	}
	
	@Override
	public void updatePose(PoseMethod method) throws StructureException, StructureAlignmentException {
		if (pose == null) pose = new PoseBS(this);
		pose.updatePose(method);
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
		else if (getBlockNum()==0) throw new StructureAlignmentException("Empty BlockSet: parent==null && getBlockNum() == 0.");
		else return blocks.get(0).size();
	}

	@Override
	public int getBlockNum() throws StructureAlignmentException {
		if (blocks == null) throw new StructureAlignmentException("Empty BlockSet: getBlocks() == null.");
		else return blocks.size();
	}

	@Override
	public int getCoreLength() throws StructureAlignmentException {
		if (coreLength == -1) updateCoreLength();
		return coreLength;
	}

	@Override
	public void updateLength() throws StructureAlignmentException {
		if(getBlockNum()==0) throw new StructureAlignmentException("Empty BlockSet: getBlockNum() == 0.");
		//Try to calculate it from the Block information
		else {
			length = 0;
			for (Block block:blocks) length += block.length();
		}
	}

	@Override
	public void updateCoreLength() throws StructureAlignmentException {
		if(getBlockNum()==0) throw new StructureAlignmentException("Empty BlockSet: getBlockNum() == 0.");
		//Try to calculate it from the Block information
		else {
			coreLength = 0;
			for (Block block:blocks) coreLength += block.getCoreLength();
		}
	}

	@Override
	public void updateCache(PoseMethod method) throws StructureAlignmentException, StructureException {
		updatePose(method);
		updateCoreLength();
		updateLength();
	}
	
}
