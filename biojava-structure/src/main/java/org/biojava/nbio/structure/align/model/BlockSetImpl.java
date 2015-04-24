package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.Pose.PoseMethod;
import org.biojava.nbio.structure.jama.Matrix;

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
	private double coverage;				//coverage of the alignment
	private double similarity;				//similarity measure for the aligned structures
	
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
		coverage = 0;
		similarity = 0;
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

	@Override
	public double getSimilarity() {
		return similarity;
	}

	@Override
	public void updateSimilarity() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double getCoverage() {
		return coverage;
	}

	@Override
	public void updateCoverage() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void updatePose(PoseMethod method) throws StructureException{
		
		//TODO is it better to catch the errors or to throw them?
		//TODO should throw error? maybe we should not check for that because a parent and an alignment should be there.
		if (parent == null) return;		//If there is not a parent the atoms cannot be obtained so nothing is done
		if (blocks == null) return;		//If the aligned Blocks are not initialized do nothing
		
		//Create a new Pose if there is not one
		if (pose == null) pose = new PoseImpl(this);
		
		//Initialize or replace the rotation and translation variables
		pose.setRotation(new ArrayList<Matrix>());
		pose.setTranslation(new ArrayList<Atom>());
		
		switch (method) {
		case REFERENCE:
			//We suppose the first molecule as reference and superimpose everything to it
			for (int i=0; i<size(); i++){
				List<Atom> atomSet1 = new ArrayList<Atom>();
				List<Atom> atomSet2 = new ArrayList<Atom>();
				for (int k=0; k<getBlockNum(); k++){
					for (int j=0; j<blocks.get(k).length(); j++){
						Integer pos1 = blocks.get(k).getAlignRes().get(0).get(j);
						Integer pos2 = blocks.get(k).getAlignRes().get(i).get(j);
						
						if (pos1==null || pos2==null) continue;
						atomSet1.add((Atom) parent.getAtomArrays().get(0)[pos1].clone());
						atomSet2.add((Atom) parent.getAtomArrays().get(i)[pos2].clone());
					}
				}
				SVDSuperimposer svd = new SVDSuperimposer(atomSet1.toArray(new Atom[0]), atomSet2.toArray(new Atom[0]));
				pose.getRotation().add(svd.getRotation());
				pose.getTranslation().add(svd.getTranslation());
			}
			break;
		case MEDIAN:
			//TODO implement this type of superimposition
			System.out.println("Not yet implemented!");
			break;
		case CONSENSUS:
			//TODO implement this type of superimposition
			System.out.println("Not yet implemented!");
			break;
		}
		
		pose.updateRMSD();
		pose.updateTMscore();
	}

}
