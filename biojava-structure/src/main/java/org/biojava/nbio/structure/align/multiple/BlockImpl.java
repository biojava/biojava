package org.biojava.nbio.structure.align.multiple;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * General implementation of a Block that supports alignments with gaps.
 * 
 * @author Aleix Lafita
 * 
 */
public class BlockImpl extends AbstractScoresCache implements Serializable, Block, Cloneable{

	private static final long serialVersionUID = -5804042669466177641L;
	
	private BlockSet parent;						//BlockSet instance
	private List<List<Integer>> alignRes;			//residues aligned as a double list of length=n and size=l (n=nr.structures; l=block length)
	private int coreLength;							//number of residues aligned without gaps (cache)
	
	/**
	 * Constructor.
	 * @param blockSet the parent BlockSet of the BlockImpl instance.
	 * @return BlockImpl a BlockImpl instance linked to its parent BlockSet.
	 */
	public BlockImpl(BlockSet blockSet) {
		
		parent = blockSet;
		if (parent!=null) parent.getBlocks().add(this);
		
		alignRes = null;
		coreLength = -1;							//Value -1 indicates not calculated.
	}
	
	
	/**
	 * Copy constructor.
	 * @param b BlockImpl object to be copied.
	 * @return BlockImpl an identical copy of the input BlockImpl object.
	 */
	public BlockImpl(BlockImpl b) {
		
		this.parent = b.parent;
		this.coreLength = b.coreLength;
		
		this.alignRes = null;
		if (b.alignRes!=null){
			//Make a deep copy of everything
			alignRes = new ArrayList<List<Integer>>();
			for (int k=0; k<b.size(); k++)
				alignRes.add(new ArrayList<Integer>(b.alignRes.get(k)));
		}
	}
	
	@Override
	public Block clone(){
		return new BlockImpl(this);
	}
	
	
	/**
	 * Clear scores and cached properties. Recursively clears member blocks.
	 */
	@Override
	public void clear() {
		super.clear();
		coreLength = -1;
	}
	
	@Override
	public String toString() {
		return "BlockImpl [alignRes=" + alignRes
				+ ", coreLength=" + coreLength + "]";
	}

	@Override
	public void setBlockSet(BlockSet parent) {
		this.parent = parent;
	}

	@Override
	public BlockSet getBlockSet() {
		return parent;
	}

	@Override
	public List<List<Integer>> getAlignRes() {
		return alignRes;
	}

	@Override
	public void setAlignRes(List<List<Integer>> alignRes) {
		this.alignRes = alignRes;
	}

	@Override
	public int length() {
		if (alignRes == null) return 0;
		if (alignRes.size() == 0) return 0;
		return alignRes.get(0).size();
	}

	@Override
	public int size() {
		return alignRes.size();
	}

	@Override
	public int getCoreLength() {
		if(coreLength == -1) updateCoreLength();
		return coreLength;
	}

	protected void updateCoreLength() {
		coreLength = 0;
		//Loop through all columns of the alignment and count how many of them do not have gaps in any structure
		for (int col=0; col<length(); col++){
			boolean core = true;
			for (int str=0; str<size(); str++){
				if (alignRes.get(str).get(col) == null){
					core = false;
					break;
				}
			}
			if (core) coreLength++;
		}
	}
}
