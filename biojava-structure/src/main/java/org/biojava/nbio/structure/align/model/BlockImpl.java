package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * General implementation of a Block that supports alignments with gaps.
 * 
 * @author Aleix Lafita
 * 
 */
public class BlockImpl implements Serializable, Block{

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
	public Object clone(){
		return new BlockImpl(this);
	}

	@Override
	public String toString() {
		return "BlockImpl [parent=" + parent + ", alignRes=" + alignRes
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
		if (alignRes == null) alignRes = new ArrayList<List<Integer>>();
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

	@Override
	public int updateCoreLength() {
		// TODO Auto-generated method stub
		//Loop through all the columns of the alignments and count how many of them do not have gaps.
		return 0;
	}
	
}
