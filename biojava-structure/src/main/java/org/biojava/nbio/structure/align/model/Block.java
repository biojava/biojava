package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * A Block is a Data Structure that stores aligned positions of a multiple alignment that fulfill the following conditions:
 * 
 * 		1- Residues are in a sequential order (increasing)
 * 		2- At least two structures have a residue in every column (no empty columns or columns with one residue = no alignment).
 *
 * It is part of a {@linkp BlockSet} instance, named as parent, which is in turn part of a {@link MultipleAlignment} instance.
 * 
 * @author Aleix Lafita
 * 
 */
public class Block implements Serializable, Cloneable{

	private static final long serialVersionUID = -5804042669466177641L;
	
	BlockSet parent;						//BlockSet instance
	List<List<Integer>> alignRes;			//residues aligned as a matrix of n*l size (n=nr.structures; l=block length)
	
	//Block Information - Utility
	int cols;							//number of aligned positions in the alignment, block length
	int rows;							//number of structures in the multiple alignment
	
	/**
	 * Constructor
	 */
	public Block(BlockSet blockSet) {
		
		parent = blockSet;
		alignRes = new ArrayList<List<Integer>>();
	}
	
	/**
	 * Copy constructor
	 */
	public Block(Block b) {
		
		this.parent = b.parent;
		this.alignRes = new ArrayList<List<Integer>>(b.alignRes);
		this.cols = b.cols;
		this.rows = b.rows;
	}
	
	/**
	 * Creates and returns a copy of this object. Uses the copy constructor.
	 */
	@Override
	public Block clone(){
		
		return new Block(this);
	}

	@Override
	public String toString() {
		return "Block [parent=" + parent + ", alignRes=" + alignRes + ", cols="
				+ cols + ", rows=" + rows + "]";
	}
	
	//Getters and Setters **************************************************************************************
	
	public BlockSet getParent() {
		return parent;
	}

	public void setParent(BlockSet parent) {
		this.parent = parent;
	}

	public List<List<Integer>> getAlignRes() {
		return alignRes;
	}

	public void setAlignRes(List<List<Integer>> alignRes) {
		this.alignRes = alignRes;
	}

	/**
	 * Number of aligned positions in the alignment, the block length.
	 */
	public int getCols() {
		return cols;
	}

	public void setCols(int cols) {
		this.cols = cols;
	}

	/**
	 * Number of structures in the multiple alignment, the block size.
	 */
	public int getRows() {
		return rows;
	}

	public void setRows(int rows) {
		this.rows = rows;
	}
}
