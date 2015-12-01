/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure.align.multiple;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;

/**
 * A general implementation of a {@link MultipleAlignment}.
 *
 * @author Aleix Lafita
 * @since 4.1.0
 * 
 */
public class MultipleAlignmentImpl extends AbstractScoresCache 
implements Serializable, MultipleAlignment, Cloneable {

	private static final long serialVersionUID = 3432043794125805139L;

	private MultipleAlignmentEnsemble parent;
	private List<BlockSet> blockSets;

	//Cache variables (can be updated)
	private int length;
	private int coreLength;
	private List<Matrix4d> pose;

	/**
	 * Default Constructor. Empty alignment. No structures assigned.
	 * 
	 * @return MultipleAlignment an empty MultipleAlignment instance.
	 */
	public MultipleAlignmentImpl() {
		this(new MultipleAlignmentEnsembleImpl());  //assign an empty ensemble.
	}

	/**
	 * Constructor linking to an existing ensemble.
	 * Automatically adds this alignment to the parent ensemble.
	 * 
	 * @param ensemble parent MultipleAlignmentEnsemble.
	 * @return MultipleAlignment an alignment instance part of an ensemble.
	 */
	public MultipleAlignmentImpl(MultipleAlignmentEnsemble ensemble) {

		super();
		parent = ensemble;
		if (parent!=null) parent.getMultipleAlignments().add(this);

		blockSets = null;
		pose = null;

		length = -1; //Value -1 reserved to indicate that has to be calculated
		coreLength = -1;
	}

	/**
	 * Copy constructor. Recursively copies member BlockSets.
	 * 
	 * @param ma MultipleAlignmentImpl to copy.
	 * @return MultipleAlignmentImpl identical copy of the alignment.
	 */
	public MultipleAlignmentImpl(MultipleAlignmentImpl ma) {

		super(ma); //Copy the scores
		parent = ma.parent;

		pose = null;
		if (ma.pose != null){
			//Make a deep copy of everything
			this.pose = new ArrayList<Matrix4d>();
			for (Matrix4d trans:ma.pose){
				Matrix4d newTrans = (Matrix4d) trans.clone();
				pose.add(newTrans);
			}
		}

		length = ma.length;
		coreLength = ma.coreLength;

		blockSets = null;
		if (ma.blockSets!=null){
			//Make a deep copy of everything
			this.blockSets = new ArrayList<BlockSet>();
			for (BlockSet bs:ma.blockSets){
				BlockSet newBS = bs.clone();
				newBS.setMultipleAlignment(this);
				this.blockSets.add(newBS);
			}
		}
	}

	@Override
	public void clear() {
		super.clear();
		length = -1;
		coreLength = -1;
		pose = null;
		for(BlockSet a : getBlockSets()) {
			a.clear();
		}
	}

	@Override
	public MultipleAlignmentImpl clone() {
		return new MultipleAlignmentImpl(this);
	}

	@Override
	public String toString() {
		String resume = "Structures:" + parent.getStructureNames() + 
				" \nAlgorithm:" + parent.getAlgorithmName() + "_" + 
				parent.getVersion() + 
				" \nBlockSets: "+ getBlockSets().size() + 
				" \nBlocks: " + getBlocks().size() +
				" \nLength: " + length() +
				" \nCore Length: "+ getCoreLength();
		for (String score:getScores()){
			resume += " \n"+score+": ";
			resume += String.format("%.2f", getScore(score));
		}
		return resume;
	}

	@Override
	public List<BlockSet> getBlockSets() {
		if (blockSets == null) blockSets = new ArrayList<BlockSet>();
		return blockSets;
	}

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
	public BlockSet getBlockSet(int index){
		return blockSets.get(index);
	}

	@Override
	public Block getBlock(int index){
		List<Block> blocks = getBlocks();
		return blocks.get(index);
	}

	@Override
	public List<Atom[]> getAtomArrays() {
		return parent.getAtomArrays();
	}

	@Override
	public int size() {
		return parent.size();
	}

	@Override
	public int length() {
		if (length < 0 ) updateLength();
		return length;
	}

	@Override
	public int getCoreLength() {
		if (coreLength < 0) updateCoreLength();
		return coreLength;
	}

	/**
	 * Force recalculation of the length (aligned columns) based on the 
	 * BlockSet lengths.
	 */
	protected void updateLength() {
		if(getBlockSets().size()==0) {
			throw new IndexOutOfBoundsException(
					"Empty MultipleAlignment: blockSets size == 0.");
		} //Otherwise try to calculate it from the BlockSet information
		else {
			length = 0;
			for (BlockSet blockSet:blockSets) length += blockSet.length();
		}
	}

	/**
	 * Force recalculation of the core length (ungapped columns) based on the 
	 * BlockSet core lengths.
	 */
	protected void updateCoreLength() {
		if(getBlockSets().size()==0) {
			throw new IndexOutOfBoundsException(
					"Empty MultipleAlignment: blockSets size == 0.");
		} //Otherwise try to calculate it from the BlockSet information
		else {
			coreLength = 0;
			for (BlockSet blockSet:blockSets) 
				coreLength += blockSet.getCoreLength();
		}
	}

	/**
	 * Updates all cached properties
	 * 
	 * @throws StructureException
	 */
	protected void updateCache() {
		updateCoreLength();
		updateLength();
	}

	@Override
	public MultipleAlignmentEnsemble getEnsemble() {
		return parent;
	}

	@Override
	public void setEnsemble(MultipleAlignmentEnsemble parent) {
		this.parent = parent;
	}
}
