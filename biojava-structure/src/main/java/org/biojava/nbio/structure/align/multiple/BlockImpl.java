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

import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentTools;

/**
 * General implementation of a {@link Block} that supports any type of
 * sequential alignment with gaps.
 *
 * @author Aleix Lafita
 * @since 4.1.0
 *
 */
public class BlockImpl extends AbstractScoresCache implements Serializable,
		Block, Cloneable {

	private static final long serialVersionUID = -5804042669466177641L;

	private BlockSet parent;
	private List<List<Integer>> alignRes;

	// CACHE
	private int coreLength;
	private List<Integer> alignResCounts;

	/**
	 * Constructor. Links also the parent to this instance, by adding the Block
	 * to the parent's list.
	 *
	 * @param blockSet
	 *            the parent BlockSet of the BlockImpl instance.
	 * @return BlockImpl a BlockImpl instance linked to its parent BlockSet.
	 */
	public BlockImpl(BlockSet blockSet) {

		parent = blockSet;
		if (parent != null)
			parent.getBlocks().add(this);
		alignRes = null;
		coreLength = -1; // Value -1 indicates not yet calculated.
		alignResCounts = null; // Value null not yet calculated
	}

	/**
	 * Copy constructor.
	 *
	 * @param b
	 *            BlockImpl object to be copied.
	 * @return BlockImpl an identical copy of the input BlockImpl object.
	 */
	public BlockImpl(BlockImpl b) {

		super(b); // This copies the cached scores
		this.parent = b.parent;
		this.coreLength = b.coreLength;

		this.alignRes = null;
		if (b.alignRes != null) {
			// Make a deep copy of everything
			alignRes = new ArrayList<List<Integer>>();
			for (int k = 0; k < b.size(); k++) {
				alignRes.add(new ArrayList<Integer>(b.alignRes.get(k)));
			}
		}
	}

	@Override
	public Block clone() {
		return new BlockImpl(this);
	}

	@Override
	public void clear() {
		super.clear();
		coreLength = -1;
		alignResCounts = null;
	}

	@Override
	public String toString() {
		return "BlockImpl [alignRes=" + alignRes + ", coreLength=" + coreLength
				+ "]";
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
		if (alignRes == null)
			return 0;
		if (alignRes.size() == 0)
			return 0;
		return alignRes.get(0).size();
	}

	@Override
	public int size() {
		return alignRes.size();
	}

	@Override
	public int getCoreLength() {
		if (coreLength == -1)
			updateCoreLength();
		return coreLength;
	}

	protected void updateCoreLength() {
		coreLength = MultipleAlignmentTools.getCorePositions(this).size();
	}

	@Override
	public int getStartIndex(int str) {
		int index = -1;
		Integer start = null;
		while (start == null && index < alignRes.get(str).size()) {
			index++;
			start = alignRes.get(str).get(index);
		}
		return index;
	}

	@Override
	public int getStartResidue(int str) {
		return alignRes.get(str).get(getStartIndex(str));
	}

	@Override
	public int getFinalIndex(int str) {
		int index = alignRes.get(str).size();
		Integer end = null;
		while (end == null && index >= 0) {
			index--;
			end = alignRes.get(str).get(index);
		}
		return index;
	}

	@Override
	public int getFinalResidue(int str) {
		return alignRes.get(str).get(getFinalIndex(str));
	}

	@Override
	public List<Integer> getAlignResCounts() {

		if (alignResCounts != null)
			return alignResCounts;

		alignResCounts = new ArrayList<Integer>(size());
		for (int s = 0; s < size(); s++) {
			int count = 0;
			for (int r = 0; r < length(); r++) {
				if (alignRes.get(s).get(r) != null)
					count++;
			}
			alignResCounts.add(count);
		}
		return alignResCounts;
	}

}
