package org.biojava.bio.structure;

import java.util.Iterator;

/**
 * Created by douglas on 1/23/15.
 */
public class AugmentedResidueRange extends ResidueRangeAndLength {

	private final AtomPositionMap map;

	public AugmentedResidueRange(String chain, ResidueNumber start, ResidueNumber end, int length, AtomPositionMap map) {
		super(chain, start, end, length);
		this.map = map;
	}

	public AugmentedResidueRange(String chain, String start, String end, int length, AtomPositionMap map) {
		super(chain, start, end, length);
		this.map = map;
	}


	/**
	 * Returns the ResidueNumber that is at position {@code positionInRange} in <em>this</em> ResidueRange.
	 * @return The ResidueNumber, or false if it does not exist or is not within this ResidueRange
	 */
	public ResidueNumber getResidue(int positionInRange) {
		return super.getResidue(positionInRange, map);
	}

	/**
	 * @return True if and only if {@code residueNumber} is within this ResidueRange
	 */
	public boolean contains(ResidueNumber residueNumber) {
		return super.contains(residueNumber, map);
	}

	/**
	 * Returns a new Iterator over every {@link ResidueNumber} in this ResidueRange.
	 * Stores the contents of {@code map} until the iterator is finished, so calling code should set the iterator to {@code null} if it did not finish.
	 */
	public Iterator<ResidueNumber> iterator() {
		return super.iterator(map);
	}

}
