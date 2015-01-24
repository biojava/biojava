package org.biojava.bio.structure;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * A chain, a start residue, and an end residue.
 *
 * Also stores a length. Because of insertion codes, this length is not necessarily {@code end − start}.
 */
public class ResidueRangeAndLength extends ResidueRange {

	private final int length;

	public ResidueRangeAndLength(String chain, ResidueNumber start, ResidueNumber end, int length) {
		super(chain, start, end);
		this.length = length;
	}

	public ResidueRangeAndLength(String chain, String start, String end, int length) {
		super(chain, start, end);
		this.length = length;
	}

	/**
	 * Returns a new Iterator over every {@link ResidueNumber} in this ResidueRange.
	 * Stores the contents of {@code map} until the iterator is finished, so calling code should set the iterator to {@code null} if it did not finish.
	 */
	@Override
	public Iterator<ResidueNumber> iterator(AtomPositionMap map) {
		return super.iterator(map, length); // just a bit faster
	}

	/**
	 * Returns a new Iterator over every {@link ResidueNumber} in this ResidueRange.
	 * Stores the contents of {@code map} until the iterator is finished, so calling code should set the iterator to {@code null} if it did not finish.
	 */
	@Override
	public Iterator<ResidueNumber> iterator(AtomPositionMap map, Integer length) {
		return super.iterator(map, length); // just a bit faster
	}

	/**
	 * Calculates the combined number of residues of the ResidueRanges in {@code rrs},
	 * <em>given that each ResidueRange has a length calculated</em>. The value, if calculated,
	 * <em>will include any alignment gaps</em>.
	 *
	 * @param rrs
	 *            A list of ResidueRanges
	 * @return The combined length
	 * @throws IllegalArgumentException
	 *             If the {@link #getLength() length} of one or more ResidueRange is null
	 * @see #getLength()
	 */
	public static int calcLength(List<ResidueRangeAndLength> rrs) {
		int l = 0;
		for (ResidueRangeAndLength rr : rrs) {
			l += rr.getLength();
		}
		return l;
	}

	/**
	 * @param s
	 *            A string of the form chain_start-end. For example: <code>A.5-100</code>.
	 * @return The unique ResidueRange corresponding to {@code s}.
	 */
	public static ResidueRangeAndLength parse(String s, AtomPositionMap map) {
		ResidueRange rr = parse(s);
		ResidueNumber start = rr.getStart();

		// get a non-null start and end
		// if it's the whole chain, choose the first and last residue numbers in the chain
		if (start==null) {
			if (rr.getChainId() == null) {
				start = map.getFirst();
			} else {
				start = map.getFirst(rr.getChainId());
			}
		}
		ResidueNumber end = rr.getEnd();
		if (end==null) { // should happen iff start==null
			if (rr.getChainId() == null) {
				end = map.getLast();
			} else {
				end = map.getLast(rr.getChainId());
			}
		}

		// now use those to calculate the length
		// if start or end is null, will throw NPE
		int length = map.calcLength(start, end);

		// to avoid confusing the user, don't add the start and end if they weren't specified in the string s
		return new ResidueRangeAndLength(rr.getChainId(), rr.getStart(), rr.getEnd(), length);
	}

	public static List<ResidueRangeAndLength> parseMultiple(List<String> ranges, AtomPositionMap map) {
		List<ResidueRangeAndLength> rrs = new ArrayList<ResidueRangeAndLength>(ranges.size());
		for (String range : ranges) {
			ResidueRangeAndLength rr = ResidueRangeAndLength.parse(range, map);
			if (rr != null) rrs.add(rr);
		}
		return rrs;
	}

	/**
	 * @param s
	 *            A string of the form chain_start-end,chain_start-end, ... For example:
	 *            <code>A.5-100,R_110-190,Z_200-250</code>.
	 * @return The unique ResidueRange corresponding to {@code s}.
	 */
	public static List<ResidueRangeAndLength> parseMultiple(String s, AtomPositionMap map) {
		String[] parts = s.split(",");
		List<ResidueRangeAndLength> list = new ArrayList<ResidueRangeAndLength>(parts.length);
		for (String part : parts) {
			list.add(parse(part, map));
		}
		return list;
	}

	/**
	 * @return The number of residues in this ResidueRange, including any alignment gaps. This value will be null if and
	 *         only if this ResidueRange was created with a null length.
	 */
	public int getLength() {
		return length;
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) {
			return true;
		}
		if (o == null || getClass() != o.getClass()) {
			return false;
		}
		if (!super.equals(o)) {
			return false;
		}
		ResidueRangeAndLength that = (ResidueRangeAndLength) o;
		return length == that.length;
	}

	@Override
	public int hashCode() {
		int result = super.hashCode();
		result = 31 * result + length;
		return result;
	}
}
