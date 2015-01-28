package org.biojava.bio.structure;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A chain, a start residue, and an end residue.
 *
 * Also stores a length. Because of insertion codes, this length is not necessarily {@code end − start}.
 */
public class ResidueRangeAndLength extends ResidueRange {
	private static final Logger logger = LoggerFactory.getLogger(ResidueRangeAndLength.class);

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
		return super.iterator(map); // just a bit faster
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
	 * Parses a residue range.
	 * 
	 * The AtomPositionMap is used to calculate the length and fill in missing
	 * information, such as for whole chains ('A:'). Supports the special chain
	 * name '_' for single-chain structures.
	 * @param s
	 *            A string of the form chain_start-end. For example: <code>A.5-100</code>.
	 * @return The unique ResidueRange corresponding to {@code s}.
	 */
	public static ResidueRangeAndLength parse(String s, AtomPositionMap map) {
		ResidueRange rr = parse(s);
		ResidueNumber start = rr.getStart();
		
		String chain = rr.getChainId();
		
		// handle special "_" chain
		if(chain == null || chain.equals("_")) {
			ResidueNumber first = map.getNavMap().firstKey();
			chain = first.getChainId();
			// Quick check for additional chains. Not guaranteed.
			if( ! map.getNavMap().lastKey().getChainId().equals(chain) ) {
				logger.warn("Multiple possible chains match '_'. Using chain {}",chain);
			}
		}

		// get a non-null start and end
		// if it's the whole chain, choose the first and last residue numbers in the chain
		if (start==null) {
			start = map.getFirst(chain);
		}
		ResidueNumber end = rr.getEnd();
		if (end==null) { // should happen iff start==null
			end = map.getLast(chain);
		}

		// Replace '_'
		start.setChainId(chain);
		end.setChainId(chain);
		
		// now use those to calculate the length
		// if start or end is null, will throw NPE
		int length = map.calcLength(start, end);

		return new ResidueRangeAndLength(chain, start, end, length);
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
