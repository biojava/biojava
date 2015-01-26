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
 * Created on 2012-11-20
 *
 */

package org.biojava.bio.structure;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A chain, a start residue, and an end residue.
 * 
 * @author dmyerstu
 * @see ResidueNumber
 * @see org.biojava.bio.structure.ResidueRangeAndLength
 */
public class ResidueRange {

	private final String chain;
	private final ResidueNumber end;
	private final ResidueNumber start;

	public static final Pattern RANGE_REGEX = Pattern.compile(
			"^\\s*([a-zA-Z0-9]+|_)" + //chain ID. Be flexible here, rather than restricting to 4-char IDs
			"(?:" + //begin range, this is a "non-capturing group"
				"(?::|_|:$|_$|$)" + //colon or underscore, could be at the end of a line, another non-capt. group.
				"(?:"+ // another non capturing group for the residue range
					"([-+]?[0-9]+[A-Za-z]?)" + // first residue
					"\\s*-\\s*" + // -
					"([-+]?[0-9]+[A-Za-z]?)" + // second residue
				")?+"+
			")?" + //end range
			"\\s*");

	/**
	 * @param s
	 *            A string of the form chain_start-end or chain.start-end. For example: <code>A.5-100</code> or <code>A_5-100</code>.
	 * @return The unique ResidueRange corresponding to {@code s}
	 */
	public static ResidueRange parse(String s) {
		ResidueNumber start = null, end = null;
		String chain = null;
		Matcher matcher = RANGE_REGEX.matcher(s);
		if (matcher.matches()) {
			try {
				chain = matcher.group(1);
				if (matcher.group(2) != null) {
					start = ResidueNumber.fromString(matcher.group(2));
					end = ResidueNumber.fromString(matcher.group(3));
					start.setChainId(chain);
					end.setChainId(chain);
				}
			} catch (IllegalStateException e) {
				throw new IllegalArgumentException("Range " + s + " was not valid", e);
			}
		}
		return new ResidueRange(chain, start, end);
	}

	/**
	 * @param s
	 *            A string of the form chain_start-end,chain_start-end, ... For example:
	 *            <code>A.5-100,R_110-190,Z_200-250</code>.
	 * @return The unique ResidueRange corresponding to {@code s}.
	 */
	public static List<ResidueRange> parseMultiple(String s) {
		String[] parts = s.split(",");
		List<ResidueRange> list = new ArrayList<ResidueRange>(parts.length);
		for (String part : parts) {
			list.add(parse(part));
		}
		return list;
	}

	public ResidueRange(String chain, String start, String end) {
		this.chain = chain;
		this.start = ResidueNumber.fromString(start);
		this.start.setChainId(chain);
		this.end = ResidueNumber.fromString(end);
		this.end.setChainId(chain);
	}

	public ResidueRange(String chain, ResidueNumber start, ResidueNumber end) {
		this.chain = chain;
		this.start = start;
		this.end = end;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		ResidueRange other = (ResidueRange) obj;
		if (chain == null) {
			if (other.chain != null) return false;
		} else if (!chain.equals(other.chain)) return false;
		if (end == null) {
			if (other.end != null) return false;
		} else if (!end.equals(other.end)) return false;
		if (start == null) {
			if (other.start != null) return false;
		} else if (!start.equals(other.start)) return false;
		return true;
	}

	public String getChainId() {
		return chain;
	}

	public ResidueNumber getEnd() {
		return end;
	}

	public ResidueNumber getStart() {
		return start;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (chain == null ? 0 : chain.hashCode());
		result = prime * result + (end == null ? 0 : end.hashCode());
		result = prime * result + (start == null ? 0 : start.hashCode());
		return result;
	}

	@Override
	public String toString() {
		return chain + "_" + start + "-" + end;
	}

	/**
	 * Returns the ResidueNumber that is at position {@code positionInRange} in <em>this</em> ResidueRange.
	 * @return The ResidueNumber, or false if it does not exist or is not within this ResidueRange
	 */
	public ResidueNumber getResidue(int positionInRange, AtomPositionMap map) {
		if (map == null) throw new IllegalArgumentException("The AtomPositionMap must be non-null");
		int i = 0;
		for (Map.Entry<ResidueNumber, Integer> entry : map.getNavMap().entrySet()) {
			if (i == positionInRange) return entry.getKey();
			if (contains(entry.getKey(), map)) {
				i++;
			}
		}
		return null;
	}

	/**
	 * @return True if and only if {@code residueNumber} is within this ResidueRange
	 */
	public boolean contains(ResidueNumber residueNumber, AtomPositionMap map) {
		if (residueNumber == null) throw new IllegalArgumentException("Can't find a null ResidueNumber");
		if (map == null) throw new IllegalArgumentException("The AtomPositionMap must be non-null");
		if (getStart() == null || getEnd() == null) throw new IllegalArgumentException("The bounds of this ResidueNumber aren't known");
		Integer pos = map.getPosition(residueNumber);
		if (pos == null) throw new IllegalArgumentException("Couldn't find residue " + residueNumber.printFull());
		Integer startPos = map.getPosition(getStart());
		if (startPos == null) throw new IllegalArgumentException("Couldn't find the start position");
		Integer endPos = map.getPosition(getEnd());
		if (endPos == null) throw new IllegalArgumentException("Couldn't find the end position");
		return pos >= startPos && pos <= endPos;
	}

	/**
	 * Returns a new Iterator over every {@link ResidueNumber} in this ResidueRange.
	 * Stores the contents of {@code map} until the iterator is finished, so calling code should set the iterator to {@code null} if it did not finish.
	 */
	public Iterator<ResidueNumber> iterator(final AtomPositionMap map) {
		return iterator(map, null);
	}

	/**
	 * Returns a new Iterator over every {@link ResidueNumber} in this ResidueRange.
	 * Stores the contents of {@code map} until the iterator is finished, so calling code should set the iterator to {@code null} if it did not finish.
	 */
	public Iterator<ResidueNumber> iterator(final AtomPositionMap map, Integer length) {
		// get the length without the side effect of setting it
		final int theLength = length==null? map.calcLength(start, end) : length;
		return new Iterator<ResidueNumber>() {
			private ResidueNumber[] residueNumbers = new ResidueNumber[map.getNavMap().size()];
			private int i = -1;
			@Override
			public boolean hasNext() {
				return i < theLength;
			}
			@Override
			public ResidueNumber next() {
				if (i == -1) {
					residueNumbers = new ResidueNumber[map.getNavMap().size()];
					int j = 0;
					for (Map.Entry<ResidueNumber,Integer> entry : map.getNavMap().entrySet()) {
						residueNumbers[j] = entry.getKey();
						if (contains(entry.getKey(), map)) {
							j++;
						}
					}
				}
				i++;
				ResidueNumber rn = residueNumbers[i];
				// let's assume we're not going to use this anymore
				// destroy array to free memory
				// we can always reconstruct
				if (i > theLength) {
					residueNumbers = null;
					i = -1;
				}
				return rn;
			}
			@Override
			public void remove() {
				throw new UnsupportedOperationException("Not modifiable");
			}
		};
	}

	/**
	 * Returns a new Iterator over every {@link ResidueNumber} in the list of ResidueRanges.
	 * Stores the contents of {@code map} until the iterator is finished, so calling code should set the iterator to {@code null} if it did not finish.
	 */
	public static Iterator<ResidueNumber> multiIterator(final AtomPositionMap map, final ResidueRange... rrs) {
		return new Iterator<ResidueNumber>() {
			private int r = 0;
			private Iterator<ResidueNumber> internal;
			@Override
			public boolean hasNext() {
				if (r == rrs.length - 1) {
					init();
					return internal.hasNext();
				}
				return true;
			}
			private void init() {
				if (internal == null) {
					internal = rrs[r].iterator(map);
				}
			}
			@Override
			public ResidueNumber next() {
				if (rrs.length == 0) throw new NoSuchElementException();
				init();
				if (!hasNext()) throw new NoSuchElementException();
				if (!internal.hasNext()) {
					r++;
					internal = rrs[r].iterator(map);
				}
				return internal.next();
			}
			@Override
			public void remove() {
				throw new UnsupportedOperationException("Not modifiable");
			}
		};
	}

	/**
	 * Returns a new Iterator over every {@link ResidueNumber} in the list of ResidueRanges.
	 * Stores the contents of {@code map} until the iterator is finished, so calling code should set the iterator to {@code null} if it did not finish.
	 */
	public static Iterator<ResidueNumber> multiIterator(AtomPositionMap map, List<? extends ResidueRange> rrs) {
		ResidueRange[] ranges = new ResidueRange[rrs.size()];
		for (int i = 0; i < rrs.size(); i++) {
			ranges[i] = rrs.get(i);
		}
		return multiIterator(map, ranges);
	}

	public static List<ResidueRange> parseMultiple(List<String> ranges) {
		List<ResidueRange> rrs = new ArrayList<ResidueRange>(ranges.size());
		for (String range : ranges) {
			ResidueRange rr = ResidueRange.parse(range);
			if (rr != null) rrs.add(rr);
		}
		return rrs;
	}

	public static List<String> toStrings(List<? extends ResidueRange> ranges) {
		List<String> list = new ArrayList<String>(ranges.size());
		for (ResidueRange range : ranges) {
			list.add(range.toString());
		}
		return list;
	}
	
	public static String toString(List<? extends ResidueRange> ranges) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < ranges.size(); i++) {
			sb.append(ranges.get(i));
			if (i < ranges.size() - 1) sb.append(",");
		}
		return sb.toString();
	}

}
