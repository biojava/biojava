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
 * A chain, a start residue, and an end residue. May also store a length value. Because of insertion codes, this length
 * is not necessarily {@code end − start}.
 * 
 * @author dmyerstu
 * @see ResidueNumber
 */
public class ResidueRange {

	private final String chain;
	private final ResidueNumber end;
	private final Integer length;
	private final ResidueNumber start;

	public static final String RANGE_REGEX = "^([a-zA-Z])+[_:](?:(-?\\d+[a-zA-Z]?)-(-?\\d+[a-zA-Z]?))?$";

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
	public static int calcLength(List<ResidueRange> rrs) {
		int l = 0;
		for (ResidueRange rr : rrs) {
			if (rr.getLength() == null) throw new IllegalArgumentException(
					"At least one ResidueRange does not have a length.");
			l += rr.getLength();
		}
		return l;
	}

	/**
	 * @param s
	 *            A string of the form chain_start-end or chain.start-end. For example: <code>A.5-100</code> or <code>A_5-100</code>.
	 * @return The unique ResidueRange corresponding to {@code s}
	 */
	public static ResidueRange parse(String s) {
		ResidueNumber start = null, end = null;
		String chain = null;
		Pattern pattern = Pattern.compile(RANGE_REGEX);
		Matcher matcher = pattern.matcher(s);
		matcher.find();
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
		return new ResidueRange(chain, start, end, null);
	}

	/**
	 * @param s
	 *            A string of the form chain_start-end. For example: <code>A.5-100</code>.
	 * @return The unique ResidueRange corresponding to {@code s}.
	 */
	public static ResidueRange parse(String s, AtomPositionMap map) {
		ResidueRange rr = parse(s);
		if (rr.getStart() == null) { // whole chain
			String chain = rr.getChainId();
			rr = new ResidueRange(chain, map.getFirst(chain), map.getLast(chain), null);
		}
		int length = map.calcLength(rr.getStart(), rr.getEnd());
		return new ResidueRange(rr.getChainId(), rr.getStart(), rr.getEnd(), length);
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

	/**
	 * @param s
	 *            A string of the form chain_start-end,chain_start-end, ... For example:
	 *            <code>A.5-100,R_110-190,Z_200-250</code>.
	 * @return The unique ResidueRange corresponding to {@code s}.
	 */
	public static List<ResidueRange> parseMultiple(String s, AtomPositionMap map) {
		String[] parts = s.split(",");
		List<ResidueRange> list = new ArrayList<ResidueRange>(parts.length);
		for (String part : parts) {
			list.add(parse(part, map));
		}
		return list;
	}

	public ResidueRange(char chain, ResidueNumber start, ResidueNumber end, Integer length) {
		this(String.valueOf(chain), start, end, length);
	}

	public ResidueRange(String chain, ResidueNumber start, ResidueNumber end, Integer length) {
		this.chain = chain;
		this.start = start;
		this.end = end;
		this.length = length;
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

	/**
	 * Returns the chain Id as a char.
	 * @deprecated Use {@link #getChainId()} instead, which does not require that chain Ids have a length of 1
	 * @throws IllegalArgumentException If the chain Id contains more than 1 character
	 */
	@Deprecated
	public char getChain() {
		if (chain.length() > 1) throw new IllegalArgumentException("Can't return full chain Id " + chain);
		return chain.charAt(0);
	}

	public String getChainId() {
		return chain;
	}

	public ResidueNumber getEnd() {
		return end;
	}

	/**
	 * @return The number of residues in this ResidueRange, including any alignment gaps. This value will be null if and
	 *         only if this ResidueRange was created with a null length.
	 */
	public Integer getLength() {
		return length;
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
	 * @return True if and only if {@code residueNumber} is within this ResidueRange
	 */
	public boolean contains(ResidueNumber residueNumber, AtomPositionMap map) {
		int pos = map.getPosition(residueNumber);
		int startPos = map.getPosition(start);
		int endPos = map.getPosition(end);
		return pos >= startPos && pos <= endPos;
	}

	/**
	 * Returns the ResidueNumber that is at position {@code positionInRange} in <em>this</em> ResidueRange.
	 * @return The ResidueNumber, or false if it does not exist or is not within this ResidueRange
	 */
	public ResidueNumber getResidue(int positionInRange, AtomPositionMap map) {
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
	 * Returns a new Iterator over every {@link ResidueNumber} in this ResidueRange.
	 * Stores the contents of {@code map} until the iterator is finished, so calling code should set the iterator to {@code null} if it did not finish.
	 */
	public Iterator<ResidueNumber> iterator(final AtomPositionMap map) {
		// get the length without the side effect of setting it
		int theLength = 0;
		if (length == null) {
			theLength = map.calcLength(start, end);
		} else {
			theLength = this.length;
		}
		final int length = theLength;
		return new Iterator<ResidueNumber>() {
			private ResidueNumber[] residueNumbers = new ResidueNumber[map.getNavMap().size()];
			private int i = -1;
			@Override
			public boolean hasNext() {
				return i < length;
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
				if (i > length) {
					residueNumbers = null;
					i = -1;
				}
				return rn;
			}
			@Override
			public void remove() {
				// do nothing since ResidueRange is not modifiable
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
				// do nothing since ResidueRange is not modifiable
			}
		};
	}

	/**
	 * Returns a new Iterator over every {@link ResidueNumber} in the list of ResidueRanges.
	 * Stores the contents of {@code map} until the iterator is finished, so calling code should set the iterator to {@code null} if it did not finish.
	 */
	public static Iterator<ResidueNumber> multiIterator(AtomPositionMap map, List<ResidueRange> rrs) {
		ResidueRange[] ranges = new ResidueRange[rrs.size()];
		for (int i = 0; i < rrs.size(); i++) {
			ranges[i] = rrs.get(i);
		}
		return multiIterator(map, ranges);
	}

	public static List<ResidueRange> parseMultiple(List<String> ranges) {
		List<ResidueRange> rrs = new ArrayList<ResidueRange>(ranges.size());
		for (String range : ranges) {
			rrs.add(ResidueRange.parse(range));
		}
		return rrs;
	}

	/**
	 * Determines whether a String is of a recognizable range format
	 */
	public static boolean looksLikeRange(String s) {
		Pattern pattern = Pattern.compile(RANGE_REGEX);
		Matcher matcher = pattern.matcher(s);
		matcher.find();
		return matcher.matches();
	}

}
