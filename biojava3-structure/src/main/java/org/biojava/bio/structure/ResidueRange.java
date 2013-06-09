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
	private final char chain;
	private final ResidueNumber end;
	private final Integer length;
	private final ResidueNumber start;

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
	 *            A string of the form chain_start-end. For example: <code>A.5-100</code>.
	 * @return The unique ResidueRange corresponding to {@code s}, or null if it is a whole chain
	 */
	public static ResidueRange parse(String s) {
		char chain = s.charAt(0);
		ResidueNumber start, end;
		if (s.length() > 2) {
			Pattern pattern = Pattern.compile("^([-]?[\\d]+[\\w]?)-([-]?[\\d]+[\\w]?)$");
			Matcher matcher = pattern.matcher(s.substring(2));
			matcher.find();
			try {
				start = ResidueNumber.fromString(matcher.group(1));
				end = ResidueNumber.fromString(matcher.group(2));
			} catch (IllegalStateException e) {
				throw new IllegalArgumentException("Range " + s + " was not valid", e);
			}
		} else {
			return null; // this is ok
		}
		start.setChainId(String.valueOf(chain));
		end.setChainId(String.valueOf(chain));
		return new ResidueRange(chain, start, end, null);
	}

	/**
	 * @param s
	 *            A string of the form chain_start-end. For example: <code>A.5-100</code>.
	 * @return The unique ResidueRange corresponding to {@code s}.
	 */
	public static ResidueRange parse(String s, AtomPositionMap map) {
		ResidueRange rr = parse(s);
		if (rr == null) { // whole chain
			String chain = s.substring(0, 1);
			rr = new ResidueRange(chain.charAt(0), map.getFirst(chain), map.getLast(chain), null);
		}
		int length = map.calcLength(rr.getStart(), rr.getEnd());
		return new ResidueRange(rr.getChain(), rr.getStart(), rr.getEnd(), length);
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
		if (chain != other.chain) return false;
		if (end == null) {
			if (other.end != null) return false;
		} else if (!end.equals(other.end)) return false;
		if (start == null) {
			if (other.start != null) return false;
		} else if (!start.equals(other.start)) return false;
		return true;
	}

	public char getChain() {
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
		result = prime * result + chain;
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

}
