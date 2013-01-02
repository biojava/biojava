/**
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
 * Created on 2012-10-18
 * Created by Douglas Myers-Turnbull
 *
 */
package org.biojava.bio.structure;

import java.util.ArrayList;
import java.util.List;

/**
 * A chain, a start {@link ResidueNumber} and an end ResidueNumber, optionally with length (number of residues). Due to
 * insertion codes, the length of a structure is not necessarily its end position minus its start. However, since
 * residue numbers are unique within any particular PDB file, each ResidueNumber in a {@link Structure} can be mapped to
 * a unique ATOM record in that Structure's PDB file; this is done by {@link AtomPositionMap}. ResidueRange can use this
 * information to determine the length of a range of residues accurately. Example use: <code>
 * AtomPositionMap map = AtomPositionMap.ofAminoAcids(cache.getAtoms("1qdm"));
 * ResidueRange range = ResidueRange.parse("A_246-85S", map);
 * System.out.println(range.getLength()); // should print 59
 * </code>
 * 
 * @author dmyerstu
 * @see ResidueNumber
 * @since 3.0.6
 */
public class ResidueRange {

	private final String chain;
	private final ResidueNumber end;
	private final Integer length;
	private final ResidueNumber start;

	/**
	 * Calculates the combined number of residues of the ResidueRanges in {@code ranges},
	 * <em>given that each ResidueRange has a length calculated</em>. Does not check for overlap between the
	 * ResidueRanges.
	 * 
	 * @param ranges
	 *            A list of ResidueRanges
	 * @return The combined length
	 * @throws IllegalArgumentException
	 *             If the {@link #getLength() length} of one or more ResidueRange is null
	 * @see #getLength()
	 */
	public static int calcLength(List<ResidueRange> ranges) {
		int length = 0;
		for (ResidueRange range : ranges) {
			if (range.getLength() == null) throw new IllegalArgumentException(
					"At least one ResidueRange does not have a length.");
			length += range.getLength();
		}
		return length;
	}

	/**
	 * @param string
	 *            A string of the form chain_start-end or chain.start-end. For example: <code>A.5-100</code> and
	 *            <code>A_5-100</code>.
	 * @return The unique ResidueRange corresponding to {@code string}. The result will have a null {@link #getLength()
	 *         length}.
	 */
	public static ResidueRange parse(String string) {
		String[] supParts = string.split("[\\._]");
		String chain = supParts[0];
		String[] parts = supParts[1].split("-");
		ResidueNumber start = ResidueNumber.fromString(parts[0]);
		start.setChainId(String.valueOf(chain));
		ResidueNumber end = ResidueNumber.fromString(parts[1]);
		end.setChainId(String.valueOf(chain));
		return new ResidueRange(chain, start, end, null);
	}

	/**
	 * @param string
	 *            A string of the form chain_start-end. For example: <code>A.5-100</code>.
	 * @return The unique ResidueRange corresponding to {@code string}. The result will have the correct
	 *         {@link #getLength() length}.
	 */
	public static ResidueRange parse(String string, AtomPositionMap map) {
		ResidueRange rr = parse(string);
		int length = map.calcLength(rr.getStart(), rr.getEnd());
		return new ResidueRange(rr.getChain(), rr.getStart(), rr.getEnd(), length);
	}

	/**
	 * @param string
	 *            A string of the form chain_start-end,chain_start-end, ... For example:
	 *            <code>A.5-100,R_110-190,Z_200-250</code>.
	 * @return The unique list of ResidueRanges corresponding to {@code string}. Each ResidueRange will have a null
	 *         {@link #getLength() length}.
	 */
	public static List<ResidueRange> parseMultiple(String string) {
		String[] parts = string.split(",");
		List<ResidueRange> list = new ArrayList<ResidueRange>(parts.length);
		for (String part : parts) {
			list.add(parse(part));
		}
		return list;
	}

	/**
	 * @param string
	 *            A string of the form chain_start-end,chain_start-end, ... For example:
	 *            <code>A.5-100,R_110-190,Z_200-250</code>.
	 * @param map
	 *            An {@link AtomPositionMap} containing one e
	 * @return The unique ResidueRange corresponding to {@code string}. Each ResidueRange will have the correct
	 *         {@link #getLength() length}.
	 */
	public static List<ResidueRange> parseMultiple(String string, AtomPositionMap map) {
		String[] parts = string.split(",");
		List<ResidueRange> list = new ArrayList<ResidueRange>(parts.length);
		for (String part : parts) {
			list.add(parse(part, map));
		}
		return list;
	}

	public ResidueRange(String chain, ResidueNumber start, ResidueNumber end, Integer length) {
		this.chain = chain;
		this.start = start;
		this.end = end;
		this.length = length;
	}

	/**
	 * Two ResidueRanges are equal iff they have the same chain, start residue, and end residue. Note that length is not
	 * included.
	 */
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

	public String getChain() {
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

}
