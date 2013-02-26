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
import java.util.HashMap;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;

/**
 * A chain, a start residue, and an end residue. May also store a length value. Because of insertion codes, this length
 * is not necessarily end-start. Since residue numbers are unique within a PDB file, each ResidueNumber in a
 * {@link Structure} can be mapped to a unique ATOM record in that Structure's PDB file. ResidueRange provides a static
 * method {@link #getAminoAcidPositions(Atom[])} that returns a {@link HashMap} that maps each {@link ResidueNumber} to
 * its position in the PDB file's ATOM records. It also provides static methods for working with this map. Note that
 * these are defined as static utility methods only for performance reasons; it is not economical to store this data in
 * each ResidueRange.
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
	 * @return The unique ResidueRange corresponding to {@code s}.
	 */
	public static ResidueRange parse(String s) {
		char chain = s.charAt(0);
		String[] parts = s.substring(2).split("-");
		ResidueNumber start = ResidueNumber.fromString(parts[0]);
		start.setChainId(String.valueOf(chain));
		ResidueNumber end = ResidueNumber.fromString(parts[1]);
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

}
