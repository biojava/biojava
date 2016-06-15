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
 * Created on 2012-12-01
 *
 */

package org.biojava.nbio.structure;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.io.mmcif.chem.ResidueType;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A map from {@link ResidueNumber ResidueNumbers} to ATOM record positions in a PDB file.
 *
 * <p>To use:
 *
 * <pre>
 * Atom[] atoms = new AtomCache().getAtoms("1w0p");
 * AtomPositionMap map = new AtomPositionMap(atoms);
 * ResidueNumber start = new ResidueNumber("A", 100, null);
 * ResidueNumber end = map.getEnd("A");
 * int pos = map.getPosition(start);
 * int length = map.calcLength(start, end);
 * </pre>
 *
 * <p>Note: The getLength() methods were introduced in BioJava 4.0.0 to replace
 * the calcLength methods. The new method returns the number of residues between
 * two residues, inclusive, whereas the previous method returned 1 less than that.
 * @author Douglas Myers-Turnbull
 */
public class AtomPositionMap {

	private static final Logger logger = LoggerFactory.getLogger(AtomPositionMap.class);

	private HashMap<ResidueNumber, Integer> hashMap;
	private TreeMap<ResidueNumber, Integer> treeMap;


	/**
	 * Used as a Predicate to indicate whether a particular Atom should be mapped
	 */
	public static interface GroupMatcher {
		boolean matches(Group group);
	}

	/**
	 * Matches CA atoms of protein groups
	 */
	public static final GroupMatcher AMINO_ACID_MATCHER = new GroupMatcher() {
		@Override
		public boolean matches(Group group) {
			if( group == null )
				return false;
			ChemComp chem = group.getChemComp();
			if(chem == null)
				return false;
			// Get polymer type
			PolymerType polyType = chem.getPolymerType();
			if( polyType == null) {
				ResidueType type = chem.getResidueType();
				if(type != null ) {
					polyType = type.getPolymerType();
				}
			}
			if( polyType == null ) {
				return false;
			}

			return PolymerType.PROTEIN_ONLY.contains(polyType)
					&& group.hasAtom(StructureTools.CA_ATOM_NAME);
		}
	};

	/**
	 * Matches all atoms
	 */
	public static final GroupMatcher ANYTHING_MATCHER = new GroupMatcher() {
		@Override
		public boolean matches(Group group) {
			return true;
		}
	};

	/**
	 * A map that is sorted by its values. Used for the treemap
	 *
	 * @author dmyerstu
	 *
	 * @param <T>
	 *            The key type
	 * @param <V>
	 *            The value type
	 */
	private static class ValueComparator<T, V extends Comparable<V>> implements Comparator<T>, Serializable {
        private static final long serialVersionUID = 1;

		private Map<T, V> map;

		public ValueComparator(Map<T, V> map) {
			this.map = map;
		}

		@Override
		public int compare(T o1, T o2) {
			return map.get(o1).compareTo(map.get(o2));
		}

	}

	/**
	 * Creates a new AtomPositionMap containing peptide alpha-carbon atoms
	 *
	 * @param atoms
	 */
	public AtomPositionMap(Atom[] atoms) {
		this(atoms, AMINO_ACID_MATCHER);
	}

	/**
	 * Creates a new AtomPositionMap containing only atoms matched by {@code matcher}.
	 *
	 * If multiple atoms are present from a group, the first atom encountered will
	 * be used.
	 * @param atoms
	 */
	public AtomPositionMap(Atom[] atoms, GroupMatcher matcher) {
		hashMap = new HashMap<ResidueNumber, Integer>();
		for (int i = 0; i < atoms.length; i++) {
			Group group = atoms[i].getGroup();
			ResidueNumber rn = group.getResidueNumber();
			if (matcher.matches(group)) {
				if (!hashMap.containsKey(rn)) {
					hashMap.put(rn, i);
				}
			}
		}
		Comparator<ResidueNumber> vc = new ValueComparator<ResidueNumber, Integer>(hashMap);
		treeMap = new TreeMap<ResidueNumber, Integer>(vc);
		treeMap.putAll(hashMap);
	}

	/**
	 * Creates a new AtomPositionMap containing representative atoms
	 * from a structure.
	 * @param s
	 */
	public AtomPositionMap(Structure s) {
		this(StructureTools.getRepresentativeAtomArray(s));
	}

	/**
	 * Calculates the number of residues of the specified chain in a given range, inclusive.
	 * @param positionA index of the first atom to count
	 * @param positionB index of the last atom to count
	 * @param startingChain Case-sensitive chain
	 * @return The number of atoms between A and B inclusive belonging to the given chain
	 */
	public int getLength(int positionA, int positionB, String startingChain) {

		int positionStart, positionEnd;
		if (positionA <= positionB) {
			positionStart = positionA;
			positionEnd = positionB;
		} else {
			positionStart = positionB;
			positionEnd = positionA;
		}

		int count = 0;
		// Inefficient search
		for (Map.Entry<ResidueNumber, Integer> entry : treeMap.entrySet()) {
			if (entry.getKey().getChainName().equals(startingChain)
					&& positionStart <= entry.getValue()
					&& entry.getValue() <= positionEnd)
			{
				count++;
			}
		}
		return count;
	}


	/**
	 * Calculates the number of residues of the specified chain in a given range.
	 * Will return a negative value if the start is past the end.
	 * @param positionStart index of the first atom to count
	 * @param positionEnd index of the last atom to count
	 * @param startingChain Case-sensitive chain
	 * @return The number of atoms from A to B inclusive belonging to the given chain
	 */
	public int getLengthDirectional(int positionStart, int positionEnd, String startingChain) {
		int count = getLength(positionStart,positionEnd,startingChain);
		if(positionStart <= positionEnd) {
			return count;
		} else {
			return -count;
		}
	}

	/**
	 * Calculates the number of atoms between two ResidueNumbers, inclusive. Both residues
	 * must belong to the same chain.
	 * @param start First residue
	 * @param end Last residue
	 * @return The number of atoms from A to B inclusive
	 * @throws IllegalArgumentException if start and end are on different chains,
	 *  or if either of the residues doesn't exist
	 */
	public int getLength(ResidueNumber start, ResidueNumber end) {
		if( ! start.getChainName().equals(end.getChainName())) {
			throw new IllegalArgumentException(String.format(
					"Chains differ between %s and %s. Unable to calculate length.",
					start,end));
		}
		Integer startPos = getPosition(start);
		Integer endPos = getPosition(end);
		if(startPos == null) {
			throw new IllegalArgumentException("Residue "+start+" was not found.");
		}
		if(endPos == null) {
			throw new IllegalArgumentException("Residue "+start+" was not found.");
		}
		return getLength(startPos, endPos, start.getChainName());
	}

	/**
	 * Calculates the number of atoms between two ResidueNumbers, inclusive. Both residues
	 * must belong to the same chain.
	 * Will return a negative value if the start is past the end.
	 * @param start First residue
	 * @param end Last residue
	 * @return The number of atoms from A to B inclusive
	 * @throws IllegalArgumentException if start and end are on different chains,
	 *  or if either of the residues doesn't exist
	 */
	public int getLengthDirectional(ResidueNumber start, ResidueNumber end) {
		if( ! start.getChainName().equals(end.getChainName())) {
			throw new IllegalArgumentException(String.format(
					"Chains differ between %s and %s. Unable to calculate length.",
					start,end));
		}
		Integer startPos = getPosition(start);
		Integer endPos = getPosition(end);
		if(startPos == null) {
			throw new IllegalArgumentException("Residue "+start+" was not found.");
		}
		if(endPos == null) {
			throw new IllegalArgumentException("Residue "+start+" was not found.");
		}
		return getLengthDirectional(startPos, endPos, start.getChainName());
	}


	public NavigableMap<ResidueNumber, Integer> getNavMap() {
		return treeMap;
	}

	/**
	 * Gets the 0-based index of residueNumber to the matched atoms
	 * @param residueNumber
	 * @return The position of the ATOM record in the PDB file corresponding to
	 * the {@code residueNumber}, or null if the residueNumber was not found
	 */
	public Integer getPosition(ResidueNumber residueNumber) {
		return hashMap.get(residueNumber);
	}

	/**
	 * @param chainId
	 * @return The first {@link ResidueNumber} of the specified chain (the one highest down in the PDB file)
	 */
	public ResidueNumber getFirst(String chainId) {
		Map.Entry<ResidueNumber,Integer> entry = treeMap.firstEntry();
		while (true) {
			if (entry.getKey().getChainName().equals(chainId)) return entry.getKey();
			entry = treeMap.higherEntry(entry.getKey());
			if (entry == null) return null;
		}
	}

	/**
	 * @param chainId
	 * @return The last {@link ResidueNumber} of the specified chain (the one farthest down in the PDB file)
	 */
	public ResidueNumber getLast(String chainId) {
		Map.Entry<ResidueNumber,Integer> entry = treeMap.lastEntry();
		while (true) {
			if (entry.getKey().getChainName().equals(chainId)) return entry.getKey();
			entry = treeMap.lowerEntry(entry.getKey());
			if (entry == null) return null;
		}
	}

	/**
	 * @return The first {@link ResidueNumber} of any chain (the one farthest down in the PDB file)
	 */
	public ResidueNumber getFirst() {
		return treeMap.firstKey();
	}

	/**
	 * @return The last {@link ResidueNumber} of any chain (the one farthest down in the PDB file)
	 */
	public ResidueNumber getLast() {
		return treeMap.lastKey();
	}

	/**
	 * Returns a list of {@link ResidueRange ResidueRanges} corresponding to this entire AtomPositionMap.
	 */
	public List<ResidueRangeAndLength> getRanges() {
		String currentChain = "";
		ResidueNumber first = null;
		ResidueNumber prev = null;
		List<ResidueRangeAndLength> ranges = new ArrayList<ResidueRangeAndLength>();
		for (ResidueNumber rn : treeMap.keySet()) {
			if (!rn.getChainName().equals(currentChain)) {
				if (first != null) {
					ResidueRangeAndLength newRange = new ResidueRangeAndLength(currentChain, first, prev, this.getLength(first, prev));
					ranges.add(newRange);
				}
				first = rn;
			}
			prev = rn;
			currentChain = rn.getChainName();
		}
		ResidueRangeAndLength newRange = new ResidueRangeAndLength(currentChain, first, prev, this.getLength(first, prev));
		ranges.add(newRange);
		return ranges;
	}

	/**
	 * Trims a residue range so that both endpoints are contained in this map.
	 * @param rr residue range
	 * @return residue range and length
	 */
	public ResidueRangeAndLength trimToValidResidues(ResidueRange rr) {
		ResidueNumber start = rr.getStart();
		ResidueNumber end = rr.getEnd();
		String chain = rr.getChainName();
		// Add chainName
		if(start.getChainName() == null) {
			start = new ResidueNumber(chain,start.getSeqNum(),start.getInsCode());
		}
		if(end.getChainName() == null) {
			end = new ResidueNumber(chain,end.getSeqNum(),end.getInsCode());
		}
		// Check that start and end are present in the map.
		// If not, try to find the next valid residue
		// (terminal residues sometimes lack CA atoms, so they don't appear)
		Integer startIndex = getPosition(start);
		if( startIndex == null) {
			// Assume that the residue numbers are sequential
			// Find startIndex such that the SeqNum is bigger than start's seqNum
			for(ResidueNumber key :  treeMap.keySet()) {
				if( !key.getChainName().equals(chain) )
					continue;
				if( start.getSeqNum() <= key.getSeqNum() ) {
					start = key;
					startIndex = getPosition(key);
					break;
				}
			}
			if( startIndex == null ) {
				logger.error("Unable to find Residue {} in AtomPositionMap, and no plausible substitute.",start);
				return null;
			} else {
				logger.warn("Unable to find Residue {}, so substituting {}.",rr.getStart(),start);
			}
		}
		Integer endIndex = getPosition(end);
		if( endIndex == null) {
			// Assume that the residue numbers are sequential
			// Find startIndex such that the SeqNum is bigger than start's seqNum
			for(ResidueNumber key :  treeMap.descendingKeySet()) {
				if( !key.getChainName().equals(chain) )
					continue;
				Integer value = getPosition(key);
				if( value < startIndex ) {
					// start is before the end!
					break;
				}
				if( end.getSeqNum() >= key.getSeqNum() ) {
					end = key;
					endIndex = value;
					break;
				}
			}
			if( endIndex == null ) {
				logger.error("Unable to find Residue {} in AtomPositionMap, and no plausible substitute.",end);
				return null;
			} else {
				logger.warn("Unable to find Residue {}, so substituting {}.",rr.getEnd(),end);
			}
		}

		// now use those to calculate the length
		// if start or end is null, will throw NPE
		int length = getLength(startIndex, endIndex,chain);

		return new ResidueRangeAndLength(chain, start, end, length);
	}
}
