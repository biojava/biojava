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

package org.biojava.bio.structure;

import org.biojava.bio.structure.io.mmcif.chem.PolymerType;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

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
 * @author dmyerstu
 */
public class AtomPositionMap {

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
			ResidueType type = group.getChemComp().getResidueType();
			return PolymerType.PROTEIN_ONLY.contains(type.getPolymerType())
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
	private static class ValueComparator<T, V extends Comparable<V>> implements Comparator<T> {

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
			if (entry.getKey().getChainId().equals(startingChain)
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
		if( ! start.getChainId().equals(end.getChainId())) {
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
		return getLength(startPos, endPos, start.getChainId());
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
		if( ! start.getChainId().equals(end.getChainId())) {
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
		return getLengthDirectional(startPos, endPos, start.getChainId());
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
			if (entry.getKey().getChainId().equals(chainId)) return entry.getKey();
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
			if (entry.getKey().getChainId().equals(chainId)) return entry.getKey();
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
			if (!rn.getChainId().equals(currentChain)) {
				if (first != null) {
					ResidueRangeAndLength newRange = new ResidueRangeAndLength(currentChain, first, prev, this.getLength(first, prev));
					ranges.add(newRange);
				}
				first = rn;
			}
			prev = rn;
			currentChain = rn.getChainId();
		}
		ResidueRangeAndLength newRange = new ResidueRangeAndLength(currentChain, first, prev, this.getLength(first, prev));
		ranges.add(newRange);
		return ranges;
	}

}
