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
 * Created on 2012-10-27
 * Created by dmyerstu
 *
 */
package org.biojava.bio.structure;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

/**
 * A map from each {@link ResidueNumber} corresponding to a {@link Group} in an array of {@link Atom Atoms} from a
 * single {@link Structure} to its position in the corresponding PDB file's ATOM records. Can map only Groups matching
 * particular criteria (see {@link AtomPositionMap.GroupMatcher GroupMatcher}). Like Structure, this is a
 * memory-intensive class. Thus client code should allow an AtomPositionMap to be released from scope in order to free
 * memory.
 * 
 * @author dmyerstu
 * @see ResidueRange
 * @since 3.0.6
 */
public class AtomPositionMap {

	/**
	 * Something that can match a subset (proper or not) of {@link Group Groups}.
	 * 
	 * @author dmyerstu
	 * @see {@link AtomPositionMap#AtomPositionMap(Atom[], GroupMatcher)}.
	 */
	public interface GroupMatcher {
		boolean matches(Group group);
	}

	/**
	 * A map that is sorted by its values.
	 * 
	 * @author dmyerstu
	 * 
	 * @param <T>
	 *            The key type
	 * @param <V>
	 *            The value type
	 */
	public class ValueComparator<T, V extends Comparable<V>> implements Comparator<T> {

		private Map<T, V> map;

		/**
		 * Creates a new ValueComparator with specified map. Note that the map is stored in the ValueComparator. Client
		 * code should thus allow the ValueComparator to exit from scope.
		 * 
		 * @param map
		 */
		public ValueComparator(Map<T, V> map) {
			this.map = map;
		}

		@Override
		public int compare(T o1, T o2) {
			return map.get(o1).compareTo(map.get(o2));
		}

	}

	/**
	 * Matches only amino acids.
	 */
	public static final GroupMatcher AMINO_ACID_MATCHER = new GroupMatcher() {
		@Override
		public boolean matches(Group group) {
			return group.getType().equals(GroupType.AMINOACID);
		}
	};

	private HashMap<ResidueNumber, Integer> hashMap;

	private TreeMap<ResidueNumber, Integer> treeMap;

	/**
	 * Factory method that constructs a new AtomPostionMap using the specified Atom array and
	 * {@link #AMINO_ACID_MATCHER}.
	 * 
	 * @param atoms
	 * @return
	 * @see #AtomPositionMap(Atom[], GroupMatcher).
	 */
	public static AtomPositionMap ofAminoAcids(Atom[] atoms) {
		return new AtomPositionMap(atoms, AMINO_ACID_MATCHER);
	}

	/**
	 * Constructs a new AtomPositionMap precisely between every Group in the Atom array for which
	 * {@link GroupMatcher#matches(Group)} is true.
	 * 
	 * @param atoms
	 * @param matcher
	 * @see GroupMatcher
	 */
	public AtomPositionMap(Atom[] atoms, GroupMatcher matcher) {
		hashMap = new HashMap<ResidueNumber, Integer>();
		for (int i = 0; i < atoms.length; i++) {
			Group group = atoms[i].getGroup();
			ResidueNumber rn = group.getResidueNumber();
			if (matcher.matches(group)) {
				if (!hashMap.containsKey(rn)) {
					hashMap.put(rn, i + 1);
				}
			}
		}
		Comparator<ResidueNumber> vc = new ValueComparator<ResidueNumber, Integer>(hashMap);
		treeMap = new TreeMap<ResidueNumber, Integer>(vc);
		treeMap.putAll(hashMap);
	}

	/**
	 * {@code positionA} and {@code positionB} must be from the same chain.
	 * 
	 * @param positionA
	 *            The {@code postionA}th ATOM in the PDB file.
	 * @param positionB
	 *            The {@code postionB}th ATOM in the PDB file
	 * @param chain
	 * @return The length between positionA and positionB.
	 */
	public int calcLength(int positionA, int positionB, String chain) {
		int positionStart, positionEnd;
		if (positionA <= positionB) {
			positionStart = positionA;
			positionEnd = positionB;
		} else {
			positionStart = positionB;
			positionEnd = positionA;
		}
		int count = 0;
		for (Map.Entry<ResidueNumber, Integer> entry : treeMap.entrySet()) {
			if (entry.getKey().getChainId().equals(chain)) {
				if (entry.getValue() == positionStart) {
					count = 0;
				}
				if (entry.getValue() == positionEnd) return count;
				count++;
			}
		}
		return -1;
	}

	/**
	 * {@code positionA} and {@code positionB} must be from the same chain.
	 * 
	 * @param positionA
	 * @param positionB
	 * @param startingChain
	 * @return The length between positionA and positionB.
	 */
	public int calcLength(ResidueNumber positionA, ResidueNumber positionB) {
		if (!positionA.getChainId().equals(positionB.getChainId())) throw new IllegalArgumentException(
				"The ResidueNumbers are not from the same chain.");
		int pA, pB;
		try {
			pA = hashMap.get(positionA);
		} catch (NullPointerException e) {
			throw new IllegalArgumentException("The specified atom " + positionA
					+ " does not exist in this AtomPositionMap.", e);
		}
		try {
			pB = hashMap.get(positionB);
		} catch (NullPointerException e) {
			throw new IllegalArgumentException("The specified atom " + positionB
					+ " does not exist in this AtomPositionMap.", e);
		}
		String chain = positionA.getChainId();
		if (pA > pB) chain = positionB.getChainId();
		return calcLength(pA, pB, chain);
	}

	/**
	 * @return A {@link NavigableMap} implementing this AtomPositionMap. Use this if the map needs to be accessed
	 *         directly for efficiency; otherwise {@link #getPosition(ResidueNumber)} is preferred.
	 */
	public NavigableMap<ResidueNumber, Integer> getNavMap() {
		return treeMap;
	}

	/**
	 * @param residueNumber
	 * @return The position in this AtomPositionMap at which {@code residueNumber} appears. Note that this depends on
	 *         how this AtomPositionMap was defined: the result may be different for different {@link GroupMatcher
	 *         GroupMatchers}.
	 */
	public int getPosition(ResidueNumber residueNumber) {
		return hashMap.get(residueNumber);
	}

}
