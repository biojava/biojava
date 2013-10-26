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

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.biojava.bio.structure.io.mmcif.chem.ResidueType;

/**
 * A map from {@link ResidueNumber ResidueNumbers} to ATOM record positions in a PDB file.
 * To use:
 * <code>
 * AtomPositionMap map = new AtomPositionMap(new AtomCache().getAtoms("1w0p"));
 * ResidueNumber start = new ResidueNumber("A", 100, null);
 * ResidueNumber end = map.getEnd("A");
 * int pos = map.getPosition(start);
 * int length = map.calcLength(start, end);
 * </code>
 * @author dmyerstu
 */
public class AtomPositionMap {

	private HashMap<ResidueNumber, Integer> hashMap;
	private TreeMap<ResidueNumber, Integer> treeMap;

	public static final Set<String> AMINO_ACID_NAMES = new TreeSet<String>();
	static {
		AMINO_ACID_NAMES.addAll(Arrays.asList(new String[] {"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}));
		AMINO_ACID_NAMES.addAll(Arrays.asList(new String[] {"ASX", "GLX", "XLE", "XAA"}));
	}

	public static interface GroupMatcher {
		boolean matches(Group group);
	}

	public static final GroupMatcher AMINO_ACID_MATCHER = new GroupMatcher() {
		@Override
		public boolean matches(Group group) {
			ResidueType type = group.getChemComp().getResidueType();
			return group.hasAtom(StructureTools.caAtomName) || AMINO_ACID_NAMES.contains(group.getPDBName()) || type == ResidueType.lPeptideLinking || type == ResidueType.glycine || type == ResidueType.lPeptideAminoTerminus || type == ResidueType.lPeptideCarboxyTerminus || type == ResidueType.dPeptideLinking || type == ResidueType.dPeptideAminoTerminus || type == ResidueType.dPeptideCarboxyTerminus;
		}
	};

	public static final GroupMatcher ANYTHING_MATCHER = new GroupMatcher() {
		@Override
		public boolean matches(Group group) {
			return true;
		}
	};
	
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
	 * Creates a new AtomPositionMap containing only amino acids C-alpha atoms. C-alpha atoms are identified somewhat liberally.
	 * @param atoms
	 */
	public AtomPositionMap(Atom[] atoms) {
		this(atoms, AMINO_ACID_MATCHER);
	}
	
	/**
	 * Creates a new AtomPositionMap containing only atoms matched by {@code matcher}.
	 * @param atoms
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
	 * This is <strong>not</em> the same as subtracting {@link #getPosition(ResidueNumber)} for {@code positionB} from {@link #getPosition(ResidueNumber)} for {@code positionA}.
	 * The latter considers only positions of ATOM entries in the PDB file and ignores chains. This method only includes ATOMs from the same chain.
	 * @param positionA
	 * @param positionB
	 * @param startingChain
	 * @return
	 */
	public int calcLength(int positionA, int positionB, char startingChain) {
		return calcLength(positionA, positionB, String.valueOf(startingChain));
	}
	
	/**
	 * This is <strong>not</em> the same as subtracting {@link #getPosition(ResidueNumber)} for {@code positionB} from {@link #getPosition(ResidueNumber)} for {@code positionA}.
	 * The latter considers only positions of ATOM entries in the PDB file and ignores chains. This method only includes ATOMs from the same chain.
	 * @param positionA
	 * @param positionB
	 * @param startingChain
	 * @return
	 */
	public int calcLength(int positionA, int positionB, String startingChain) {
		int positionStart, positionEnd;
		if (positionA <= positionB) {
			positionStart = positionA;
			positionEnd = positionB;
		} else {
			positionStart = positionB;
			positionEnd = positionA;
		}
		return calcLengthDirectional(positionStart, positionEnd, startingChain);
	}

	/**
	 * Calculates the distance between {@code positionStart} and {@code positionEnd}. Will return a negative value if the start is past the end.
	 * @param start
	 * @param end
	 * @return
	 */
	public int calcLengthDirectional(ResidueNumber start, ResidueNumber end) {
		return calcLengthDirectional(getPosition(start), getPosition(end), start.getChainId());
	}

	/**
	 * Calculates the distance between {@code positionStart} and {@code positionEnd}. Will return a negative value if the start is past the end.
	 * @param positionStart
	 * @param positionEnd
	 * @param startingChain
	 * @return
	 */
	public int calcLengthDirectional(int positionStart, int positionEnd, char startingChain) {
		return calcLengthDirectional(positionStart, positionEnd, String.valueOf(startingChain));
	}
	
	/**
	 * Calculates the distance between {@code positionStart} and {@code positionEnd}. Will return a negative value if the start is past the end.
	 * @param positionStart
	 * @param positionEnd
	 * @param startingChain
	 * @return
	 */
	public int calcLengthDirectional(int positionStart, int positionEnd, String startingChain) {
		int count = 0;
		for (Map.Entry<ResidueNumber, Integer> entry : treeMap.entrySet()) {
			if (entry.getKey().getChainId().equals(startingChain)) {
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
	 * Convenience method for {@link #calcLength(int, int, String)}.
	 * @param positionA
	 * @param positionB
	 * @return
	 * @see #calcLength(int, int, char)
	 */
	public int calcLength(ResidueNumber positionA, ResidueNumber positionB) {
		int pA = hashMap.get(positionA);
		int pB = hashMap.get(positionB);
		String chain = positionA.getChainId();
		if (pA > pB) chain = positionB.getChainId();
		return calcLength(pA, pB, chain);
	}

	public NavigableMap<ResidueNumber, Integer> getNavMap() {
		return treeMap;
	}

	/**
	 * @param residueNumber
	 * @return The position of the ATOM record in the PDB file corresponding to the {@code residueNumber}
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
	 * @param chainId
	 * @return The first {@link ResidueNumber} of any chain (the one farthest down in the PDB file)
	 */
	public ResidueNumber getFirst() {
		return treeMap.firstKey();
	}

	/**
	 * @param chainId
	 * @return The last {@link ResidueNumber} of any chain (the one farthest down in the PDB file)
	 */
	public ResidueNumber getLast() {
		return treeMap.lastKey();
	}

}
