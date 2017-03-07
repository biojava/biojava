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
 */

package org.biojava.nbio.structure.symmetry.core;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 *
 * @author Peter
 */
public class PermutationGroup implements Iterable<List<Integer>> {
	List<List<Integer>> permutations = new ArrayList<List<Integer>>();

	public void addPermutation(List<Integer> permutation) {
		if (!permutations.contains(permutation)) {
			permutations.add(permutation);
		}
	}

	public List<Integer> getPermutation(int index) {
		return permutations.get(index);
	}

	public int getOrder() {
		return permutations.size();
	}


	/**
	 * Starts with an incomplete set of group generators in `permutations` and
	 * expands it to include all possible combinations.
	 *
	 * Ways to complete group:
	 * - combinations of permutations pi x pj
	 * - combinations with itself p^k
	 *
	 */
	public void completeGroup() {
		// Copy initial set to allow permutations to grow
		List<List<Integer>> gens = new ArrayList<List<Integer>>(permutations);
		// Keep HashSet version of permutations for fast lookup.
		Set<List<Integer>> known = new HashSet<List<Integer>>(permutations);
		//breadth-first search through the map of all members
		List<List<Integer>> currentLevel = new ArrayList<List<Integer>>(permutations);
		while( currentLevel.size() > 0) {
			List<List<Integer>> nextLevel = new ArrayList<List<Integer>>();
			for( List<Integer> p : currentLevel) {
				for(List<Integer> gen : gens) {
					List<Integer> y = combine(p,gen);
					if(!known.contains(y)) {
						nextLevel.add(y);
						//bypass addPermutation(y) for performance
						permutations.add(y);
						known.add(y);
					}
				}
			}
			currentLevel = nextLevel;
		}
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Permutation Group: " + permutations.size() + " permutation");
		for (List<Integer> permutation : permutations) {
			sb.append(permutation.toString());
		}
		return sb.toString();
	}

	public static List<Integer> combine(List<Integer> permutation1, List<Integer> permutation2) {
		List<Integer> intermediate = new ArrayList<Integer>(permutation1.size());
		for (int i = 0, n = permutation1.size(); i < n; i++) {
			intermediate.add(permutation2.get(permutation1.get(i)));
		}
		return intermediate;
	}

	public static int getOrder(List<Integer> permutation) {
		List<Integer> copy = new ArrayList<Integer>(permutation);
		for (int i = 0, n = permutation.size(); i < n; i++) {
			copy = combine(copy, permutation);
			if (copy.equals(permutation)) {
				return i + 1;
			}
		}
		return 0;
	}

	public String getGroupTable() {
		StringBuilder builder = new StringBuilder();
		builder.append("  |");
		for (int i = 0; i < getOrder(); i++) {
			builder.append(" ");
			builder.append(i);
		}
		builder.append("\n");
		builder.append("---");
		for (int i = 0; i < getOrder(); i++) {
			builder.append("--");
		}
		builder.append("\n");
		for (int i = 0; i < getOrder(); i++) {
			builder.append(i);
			builder.append(" |");
			for (int j = 0; j < getOrder(); j++) {
				builder.append(" ");
				builder.append(permutations.indexOf(combine(permutations.get(i), permutations.get(j))));
			}
			builder.append("\n");
		}
		return builder.toString();
	}

	@Override
	public int hashCode() {
	    return getGroupTable().hashCode();
	}

	@Override
	public Iterator<List<Integer>> iterator() {
		return permutations.iterator();
	}

	@Override
	public boolean equals(Object obj) {
	    if (this == obj) {
		return true;
	    }
	    if (obj == null) {
		return false;
	    }
	    if (this.getClass() == obj.getClass()) {
		return permutations.equals(((PermutationGroup)obj).permutations);
	    }
	    return false;
	}

}
