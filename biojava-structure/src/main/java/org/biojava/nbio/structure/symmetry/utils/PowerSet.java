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
package org.biojava.nbio.structure.symmetry.utils;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;

/**
 * In mathematics, the power set (or powerset) of any set S, written P(S), is
 * the set of all subsets of S, including the empty set and S itself.
 * <p>
 * Code taken from StackOverflow best answer in:
 * http://stackoverflow.com/questions/4640034/calculating-all-of-the-subsets
 * -of-a-set-of-numbers. HashSet changed to LinkedHashSet for the consistent order
 * of the subsets and easier testing.
 * <p>
 * Currently used to calculate the possible LOCAL symmetries in
 * {@link QuatSymmetryDetector}.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class PowerSet<T> {

	public PowerSet() {
	}

	/**
	 * @return the set of power Sets of the original Set
	 */
	public Set<Set<T>> powerSet(Set<T> originalSet) {
		Set<Set<T>> sets = new LinkedHashSet<Set<T>>();
		if (originalSet.isEmpty()) {
			sets.add(new LinkedHashSet<T>());
			return sets;
		}
		List<T> list = new ArrayList<T>(originalSet);
		T head = list.get(0);
		Set<T> rest = new LinkedHashSet<T>(list.subList(1, list.size()));
		for (Set<T> set : powerSet(rest)) {
			Set<T> newSet = new LinkedHashSet<T>();
			newSet.add(head);
			newSet.addAll(set);
			sets.add(newSet);
			sets.add(set);
		}
		return sets;
	}
}
