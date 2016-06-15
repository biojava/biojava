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
package org.biojava.nbio.structure.symmetry.internal;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * A ResidueGroup is a set of residues that are part of a maximally connected
 * component of the self-Alignment Graph in symmetry analysis.
 * <p>
 * This class provides an interface for comparing and combining them to refine
 * self-Alignments into consistent MultipleAlignments of repeats.
 *
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class ResidueGroup {

	private final List<Integer> residues;

	/**
	 * Create a ResidueGroup object from a maximally connected component.
	 *
	 * @param component
	 *            set of residues connected
	 */
	public ResidueGroup(Set<Integer> component) {
		// Transform component into sorted List of residues
		residues = new ArrayList<Integer>(component);
		Collections.sort(residues);
	}

	/**
	 * The order of symmetry of the group is the number of connected residues.
	 *
	 * @return size of residues List
	 */
	public int order() {
		return residues.size();
	}

	/**
	 * Determine if two Residuegroups (maximally connected components of the
	 * alignment Graph) are compatible, based in the following criterion:
	 *
	 * <pre>
	 * Two maximally connected components of the self-alignment Graph are
	 * compatible if they can be combined in a consistent multiple alignment
	 * of repeats, i.e.there exists one residue in c1 between each sorted
	 * pair of residues in c2.
	 * </pre>
	 *
	 * Compatibility is an intransitive relation, which means that for three
	 * ResidueGroups {A,B,C}, if A is compatible with B and B is compatible with
	 * C, then A is not necessarily compatible with C.
	 *
	 * @param c2
	 *            second maximally connected component
	 * @return true if compatible, false otherwise
	 */
	public boolean isCompatible(ResidueGroup other) {

		// Same order needed is necessary
		if (this.order() != other.order())
			return false;

		// Use the method of the smallest ResidueGroup
		if (this.residues.get(0) > other.residues.get(0))
			return other.isCompatible(this);

		// Check for intercalation of residues
		for (int i = 0; i < order() - 1; i++) {
			if (other.residues.get(i) > residues.get(i + 1))
				return false;
			if (residues.get(i) > other.residues.get(i + 1))
				return false;
		}

		return true;
	}

	/**
	 * Combine the ResidueGroup with the alignment block.
	 *
	 * @param alignRes
	 *            the alignment block, will be modified.
	 */
	public void combineWith(List<List<Integer>> alignRes) {
		for (int i = 0; i < order(); i++)
			alignRes.get(i).add(residues.get(i));
	}

}
