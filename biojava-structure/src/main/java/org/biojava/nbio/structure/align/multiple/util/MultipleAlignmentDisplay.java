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
package org.biojava.nbio.structure.align.multiple.util;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Matrix4d;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Utility functions to generalize the visualization of MultipleAlignments in
 * molecular viewers. The methods return different types of selectors for the
 * aligned residues in the alignment.
 *
 * @author Andreas Prlic
 * @author Aleix Lafita
 * @author Spencer Bliven
 * @since 4.2.0
 *
 */
public class MultipleAlignmentDisplay {

	private static final Logger logger = LoggerFactory
			.getLogger(MultipleAlignmentDisplay.class);

	/**
	 * New structures are downloaded if they were not cached in the alignment
	 * and they are entirely transformed here with the superposition information
	 * in the Multiple Alignment.
	 *
	 * @param multAln
	 * @return list of transformed AtomArrays
	 * @throws StructureException
	 */
	public static List<Atom[]> getRotatedAtoms(MultipleAlignment multAln)
			throws StructureException {

		int size = multAln.size();

		List<Atom[]> atomArrays = multAln.getAtomArrays();
		for (int i = 0; i < size; i++) {
			if (atomArrays.get(i).length < 1)
				throw new StructureException(
						"Length of atoms arrays is too short! Size: "
								+ atomArrays.get(i).length);
		}

		List<Atom[]> rotatedAtoms = new ArrayList<Atom[]>();

		// TODO implement independent BlockSet superposition of the structure
		List<Matrix4d> transf = multAln.getBlockSet(0).getTransformations();

		if (transf == null) {

			logger.error("Alignment Transformations are not calculated. "
					+ "Superimposing to first structure as reference.");

			multAln = multAln.clone();
			MultipleSuperimposer imposer = new ReferenceSuperimposer();
			imposer.superimpose(multAln);
			transf = multAln.getBlockSet(0).getTransformations();
			assert (transf != null);
		}

		// Rotate the atom coordinates of all the structures
		for (int i = 0; i < size; i++) {
			// TODO handle BlockSet-level transformations
			// make sure this method has the same behavior as the other display.
			// -SB 2015-06

			// Assume all atoms are from the same structure
			Structure displayS = atomArrays.get(i)[0].getGroup().getChain()
					.getStructure().clone();
			
			// Get all the atoms and include ligands and hetatoms
			Atom[] rotCA = StructureTools.getRepresentativeAtomArray(displayS);
			List<Group> hetatms = StructureTools.getUnalignedGroups(rotCA);
			int index = rotCA.length;
			rotCA = Arrays.copyOf(rotCA, rotCA.length + hetatms.size());
			for (Group g : hetatms) {
				rotCA[index] = g.getAtom(0);
				index++;
			}

			// Transform the structure to ensure a full rotation in the display
			Calc.transform(displayS, transf.get(i));
			rotatedAtoms.add(rotCA);
		}

		return rotatedAtoms;
	}
}
