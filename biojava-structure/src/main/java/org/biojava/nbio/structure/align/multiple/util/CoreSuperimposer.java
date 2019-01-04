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

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.geometry.SuperPositionSVD;
import org.biojava.nbio.structure.geometry.SuperPositions;

/**
 * Superimposes the core aligned residues of every structure in a
 * {@link MultipleAlignment} onto a reference structure. This method
 * can eliminate the pairwise similarities of some structures to the
 * reference, when doing the superposition, taking into account only
 * those shared parts between the structures.
 * <p>
 * Performs a global superposition of the MultipleAlignment in case
 * there is only one {@link BlockSet}, and a superposition for every
 * BlockSet in case there is more than one (flexible alignment).
 * <p>
 * This class uses the {@link SuperPositionSVD} algorithm.
 *
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class CoreSuperimposer implements MultipleSuperimposer {

	private int reference;

	/**
	 * Default Constructor.
	 * Uses the first structure as the reference.
	 */
	public CoreSuperimposer() {
		this(0);
	}

	/**
	 * Constructor using a specified structure as reference.
	 *
	 * @param reference Index of the structure to use as a reference
	 * 			(it has to be > 0)
	 */
	public CoreSuperimposer(int reference) {
		if (reference<0) {
			throw new IllegalArgumentException(
					"reference index has to be positive, but was "+reference);
		}
		this.reference = reference;
	}

	@Override
	public void superimpose(MultipleAlignment alignment)
			throws StructureException {

		//Check for inconsistencies in the alignment
		if(alignment.getEnsemble() == null) {
			throw new NullPointerException("No ensemble set for this alignment."
					+ " Structure information cannot be obtained.");
		}
		if (alignment.size() < 1) {
			throw new IndexOutOfBoundsException(
					"No aligned structures, alignment size == 0.");
		}
		if (alignment.getCoreLength() < 1){
			throw new IndexOutOfBoundsException(
					"Alignment too short, core alignment length < 1.");
		}

		List<Atom[]> atomArrays = alignment.getAtomArrays();
		if (atomArrays.size() <= reference) {
			throw new IndexOutOfBoundsException(String.format(
					"Invalid reference structure: requested %d but "
							+ "only %d structures.",
							reference,atomArrays.size()));
		}

		alignment.clear();

		//Calculate BlockSet transformations
		for (BlockSet bs:alignment.getBlockSets()){

			//Block transformations
			List<Matrix4d> transforms =
					new ArrayList<Matrix4d>(atomArrays.size());

			//Loop through structures
			for (int i=0; i<atomArrays.size(); i++){

				if( i == reference) {
					//Identity operation
					Matrix4d ident = new Matrix4d();
					ident.setIdentity();
					transforms.add(ident);
					continue;
				}

				Atom[] ref = atomArrays.get(reference);
				Atom[] curr = atomArrays.get(i);

				List<Atom> atomSet1 = new ArrayList<Atom>();
				List<Atom> atomSet2 = new ArrayList<Atom>();

				for( Block blk : bs.getBlocks() ) {

					if( blk.size() != atomArrays.size()) {
						throw new IllegalStateException(String.format(
								"Mismatched block length. Expected %d "
										+ "structures, found %d.",
										atomArrays.size(),blk.size() ));
					}

					List<Integer> corePositions =
							MultipleAlignmentTools.getCorePositions(blk);

					//Loop through all aligned residues
					for (int j=0; j<blk.length(); j++){
						//Check that the position is in the core
						if (!corePositions.contains(j)) continue;

						Integer pos1 = blk.getAlignRes().get(reference).get(j);
						Integer pos2 = blk.getAlignRes().get(i).get(j);

						if (pos1==null || pos2==null) continue;

						atomSet1.add(ref[pos1]);
						atomSet2.add(curr[pos2]);
					}
				}
				Atom[] array1 = atomSet1.toArray(new Atom[atomSet1.size()]);
				Atom[] array2 = atomSet2.toArray(new Atom[atomSet2.size()]);

				array2 = StructureTools.cloneAtomArray(array2);

				//From the superimposer we obtain the rotation and translation
				Matrix4d trans = SuperPositions.superpose(Calc.atomsToPoints(array1),
						Calc.atomsToPoints(array2));
				transforms.add(trans);
			}
			//Set transformation of the BlockSet
			bs.setTransformations(transforms);
		}
	}
}
