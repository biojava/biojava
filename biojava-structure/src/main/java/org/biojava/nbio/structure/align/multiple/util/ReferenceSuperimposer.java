package org.biojava.nbio.structure.align.multiple.util;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;

/**
 * Superimposes each structure in a {@link MultipleAlignment} onto a reference 
 * structure. 
 * <p>
 * Performs a global superposition of the MultipleAlignment in case 
 * there is only one {@link BlockSet}, and a superposition for every BlockSet 
 * in case there is more than one (flexible alignment).
 * <p>
 * This class uses the {@link SVDSuperimposer} algorithm.
 * 
 * @author Spencer Bliven
 * @author Aleix Lafita
 * @since 4.1.0
 *
 */
public class ReferenceSuperimposer implements MultipleSuperimposer {

	private int reference;

	/**
	 * Default Constructor. 
	 * Uses the first structure as the reference.
	 */
	public ReferenceSuperimposer() {
		this(0);
	}

	/**
	 * Constructor using a specified structure as reference.
	 * 
	 * @param reference Index of the structure to use as a reference 
	 * 			(it has to be > 0)
	 */
	public ReferenceSuperimposer(int reference) {
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
					//Loop through all aligned residues
					for (int j=0; j<blk.length(); j++){
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
				SVDSuperimposer svd = new SVDSuperimposer(array1, array2);
				Calc.transform(array2, svd.getTransformation());
				transforms.add(svd.getTransformation());
			}
			//Set transformation of the BlockSet
			bs.setTransformations(transforms);
		}
	}
}
