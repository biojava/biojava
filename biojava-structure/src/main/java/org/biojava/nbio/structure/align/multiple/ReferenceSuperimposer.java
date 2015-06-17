package org.biojava.nbio.structure.align.multiple;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;

/**
 * Superimposes each structure in a {@link MultipleAlignment} onto a reference structure.
 * Performs a global superposition of the MultipleAlignment in case there is only one 
 * {@link BlockSet}, and a superposition for every BlockSet in case there are more than 
 * one (flexible alignment).
 * <p>
 * This class uses the {@link SVDSuperimposer} algorithm.
 * 
 * @author Spencer Bliven
 * @author Aleix Lafita
 *
 */
public class ReferenceSuperimposer implements MultipleSuperimposer {
	
	private int reference;
	
	/**
	 * Default Constructor. Uses the first structure as the reference.
	 */
	public ReferenceSuperimposer() {
		this(0);
	}
	
	/**
	 * Constructor using a specified structure as reference.
	 * @param reference Index of the structure to use as a reference (it has to be >0)
	 */
	public ReferenceSuperimposer(int reference) {
		if (reference<0) throw new IllegalArgumentException("reference index has to be positive, but was "+reference);
		this.reference = reference;
	}
	
	/**
	 * Superimpose all structures in a MultipleAlignment to the reference structure. The superposition is done
	 * for all individual BlockSets. If there is only one BlockSet, the superposition is alse set in the 
	 * MultipleAlignment as the global superposition.
	 * <p>
	 * This method only calculates and sets the transformation 4D Matrices. If any score is needed
	 * it should be calculated and set separately afterwards with {@link MultipleAlignmentScorer}.
	 * 
	 * @param alignment MultipleAlignment object to superimpose.
	 */
	@Override
	public void superimpose(MultipleAlignment alignment) throws StructureException {
		
		//Check for inconsistencies in the alignment
		if(alignment.getEnsemble() == null) throw new NullPointerException("No ensemble set for this alignment");
		List<Atom[]> atomArrays = alignment.getEnsemble().getAtomArrays();
		
		if (alignment.getBlockSets().size() < 1)
			throw new IndexOutOfBoundsException("No aligned residues, number of Blocks == 0.");
		
		if (atomArrays.size() <= reference) {
			throw new IndexOutOfBoundsException(String.format(
					"Invalid reference structure: requested %d but only %d structures.",
					reference,atomArrays.size()) );
		}
		
		//Clear all the Cached information of the MultipleAlignment that depends on the superposition
		alignment.clear();

		//Loop through all BlockSets and calculate the transformations
		for (BlockSet bs:alignment.getBlockSets()){
			
			//Block transformations
			List<Matrix4d> transforms = new ArrayList<Matrix4d>(atomArrays.size());
	
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
	
				//Loop through all the Blocks that define the aligned positions
				for( Block blk : bs.getBlocks() ) {
					if( blk.size() != atomArrays.size()) {
						throw new IllegalStateException(String.format(
								"Mismatched block length. Expected %d structures, found %d.",
								atomArrays.size(),blk.size() ));
					}
					//Loop through all the columns of each Block (aligned residues)
					for (int j=0; j<blk.length(); j++){
						Integer pos1 = blk.getAlignRes().get(0).get(j);
						Integer pos2 = blk.getAlignRes().get(i).get(j);
	
						if (pos1==null || pos2==null) continue;
						atomSet1.add(ref[pos1]);
						atomSet2.add(curr[pos2]);
					}
				}
				Atom[] array1 = atomSet1.toArray(new Atom[atomSet1.size()]);
				Atom[] array2 = atomSet2.toArray(new Atom[atomSet2.size()]);
	
				array2 = StructureTools.cloneAtomArray(array2);
	
				//From the superimposer we obtain the rotation and translation information
				SVDSuperimposer svd = new SVDSuperimposer(array1, array2);
				Calc.transform(array2, svd.getTransformation());
				transforms.add(svd.getTransformation());
			}
			//Set transformation of the BlockSet
			bs.setTransformations(transforms);
		}
		//Set transformation of the MultipleAlignment
		if (alignment.getBlockSets().size() == 1) 
			alignment.setTransformations(alignment.getBlockSets().get(0).getTransformations());
	}
}
