package org.biojava.nbio.structure.align.superimpose;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.model.Block;
import org.biojava.nbio.structure.align.model.BlockSet;
import org.biojava.nbio.structure.align.model.MultipleAlignment;
import org.biojava.nbio.structure.align.model.StructureAlignmentException;

/**
 * Superimposes each structure in a MultipleAlignment onto a reference structure.
 * 
 * Performs a global superposition over all block sets.
 * @author Spencer Bliven, Aleix Lafita
 *
 */
public class ReferenceSuperimposer implements MultipleSuperimposer {
	private int reference;
	public ReferenceSuperimposer() {
		this(0);
	}
	/**
	 * @param reference Index of the structure to use as a reference
	 */
	public ReferenceSuperimposer(int reference) {
		this.reference = reference;
	}
	
	/**
	 * Superimpose all structures in a BlockSet to the reference structure
	 * @ param bs
	 * @return
	 * @throws StructureAlignmentException 
	 */
	@Override
	public void superimpose(MultipleAlignment alignment) throws StructureException, StructureAlignmentException {
		if(alignment.getEnsemble() == null ) {
			throw new NullPointerException("No ensemble set for this alignment");
		}
		
		List<Atom[]> atomArrays = alignment.getEnsemble().getAtomArrays();
		List<Block> blocks = alignment.getBlocks();
		
		if(blocks.size() < 1)
			throw new StructureAlignmentException("No aligned residues");
		
		if( reference < 0 || atomArrays.size() <= reference ) {
			throw new IndexOutOfBoundsException(String.format(
					"Invalid reference structure: requested %d but only %d structures.",
					reference,atomArrays.size()) );
		}


		// Output variables
		List<Matrix4d> transforms = new ArrayList<Matrix4d>(atomArrays.size());
		double totalDistanceSq = 0;
		double totalLength = 0;
		double totalTM = 0;

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
			for( Block blk : blocks ) {
				if( blk.size() != atomArrays.size()) {
					throw new StructureAlignmentException(String.format(
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

			transforms.add(svd.getTransformation());
			double rmsd = SVDSuperimposer.getRMS(array1, array2);
			totalDistanceSq += rmsd*rmsd*array1.length;
			totalLength += array1.length;
			
			totalTM += SVDSuperimposer.getTMScore(array1, array2, ref.length,curr.length);
		}
		
		// Set output parameters
		alignment.setTransformations(transforms);
		for(BlockSet bs : alignment.getBlockSets()) {
			bs.setTransformations(null);
		}
		//TODO should these be standard scores? If not, calculate standard scores -SB 2015-05
		alignment.putScore("RefRMSD", totalDistanceSq/totalLength);
		alignment.putScore("RefTMScore", totalTM/atomArrays.size());
	}

}
