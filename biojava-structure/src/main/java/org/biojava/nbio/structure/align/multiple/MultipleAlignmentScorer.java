package org.biojava.nbio.structure.align.multiple;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;

/**
 * Utility class for calculating common scores over multiple alignments.
 * 
 * @author Spencer Bliven
 *
 */
public class MultipleAlignmentScorer {
	// Names for commonly used properties
	//	public static final String SCORE_TMSCORE = "TMScore";
	//	public static final String SCORE_PROBABILITY = "Probability";
	//	public static final String SCORE_CE = "CEScore";
	//	public static final String SCORE_RMSD = "RMSD";

	public static final String CEMC_SCORE = "CEMCscore";
	public static final String SCORE_REF_RMSD = "RefRMSD";
	public static final String SCORE_REF_TMSCORE = "RefTMScore";

	public static void calculateScores(MultipleAlignment alignment) throws StructureAlignmentException, StructureException {
		List<Atom[]> transformed = transformAtoms(alignment);
		alignment.putScore(SCORE_REF_RMSD, getRefRMSD(transformed,0));
		List<Integer> lengths = new ArrayList<Integer>(alignment.size());
		for(Atom[] atoms : alignment.getEnsemble().getAtomArrays()) {
			lengths.add(atoms.length);
		}
		alignment.putScore(SCORE_REF_TMSCORE, getRefTMScore(transformed,lengths,0));
	}

	/**
	 * Transforms atoms according to the superposition stored in the alignment.
	 * <p>
	 * For each structure in the alignment, returns an atom for each
	 * representative atom in the aligned columns, omitting unaligned residues
	 * (i.e. an array of length <tt>alignment.getLength()</tt>).
	 * All blocks are concatenated together, so Atoms may not appear in the
	 * same order as in their parent structure. If the alignment blocks contain
	 * null residue (gaps), then the returned array will also contain gaps.
	 * @return
	 * @throws StructureAlignmentException 
	 */
	public static List<Atom[]> transformAtoms(MultipleAlignment alignment) throws StructureAlignmentException {
		if(alignment.getEnsemble() == null ) {
			throw new NullPointerException("No ensemble set for this alignment");
		}

		List<Atom[]> atomArrays = alignment.getEnsemble().getAtomArrays();
		List<Atom[]> transformed = new ArrayList<Atom[]>(atomArrays.size());

		//Loop through structures
		for (int i=0; i<atomArrays.size(); i++){

			Matrix4d transform = null;
			if( alignment.getTransformations() != null) {
				transform = alignment.getTransformations().get(i);
			}
			Atom[] curr = atomArrays.get(i); // all CA atoms from structure

			//Concatenated list of all blocks for this structure
			Atom[] transformedAtoms = new Atom[alignment.length()];
			int transformedAtomsLength = 0;

			// Each blockset gets transformed independently
			for( BlockSet bs : alignment.getBlockSets()) {

				Atom[] blocksetAtoms = new Atom[bs.length()];

				for( Block blk : bs.getBlocks() ) {
					if( blk.size() != atomArrays.size()) {
						throw new StructureAlignmentException(String.format(
								"Mismatched block length. Expected %d structures, found %d.",
								atomArrays.size(),blk.size() ));
					}
					//Extract aligned atoms
					for (int j=0; j<blk.length(); j++){
						Integer alignedPos = blk.getAlignRes().get(i).get(j);
						if (alignedPos != null) {
							blocksetAtoms[j] = (Atom) curr[alignedPos].clone();
						}
					}
				}
				
				// transform according to (1) the blockset matrix, or (2) the alignment matrix
				Matrix4d blockTrans = null;
				if(bs.getTransformations() != null)
					blockTrans = bs.getTransformations().get(i);
				if(blockTrans == null) {
					blockTrans = transform;
				}

				for(Atom a : blocksetAtoms) {
					if (a!=null) Calc.transform(a, blockTrans);
					transformedAtoms[transformedAtomsLength] = a;
					transformedAtomsLength++;
				}
			}
			assert(transformedAtomsLength == alignment.length());

			transformed.add(transformedAtoms);
		}
		return transformed;
	}
	

	public static double getRefRMSD(MultipleAlignment alignment, int reference) throws StructureAlignmentException {
		List<Atom[]> transformed = transformAtoms(alignment);
		return getRefRMSD(transformed,reference);
	}
	/**
	 * Calculates the average RMSD from all structures to a reference structure,
	 * given a set of superimposed atoms.
	 * <p>
	 * For ungapped alignments, this is just the sqroot of the average distance
	 * from an atom to the aligned atom from the reference. Thus, for R structures
	 * aligned over C columns (with structure 0 as the reference), we have <br/>
	 * <pre>RefRMSD = \sqrt{  1/(C*(R-1)) * \sum_{r=1}^{R-1} \sum_{j=0}^{C-1} (atom[1][c]-atom[r][c])^2  }</pre>
	 * 
	 * <p>
	 * For gapped alignments, null atoms are omitted from consideration, so that
	 * the RMSD is the average over all columns with non-null reference of the
	 * average RMSD within the non-null elements of the column.
	 * @param transformed
	 * @param reference
	 * @return
	 * @throws StructureAlignmentException
	 */
	public static double getRefRMSD(List<Atom[]> transformed, int reference) throws StructureAlignmentException {

		double sumSqDist = 0;
		int totalLength = 0;
		for(int c=0;c<transformed.get(reference).length;c++) {
			Atom refAtom = transformed.get(reference)[c];
			if(refAtom == null)
				continue;

			double nonNullSqDist = 0;
			int nonNullLength = 0;
			for(int r=0;r<transformed.size();r++) {
				if(r == reference)
					continue;
				Atom atom = transformed.get(r)[c];
				if(atom != null) {
					nonNullSqDist += Calc.getDistanceFast(refAtom, atom);
					nonNullLength++;
				}
			}

			if(nonNullLength > 0) {
				totalLength++;
				sumSqDist += nonNullSqDist/nonNullLength;
			}
		}
		return Math.sqrt(sumSqDist/totalLength);
	}


	public static double getRefTMScore(MultipleAlignment alignment, int reference)
			throws StructureAlignmentException, StructureException
	{
		List<Atom[]> transformed = transformAtoms(alignment);
		List<Integer> lengths = new ArrayList<Integer>(alignment.size());
		for(Atom[] atoms : alignment.getEnsemble().getAtomArrays()) {
			lengths.add(atoms.length);
		}
		return getRefTMScore(transformed,lengths,reference);
	}
	/**
	 * Calculates the average TMScore from all structures to a reference structure,
	 * given a set of superimposed atoms.
	 * @param transformed Arrays of aligned atoms, after superposition
	 * @param lengths lengths of the full input structures
	 * @param reference Index of the reference structure
	 * @return
	 * @throws StructureAlignmentException
	 * @throws StructureException 
	 */
	public static double getRefTMScore(List<Atom[]> transformed, List<Integer> lengths, int reference)
			throws StructureAlignmentException, StructureException {


		if(transformed.size() != lengths.size()) {
			throw new IllegalArgumentException("Input sizes differ");
		}
		double sumTM = 0;
		int len = transformed.get(reference).length;
		for(int r=0;r<transformed.size();r++) {
			if(r==reference)
				continue;
			//remove nulls from both arrays
			Atom[] ref = new Atom[len];
			Atom[] aln = new Atom[len];
			int nonNullLen = 0;
			for(int c=0;c<len;c++) {
				if( transformed.get(reference)[c] != null && transformed.get(r)[c] != null) {
					ref[nonNullLen] = transformed.get(reference)[c];
					aln[nonNullLen] = transformed.get(r)[c];
					nonNullLen++;
				}
			}
			//truncate nulls
			if(nonNullLen<len) {
				ref = Arrays.copyOf(ref, nonNullLen);
				aln = Arrays.copyOf(aln, nonNullLen);
			}
			sumTM += SVDSuperimposer.getTMScore(ref, aln, lengths.get(reference), lengths.get(r));
		}
		return sumTM/transformed.size();
	}

}
