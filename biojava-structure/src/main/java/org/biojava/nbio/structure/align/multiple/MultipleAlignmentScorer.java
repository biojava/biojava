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
 * Utility class for calculating common scores of {@link MultipleAlignment}s.
 * 
 * @author Spencer Bliven
 * @author Aleix Lafita
 *
 */
public class MultipleAlignmentScorer {
	
	//Names for commonly used properties
	public static final String PROBABILITY = "Probability";  		//AFPChain conversion
	public static final String CE_SCORE = "CEscore";			    //AFPChain conversion
	public static final String RMSD = "RMSD";
	public static final String AVG_TMSCORE = "Avg-TMscore";
	public static final String CEMC_SCORE = "CEMCscore";
	public static final String REF_RMSD = "Ref-RMSD";
	public static final String REF_TMSCORE = "Ref-TMscore";

	/**
	 * Calculates and puts the RMSD and the average TM-Score of the MultipleAlignment.
	 * 
	 * @param alignment
	 * @throws StructureException
	 * @see #getAvgTMScore(MultipleAlignment)
	 * @see #getRMSD(MultipleAlignment)
	 */
	public static void calculateScores(MultipleAlignment alignment) throws StructureException {
		
		//Put RMSD
		List<Atom[]> transformed = transformAtoms(alignment);
		alignment.putScore(RMSD, getRMSD(transformed));
		
		//Put TM-Score
		List<Integer> lengths = new ArrayList<Integer>(alignment.size());
		for(Atom[] atoms : alignment.getEnsemble().getAtomArrays()) {
			lengths.add(atoms.length);
		}
		alignment.putScore(AVG_TMSCORE, getAvgTMScore(transformed,lengths));
	}
	
	/**
	 * Transforms atoms according to the superposition stored in the alignment.
	 * <p>
	 * For each structure in the alignment, returns an atom for each
	 * representative atom in the aligned columns, omitting unaligned residues
	 * (i.e. an array of length <tt>alignment.length()</tt> ).
	 * <p>
	 * All blocks are concatenated together, so Atoms may not appear in the
	 * same order as in their parent structure. If the alignment blocks contain
	 * null residues (gaps), then the returned array will also contain null Atoms.
	 * 
	 * @param alignment MultipleAlignment
	 * @return
	 */
	public static List<Atom[]> transformAtoms(MultipleAlignment alignment) {
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
						throw new IllegalStateException(String.format(
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
	
	/**
	 * Calculates the RMSD of all-to-all structure comparisons (distances) of the
	 * given MultipleAlignment. <p>
	 * Complexity: T(n) = O(n^2), if n=number of structures.
	 * <p>
	 * The formula used is just the sqroot of the average distance
	 * of all possible pairs of atoms. Thus, for R structures
	 * aligned over C columns without gaps, we have
	 * <pre>RMSD = \sqrt{  1/(C*(R*(R-1)/2)) * \sum_{r1=1}^{R-1} \sum_{r2=r1+1}^{R-1} \sum_{j=0}^{C-1} (atom[r1][c]-atom[r2][c])^2  }</pre>
	 * 
	 * @param alignment
	 * @return double RMSD
	 */
	public static double getRMSD(MultipleAlignment alignment) {
		List<Atom[]> transformed = transformAtoms(alignment);
		return getRMSD(transformed);
	}
	/**
	 * Calculates the RMSD of all-to-all structure comparisons (distances),
	 * given a set of superimposed atoms.
	 * 
	 * @param transformed
	 * @return double RMSD
	 * @see #getRMSD(MultipleAlignment)
	 */
	private static double getRMSD(List<Atom[]> transformed) {

		double sumSqDist = 0;
		int comparisons = 0;
		
		for (int r1=0; r1<transformed.size(); r1++){
			for(int c=0;c<transformed.get(r1).length;c++) {
				Atom refAtom = transformed.get(r1)[c];
				if(refAtom == null) continue;

				double nonNullSqDist = 0;
				int nonNullLength = 0;
				for(int r2=r1+1;r2<transformed.size();r2++) {
					Atom atom = transformed.get(r2)[c];
					if(atom != null) {
						nonNullSqDist += Calc.getDistanceFast(refAtom, atom);
						nonNullLength++;
					}
				}
	
				if(nonNullLength > 0) {
					comparisons++;
					sumSqDist += nonNullSqDist/nonNullLength;
				}
			}
		}
		return Math.sqrt(sumSqDist/comparisons);
	}
	
	public static double getRefRMSD(MultipleAlignment alignment, int reference) {
		List<Atom[]> transformed = transformAtoms(alignment);
		return getRefRMSD(transformed,reference);
	}
	/**
	 * Calculates the average RMSD from all structures to a reference structure,
	 * given a set of superimposed atoms.<p>
	 * Complexity: T(n) = O(n), if n=number of structures.
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
	 * 
	 * @param transformed
	 * @param reference
	 * @return
	 */
	private static double getRefRMSD(List<Atom[]> transformed, int reference) {

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
	
	/**
	 * Calculates the average TMScore all the possible pairwise structure comparisons of the
	 * given MultipleAlignment. <p>
	 * Complexity: T(n) = O(n^2), if n=number of structures.
	 * 
	 * @param alignment
	 * @return double Average TMscore
	 * @throws StructureException
	 */
	public static double getAvgTMScore(MultipleAlignment alignment) throws StructureException {
		List<Atom[]> transformed = transformAtoms(alignment);
		List<Integer> lengths = new ArrayList<Integer>(alignment.size());
		for(Atom[] atoms : alignment.getEnsemble().getAtomArrays()) {
			lengths.add(atoms.length);
		}
		return getAvgTMScore(transformed,lengths);
	}
	/**
	 * Calculates the average TMScore all the possible pairwise structure comparisons of the
	 * given a set of superimposed Atoms and the original structure lengths.<p>
	 * Complexity: T(n) = O(n^2), if n=number of structures.
	 * 
	 * @param transformed
	 * @param lengths lengths of the structures in residue number
	 * @return double Average TMscore
	 * @throws StructureException
	 */
	private static double getAvgTMScore(List<Atom[]> transformed, List<Integer> lengths) throws StructureException {

		if(transformed.size() != lengths.size()) throw new IllegalArgumentException("Input sizes differ.");
		double sumTM = 0;
		int comparisons = 0;
		
		for (int r1=0; r1<transformed.size(); r1++){
			for(int r2=r1+1;r2<transformed.size();r2++) {
				int len = transformed.get(r1).length;
				//Remove nulls from both arrays
				Atom[] ref = new Atom[len];
				Atom[] aln = new Atom[len];
				int nonNullLen = 0;
				for(int c=0;c<len;c++) {
					if( transformed.get(r1)[c] != null && transformed.get(r2)[c] != null) {
						ref[nonNullLen] = transformed.get(r1)[c];
						aln[nonNullLen] = transformed.get(r2)[c];
						nonNullLen++;
					}
				}
				//Truncate nulls
				if(nonNullLen<len) {
					ref = Arrays.copyOf(ref, nonNullLen);
					aln = Arrays.copyOf(aln, nonNullLen);
				}
				sumTM += SVDSuperimposer.getTMScore(ref, aln, lengths.get(r1), lengths.get(r2));
				comparisons++;
			}
		}
		return sumTM/comparisons;
	}

	/**
	 * Calculates the average TMScore from all structures to a reference structure,
	 * given a set of superimposed atoms.<p>
	 * Complexity: T(n) = O(n), if n=number of structures.
	 * 
	 * @param alignment
	 * @param reference Index of the reference structure
	 * @return
	 * @throws StructureException 
	 */
	public static double getRefTMScore(MultipleAlignment alignment, int reference) throws StructureException {
		List<Atom[]> transformed = transformAtoms(alignment);
		List<Integer> lengths = new ArrayList<Integer>(alignment.size());
		for(Atom[] atoms : alignment.getEnsemble().getAtomArrays()) {
			lengths.add(atoms.length);
		}
		return getRefTMScore(transformed,lengths,reference);
	}
	/**
	 * Calculates the average TMScore from all structures to a reference structure,
	 * given a set of superimposed atoms.<p>
	 * Complexity: T(n) = O(n), if n=number of structures.
	 * 
	 * @param transformed Arrays of aligned atoms, after superposition
	 * @param lengths lengths of the full input structures
	 * @param reference Index of the reference structure
	 * @return
	 * @throws StructureException 
	 */
	private static double getRefTMScore(List<Atom[]> transformed, List<Integer> lengths, int reference) throws StructureException {

		if(transformed.size() != lengths.size()) throw new IllegalArgumentException("Input sizes differ");
		double sumTM = 0;
		int comparisons = 0;
		
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
			comparisons++;
		}
		return sumTM/comparisons;
	}
	
	/**
	 * Calculates the CEMC score, specific for the MultipleAlignment algorithm.
	 * The score function is modified from the original CEMC paper, making it
	 * continuous and differentiable.<p>
	 * Complexity: T(n) = O(n^2), if n=number of structures.
	 * 
	 * @param alignment
	 * @return
	 * @throws StructureException 
	 */
	public static double getCEMCScore(MultipleAlignment alignment) throws StructureException {
		//Transform Atoms
		List<Atom[]> transformed = transformAtoms(alignment);
		//Calculate d0
		int minLen = Integer.MAX_VALUE;
		for(Atom[] atoms : alignment.getEnsemble().getAtomArrays())
			if (atoms.length < minLen) minLen = atoms.length;
		double d0 =  1.24 * Math.cbrt((minLen) - 15.) - 1.8;	//d0 is calculated as in the TM-score
		//Calculate gaps
		int gaps = alignment.length()-alignment.getCoreLength();
		return getCEMCScore(transformed, d0, gaps);
	}
	
	private static double getCEMCScore(List<Atom[]> transformed, double d0, int gaps) throws StructureException {

		int length = transformed.get(0).length;
		double[] colDistances = new double[length];
		for (int i=0; i<length; i++) colDistances[i] = -1;
		double scoreMC = 0.0;
		
		//Calculate the average column distances
		for (int r1=0; r1<transformed.size(); r1++){
			for(int c=0;c<transformed.get(r1).length;c++) {
				Atom refAtom = transformed.get(r1)[c];
				if(refAtom == null) continue;

				double nonNullDist = 0;
				int nonNullLength = 0;
				for(int r2=r1+1;r2<transformed.size();r2++) {
					Atom atom = transformed.get(r2)[c];
					if(atom != null) {
						nonNullDist += Calc.getDistance(refAtom, atom);
						nonNullLength++;
					}
				}
				if(nonNullLength > 0) {
					if (colDistances[c] == -1) colDistances[c] = 0;
					colDistances[c] += nonNullDist/nonNullLength;
				}
			}
		}
		
		//Loop through all the columns
		for (int col=0; col<length; col++){
			if (colDistances[col]==-1) continue;
			double d1 = colDistances[col];
			double colScore = 40.0/(1+(d1*d1)/(d0*d0))-20.0;  //d0=5, M=20
			scoreMC += colScore;
		}
		//Apply the Gap penalty and return
		return scoreMC - gaps*10.0;
	}
}
