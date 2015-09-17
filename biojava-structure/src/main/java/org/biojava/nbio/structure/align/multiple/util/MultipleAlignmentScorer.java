package org.biojava.nbio.structure.align.multiple.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * Utility class for calculating common scores of {@link MultipleAlignment}s.
 * 
 * @author Spencer Bliven
 * @author Aleix Lafita
 * @since 4.1.0
 *
 */
public class MultipleAlignmentScorer {
	
	//Names for commonly used scores
	public static final String PROBABILITY = "Probability";
	public static final String CE_SCORE = "CE-score";
	public static final String RMSD = "RMSD";
	public static final String AVGTM_SCORE = "AvgTM-score";
	public static final String MC_SCORE = "MC-score";
	public static final String REF_RMSD = "Ref-RMSD";
	public static final String REFTM_SCORE = "RefTM-score";

	/**
	 * Calculates and puts the RMSD and the average TM-Score 
	 * of the MultipleAlignment.
	 * 
	 * @param alignment
	 * @throws StructureException
	 * @see #getAvgTMScore(MultipleAlignment)
	 * @see #getRMSD(MultipleAlignment)
	 */
	public static void calculateScores(MultipleAlignment alignment) 
			throws StructureException {
		
		//Put RMSD
		List<Atom[]> trans = MultipleAlignmentTools.transformAtoms(alignment);
		alignment.putScore(RMSD, getRMSD(trans));
		
		//Put AvgTM-Score
		List<Integer> lengths = new ArrayList<Integer>(alignment.size());
		for(Atom[] atoms : alignment.getAtomArrays()) {
			lengths.add(atoms.length);
		}
		alignment.putScore(AVGTM_SCORE, getAvgTMScore(trans,lengths));
	}
	
	/**
	 * Calculates the RMSD of all-to-all structure comparisons (distances) of 
	 * the given MultipleAlignment. <p>
	 * Complexity: T(n,l) = O(l*n^2), if n=number of structures and 
	 * l=alignment length.
	 * <p>
	 * The formula used is just the sqroot of the average distance
	 * of all possible pairs of atoms. Thus, for R structures
	 * aligned over C columns without gaps, we have
	 * <pre>RMSD = \sqrt{1/(C*(R*(R-1)/2)) * \sum_{r1=1}^{R-1} 
	 * \sum_{r2=r1+1}^{R-1} \sum_{j=0}^{C-1} (atom[r1][c]-atom[r2][c])^2}
	 * </pre>
	 * 
	 * @param alignment
	 * @return double RMSD
	 */
	public static double getRMSD(MultipleAlignment alignment) {
		List<Atom[]> trans = MultipleAlignmentTools.transformAtoms(alignment);
		return getRMSD(trans);
	}
	
	/**
	 * Calculates the RMSD of all-to-all structure comparisons (distances),
	 * given a set of superimposed atoms.
	 * 
	 * @param transformed
	 * @return double RMSD
	 * @see #getRMSD(MultipleAlignment)
	 */
	public static double getRMSD(List<Atom[]> transformed) {

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
	
	/**
	 * /**
	 * Calculates the average RMSD from all structures to a reference s
	 * tructure, given a set of superimposed atoms.<p>
	 * Complexity: T(n,l) = O(l*n), if n=number of structures and 
	 * l=alignment length.
	 * <p>
	 * For ungapped alignments, this is just the sqroot of the average 
	 * distance from an atom to the aligned atom from the reference. 
	 * Thus, for R structures aligned over C columns (with structure 0 as 
	 * the reference), we have: 
	 * <pre>RefRMSD = \sqrt{ 1/(C*(R-1)) * \sum_{r=1}^{R-1} \sum_{j=0}^{C-1} 
	 * (atom[1][c]-atom[r][c])^2 }</pre>
	 * <p>
	 * For gapped alignments, null atoms are omitted from consideration, 
	 * so that the RMSD is the average over all columns with non-null 
	 * reference of the average RMSD within the non-null elements of the 
	 * column.
	 * 
	 * @param alignment MultipleAlignment
	 * @param ref reference structure index
	 * @return
	 */
	public static double getRefRMSD(MultipleAlignment alignment, int ref){
		List<Atom[]> trans = MultipleAlignmentTools.transformAtoms(alignment);
		return getRefRMSD(trans, ref);
	}
	
	/**
	 * Calculates the average RMSD from all structures to a reference s
	 * tructure, given a set of superimposed atoms.<p>
	 * Complexity: T(n,l) = O(l*n), if n=number of structures and 
	 * l=alignment length.
	 * <p>
	 * For ungapped alignments, this is just the sqroot of the average 
	 * distance from an atom to the aligned atom from the reference. 
	 * Thus, for R structures aligned over C columns (with structure 0 as 
	 * the reference), we have: 
	 * <pre>RefRMSD = \sqrt{ 1/(C*(R-1)) * \sum_{r=1}^{R-1} \sum_{j=0}^{C-1} 
	 * (atom[1][c]-atom[r][c])^2 }</pre>
	 * <p>
	 * For gapped alignments, null atoms are omitted from consideration, 
	 * so that the RMSD is the average over all columns with non-null 
	 * reference of the average RMSD within the non-null elements of the 
	 * column.
	 * 
	 * @param transformed
	 * @param reference
	 * @return
	 */
	public static double getRefRMSD(List<Atom[]> transformed, int reference) {

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
	 * Calculates the average TMScore of all the possible pairwise structure 
	 * comparisons of the given alignment. <p>
	 * Complexity: T(n,l) = O(l*n^2), if n=number of structures and 
	 * l=alignment length.
	 * 
	 * @param alignment
	 * @return double Average TMscore
	 * @throws StructureException
	 */
	public static double getAvgTMScore(MultipleAlignment alignment) 
			throws StructureException {
		
		List<Atom[]> trans = MultipleAlignmentTools.transformAtoms(alignment);
		
		List<Integer> lengths = new ArrayList<Integer>(alignment.size());
		for(Atom[] atoms : alignment.getAtomArrays()){
			lengths.add(atoms.length);
		}
		return getAvgTMScore(trans, lengths);
	}
	/**
	 * Calculates the average TMScore all the possible pairwise structure 
	 * comparisons of the given a set of superimposed Atoms and the original 
	 * structure lengths.<p>
	 * Complexity: T(n,l) = O(l*n^2), if n=number of structures and 
	 * l=alignment length.
	 * 
	 * @param transformed aligned Atoms transformed
	 * @param lengths lengths of the structures in residue number
	 * @return double Average TMscore
	 * @throws StructureException
	 */
	public static double getAvgTMScore(List<Atom[]> transformed, 
			List<Integer> lengths) throws StructureException {

		if(transformed.size() != lengths.size()) 
			throw new IllegalArgumentException("Input sizes differ.");
		
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
					if(transformed.get(r1)[c] != null 
							&& transformed.get(r2)[c] != null) {
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
				sumTM += SVDSuperimposer.getTMScore(ref, aln, 
						lengths.get(r1), lengths.get(r2));
				comparisons++;
			}
		}
		return sumTM/comparisons;
	}

	/**
	 * Calculates the average TMScore from all structures to a reference 
	 * structure, given a set of superimposed atoms.<p>
	 * Complexity: T(n,l) = O(l*n), if n=number of structures and 
	 * l=alignment length.
	 * 
	 * @param alignment
	 * @param reference Index of the reference structure
	 * @return
	 * @throws StructureException 
	 */
	public static double getRefTMScore(MultipleAlignment alignment, int ref) 
			throws StructureException {
		
		List<Atom[]> trans = MultipleAlignmentTools.transformAtoms(alignment);
		
		List<Integer> lengths = new ArrayList<Integer>(alignment.size());
		for(Atom[] atoms : alignment.getAtomArrays()) {
			lengths.add(atoms.length);
		}
		return getRefTMScore(trans,lengths,ref);
	}
	/**
	 * Calculates the average TMScore from all structures to a reference 
	 * structure, given a set of superimposed atoms.<p>
	 * Complexity: T(n,l) = O(l*n^2), if n=number of structures and 
	 * l=alignment length.
	 * 
	 * @param transformed Arrays of aligned atoms, after superposition
	 * @param lengths lengths of the full input structures
	 * @param reference Index of the reference structure
	 * @return
	 * @throws StructureException 
	 */
	public static double getRefTMScore(List<Atom[]> transformed, 
			List<Integer> lengths, int reference) throws StructureException {

		if(transformed.size() != lengths.size()) 
			throw new IllegalArgumentException("Input sizes differ");
		
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
				if(transformed.get(reference)[c] != null 
						&& transformed.get(r)[c] != null) {
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
			sumTM += SVDSuperimposer.getTMScore(ref, aln, 
					lengths.get(reference), lengths.get(r));
			comparisons++;
		}
		return sumTM/comparisons;
	}
	
	/**
	 * Calculates the MC score, specific for the MultipleAlignment algorithm.
	 * The score function is modified from the original CEMC paper, making it
	 * continuous and differentiable.
	 * <p>
	 * The maximum score of a match is 20, and the penalties for gaps are part 
	 * of the input. The half-score distance, d0, is chosen as in the TM-score.
	 * <p>
	 * Complexity: T(n,l) = O(l*n^2), if n=number of structures and 
	 * l=alignment length.
	 * 
	 * @param alignment
	 * @param gapOpen penalty for gap opening (reasonable values are in the range
	 *  (1.0-20.0)
	 * @param gapExtension penalty for extending a gap (reasonable values are 
	 * in the range (0.5-10.0)
	 * @param dCutoff the distance cutoff
	 * @return the value of the score
	 * @throws StructureException 
	 */
	public static double getMCScore(MultipleAlignment alignment, 
			double gapOpen, double gapExtension, double dCutoff) 
					throws StructureException {
		
		List<Atom[]> trans = MultipleAlignmentTools.transformAtoms(alignment);
		
		//Calculate d0: same as the one in TM score
		int minLen = Integer.MAX_VALUE;
		for(Atom[] atoms : alignment.getAtomArrays())
			if (atoms.length < minLen) minLen = atoms.length;
		double d0 =  1.24 * Math.cbrt((minLen) - 15.) - 1.8;
		
		//Calculate the distance cutoff penalization
		double A = 20.0/(1+(dCutoff*dCutoff)/(d0*d0));
		
		return getMCScore(trans, d0, gapOpen, gapExtension, A);
	}
	
	/**
	 * Calculates the MC score, specific for the MultipleAlignment algorithm.
	 * The score function is modified from the original CEMC paper, making it
	 * continuous and differentiable.
	 * <p>
	 * The maximum score of a match is 20, and the penalties for gaps are part
	 * of the input.
	 * <p>
	 * Complexity: T(n,l) = O(l*n^2), if n=number of structures and 
	 * l=alignment length.
	 * 
	 * @param transformed List of transformed Atom arrays
	 * @param d0 parameter for the half-score distance
	 * @param gapOpen penalty for gap opening (reasonable values are in the 
	 * range (1.0-20.0)
	 * @param gapExtension penalty for extending a gap (reasonable values 
	 * are in the range (0.5-10.0)
	 * @param A the distance cutoff penalization
	 * @return the value of the score
	 * @throws StructureException 
	 */
	private static double getMCScore(List<Atom[]> trans, double d0, 
			double gapOpen, double gapExtension, double A) 
					throws StructureException {

		int size = trans.size();
		int length = trans.get(0).length;
		Matrix residueDistances = new Matrix(size,length,-1);
		double scoreMC = 0.0;
		int openGaps = 0;
		int extensionGaps = 0;
		
		//Calculate the average residue distances
		for (int r1=0; r1<size; r1++){
			boolean gapped = false;
			for(int c=0;c<trans.get(r1).length;c++) {
				Atom refAtom = trans.get(r1)[c];
				//Calculate the gap extension and opening on the fly
				if(refAtom == null) {
					if (gapped) extensionGaps++;
					else {
						gapped = true;
						openGaps++;
					}
					continue;
				} else gapped = false;

				for(int r2=r1+1;r2<size;r2++) {
					Atom atom = trans.get(r2)[c];
					if(atom != null) {
						double distance = Calc.getDistance(refAtom, atom);
						if (residueDistances.get(r1, c) == -1) {
							residueDistances.set(r1, c, 1+distance);
						} else {
							residueDistances.set(r1, c, 
									residueDistances.get(r1, c)+distance);
						}
						if (residueDistances.get(r2, c) == -1) {
							residueDistances.set(r2, c, 1+distance);
						} else {
							residueDistances.set(r2, c, 
									residueDistances.get(r2, c)+distance);
						}
					}
				}
			}
		}
		
		for(int c=0;c<length;c++) {
			int nonNullRes = 0;
			for(int r=0;r<size;r++) {
				if (residueDistances.get(r, c) != -1) nonNullRes++;
			}
			for(int r=0;r<size;r++) {
				if (residueDistances.get(r, c) != -1) {
					residueDistances.set(r, c, 
							residueDistances.get(r, c)/nonNullRes);
				}
			}
		}
		
		//Sum all the aligned residue scores
		for(int row=0;row<size;row++) {
			for (int col=0; col<length; col++){
				if (residueDistances.get(row,col)==-1) continue;
				double d1 = residueDistances.get(row,col);
				double resScore = 20.0/(1+(d1*d1)/(d0*d0));
				scoreMC += resScore - A;
			}
		}
		
		//Apply the Gap penalty and return
		return scoreMC - (openGaps*gapOpen + extensionGaps*gapExtension);
	}
}
