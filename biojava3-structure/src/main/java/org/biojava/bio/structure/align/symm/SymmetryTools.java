/**
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
 * @since 3.0.8
 */
package org.biojava.bio.structure.align.symm;


import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.helper.AlignTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.jama.Matrix;

/**
 * Tools related to symmetry-detection.
 */
public class SymmetryTools {

	// there won't be an instance of this
	private SymmetryTools(){}

//	public static void showMatrix(Matrix m, String string) {
//		ScaleableMatrixPanel smp = new ScaleableMatrixPanel();
//		JFrame frame = new JFrame();
//
//		smp.setMatrix((Matrix)m.clone());
//		//smp.getMatrixPanel().setScale(0.8f);
//
//		frame.setTitle(string);
//		frame.getContentPane().add(smp);
//		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
//		frame.pack();
//		frame.setVisible(true);
//	}

	/**
	 * Returns the "reset value" for graying out the main diagonal. If we're blanking out the main diagonal, this value is always Integer.MIN_VALUE. This is possible if {@code gradientPolyCoeff = {Integer.MIN_VALUE}} and {@code gradientExpCoeff = 0}.
	 * @param unpenalizedScore
	 * @param nResFromMainDiag
	 * @param gradientPolyCoeff
	 * @param gradientExpCoeff
	 * @return
	 */
	private static double getResetVal(double unpenalizedScore, double nResFromMainDiag, double[] gradientPolyCoeff, double gradientExpCoeff) {
		if (unpenalizedScore == Double.NaN) return 0; // what else?
		double updateVal = unpenalizedScore; // notice that we can actually return a positive value if this is high enough
		updateVal -= gradientExpCoeff * Math.pow(Math.E, -nResFromMainDiag);
		for (int p = 0; p < gradientPolyCoeff.length; p++) {
			updateVal -= gradientPolyCoeff[gradientPolyCoeff.length-1-p] * Math.pow(nResFromMainDiag, -p);
		}
		//System.out.println("For unpenalized " + unpenalizedScore + " and " + nResFromMainDiag + " residues from diagonal: " + (updateVal-unpenalizedScore));
		return updateVal;
	}

	public static Matrix grayOutCEOrig(Atom[] ca2, int rows, int cols,
			CECalculator calculator, Matrix origM, int blankWindowSize, double[] gradientPolyCoeff, double gradientExpCoeff) {

		if ( origM == null)
			origM =   new Matrix( calculator.getMatMatrix());

		// symmetry hack, disable main diagonal

		for ( int i = 0 ; i< rows; i++){
			for ( int j = 0 ; j < cols ; j++){
				int diff = Math.abs(i-j);

				double resetVal = getResetVal(origM.get(i, j), diff, gradientPolyCoeff, gradientExpCoeff);

				if ( diff < blankWindowSize ){
					origM.set(i,j, origM.get(i, j) + resetVal);

				}
				int diff2 = Math.abs(i-(j-ca2.length/2)); // other side

				double resetVal2 = getResetVal(origM.get(i, j), diff2, gradientPolyCoeff, gradientExpCoeff);

				if ( diff2 < blankWindowSize ){
					origM.set(i,j, origM.get(i, j) + resetVal2);

				}
			}
		}
		return origM;
	}

	public static Matrix grayOutPreviousAlignment(AFPChain afpChain, Atom[] ca2,
			int rows, int cols, CECalculator calculator, Matrix max, int blankWindowSize, double[] gradientPolyCoeff, double gradientExpCoeff) {

		max =  grayOutCEOrig(ca2, rows, cols, calculator, max,  blankWindowSize, gradientPolyCoeff, gradientExpCoeff);

		double[][] dist1 = calculator.getDist1();
		double[][] dist2 = calculator.getDist2();

		int[][][] optAln = afpChain.getOptAln();
		int blockNum = afpChain.getBlockNum();

		int[] optLen = afpChain.getOptLen();

		// ca2 is circularly permutated
		int breakPoint = ca2.length / 2;
		for (int bk = 0; bk < blockNum; bk++)       {

			for ( int i=0;i< optLen[bk];i++){
				int pos1 = optAln[bk][0][i];
				int pos2 = optAln[bk][1][i];

				int dist = blankWindowSize/2 ;
				int start1 = Math.max(pos1-dist,0);
				int start2 = Math.max(pos2-dist,0);
				int end1 = Math.min(pos1+dist, rows-1);
				int end2 = Math.min(pos2+dist, cols-1);

				for ( int i1 = start1; i1< end1 ; i1++){

					for ( int k=0; k < blankWindowSize/2 ; k ++){
						if ( i1-k >= 0) {
							double resetVal = getResetVal(max.get(i1-k, i1-k), 0, gradientPolyCoeff, gradientExpCoeff);
							dist1[i1-k][i1-k] = resetVal;
						} else if ( i1+k < rows) {
							double resetVal = getResetVal(max.get(i1+k, i1+k), 0, gradientPolyCoeff, gradientExpCoeff);
							dist1[i1+k][i1+k] = resetVal;
						}

					}

					for ( int j2 = start2 ; j2 < end2 ; j2++){
						double resetVal = getResetVal(max.get(i1, j2), Math.abs(i1-j2), gradientPolyCoeff, gradientExpCoeff);
						max.set(i1,j2,resetVal);
						if ( j2 < breakPoint) {
							double resetVal2 = getResetVal(max.get(i1, j2+breakPoint), Math.abs(i1-(j2+breakPoint)), gradientPolyCoeff, gradientExpCoeff);
							max.set(i1,j2+breakPoint,resetVal2);
						} else {
							double resetVal2 = getResetVal(max.get(i1, j2-breakPoint), Math.abs(i1-(j2-breakPoint)), gradientPolyCoeff, gradientExpCoeff);
							max.set(i1,j2-breakPoint,resetVal2);
						}
						for ( int k=0; k <blankWindowSize/2 ; k ++){							
							if ( j2-k >=0) {
								double resetVal2 = getResetVal(max.get(j2-k, j2-k), 0, gradientPolyCoeff, gradientExpCoeff);
								dist2[j2-k][j2-k] = resetVal2;
							} else if ( j2+k < cols) {
								double resetVal2 = getResetVal(max.get(j2+k, j2+k), 0, gradientPolyCoeff, gradientExpCoeff);
								dist2[j2+k][j2+k] = resetVal2;
							}
						}
					}
				}

			}
		}
		calculator.setDist1(dist1);
		calculator.setDist2(dist2);
		return max;

	}

	public Matrix  getDkMatrix(Atom[] ca1, Atom[] ca2,int fragmentLength,
			double[] dist1, double[] dist2, int rows, int cols) {
		Matrix diffDistMax =  Matrix.identity(ca1.length, ca2.length);

		for ( int i = 0 ; i< rows; i++){
			double score1 = 0;
			for ( int x=0 ; x < fragmentLength ; x++){
				score1 += dist1[i+x];
			}
			for ( int j = 0 ; j < cols ; j++){
				double score2 = 0;
				for ( int y=0 ; y < fragmentLength ; y++){
					score2 += dist2[j+y];
				}

				// if the intramolecular distances are very similar
				// the two scores should be similar, i.e. the difference is close to 0
				diffDistMax.set(i,j, Math.abs(score1-score2));
			}
		}


		// symmetry hack, disable main diagonal

		for ( int i = 0 ; i< rows; i++){
			for ( int j = 0 ; j < cols ; j++){
				int diff = Math.abs(i-j);

				if ( diff < 15 ){
					diffDistMax.set(i,j, 99);
				}
				int diff2 = Math.abs(i-(j-ca2.length/2));
				if ( diff2 < 15 ){
					diffDistMax.set(i,j, 99);
				}
			}
		}
		return diffDistMax;

	}

	public static Matrix blankOutPreviousAlignment(AFPChain afpChain, Atom[] ca2,
			int rows, int cols, CECalculator calculator, Matrix max, int blankWindowSize) {
		return grayOutPreviousAlignment(afpChain, ca2, rows, cols, calculator, max, blankWindowSize, new double[] {Integer.MIN_VALUE}, 0.0);

	}

	public static Matrix blankOutCEOrig(Atom[] ca2, int rows, int cols,
			CECalculator calculator, Matrix origM, int blankWindowSize) {
		return grayOutCEOrig(ca2, rows, cols, calculator, origM, blankWindowSize, new double[] {Integer.MIN_VALUE}, 0.0);
	}

	public static Atom[] cloneAtoms(Atom[] ca2) throws StructureException{
		// we don't want to rotate input atoms, do we?
		Atom[] ca2clone = new Atom[ca2.length];

		int pos = 0;
		for (Atom a : ca2){
			Group g = (Group) a.getGroup().clone(); // works because each group has only a CA atom

			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos++;
		}

		return ca2clone;
	}

	public static Matrix getDkMatrix(Atom[] ca1, Atom[] ca2, int k, int fragmentLength) {
		double[] dist1 = AlignTools.getDiagonalAtK(ca1, k);

		double[] dist2 = AlignTools.getDiagonalAtK(ca2, k);

		int rows = ca1.length - fragmentLength - k + 1;
		int cols = ca2.length - fragmentLength - k + 1;

		// Matrix that tracks similarity of a fragment of length fragmentLength
		// starting a position i,j.

		Matrix m2 = new Matrix(rows,cols); 

		for ( int i = 0 ; i< rows; i++){
			double score1 = 0;
			for ( int x=0 ; x < fragmentLength ; x++){
				score1 += dist1[i+x];
			}
			for ( int j = 0 ; j < cols ; j++){
				double score2 = 0;
				for ( int y=0 ; y < fragmentLength ; y++){
					score2 += dist2[j+y];
				}	

				// if the intramolecular distances are very similar
				// the two scores should be similar, i.e. the difference is close to 0
				m2.set(i,j, Math.abs(score1-score2));
			}
		}
		return m2;
	}


	public static boolean[][] blankOutBreakFlag(AFPChain afpChain,
			Atom[] ca2, int rows, int cols, CECalculator calculator,
			boolean[][] breakFlag, int blankWindowSize) {


		int[][][] optAln = afpChain.getOptAln();
		int blockNum = afpChain.getBlockNum();

		int[] optLen = afpChain.getOptLen();

		// ca2 is circularly permutated at this point.
		int breakPoint = ca2.length / 2;

		for(int bk = 0; bk < blockNum; bk ++)       {

			//Matrix m= afpChain.getBlockRotationMatrix()[bk];
			//Atom shift = afpChain.getBlockShiftVector()[bk];
			for ( int i=0;i< optLen[bk];i++){
				int pos1 = optAln[bk][0][i];
				int pos2 = optAln[bk][1][i];
				// blank out area around these positions...

				int dist = blankWindowSize ;
				int start1 = Math.max(pos1-dist,0);
				int start2 = Math.max(pos2-dist,0);
				int end1 = Math.min(pos1+dist, rows-1);
				int end2 = Math.min(pos2+dist, cols-1);

				//System.out.println(pos1 + "  " + pos2 + " " + start1 + " " + end1 + " " + start2 + " " + end2);

				for ( int i1 = start1; i1< end1 ; i1++){

					for ( int j2 = start2 ; j2 < end2 ; j2++){
						//System.out.println(i1 + " " + j2 + " (***)");
						breakFlag[i1][j2] = true;
						if ( j2 < breakPoint) {
							breakFlag[i1][j2+ breakPoint ] = true;
						}
					}
				}

			}
		}

		return breakFlag;
	}

	/**
	 * Returns the <em>magnitude</em> of the angle between the first and second blocks of {@code afpChain}, measured in degrees. This is always a positive value (unsigned).
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @return
	 */
	public static double getAngle(AFPChain afpChain, Atom[] ca1, Atom[] ca2) {
		Matrix rotation = afpChain.getBlockRotationMatrix()[0];
		return Math.acos(rotation.trace() - 1) * 180/Math.PI;
	}
	

}
