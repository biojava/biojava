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
 * Created on Mar 9, 2010
 * Author: Spencer Bliven
 *
 */

package org.biojava.nbio.structure.align.ce;


import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.nbio.structure.align.util.ConfigurationException;
import org.biojava.nbio.structure.geometry.Matrices;
import org.biojava.nbio.structure.geometry.SuperPositions;
import org.biojava.nbio.structure.jama.Matrix;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Matrix4d;

/**
 * A wrapper for {@link CeMain} which sets default parameters to be appropriate for finding
 * circular permutations.
 * <p>
 * A circular permutation consists of a single cleavage point and rearrangement
 * between two structures, for example:
 * <pre>
 * ABCDEFG
 * DEFGABC
 * </pre>
 * @author Spencer Bliven.
 *
 */
public class CeCPMain extends CeMain {
	private static boolean debug = false;

	public static final String algorithmName = "jCE Circular Permutation";

	/**
	 *  version history:
	 *  1.5 - Added more parameters to the command line, including -maxOptRMSD
	 *  1.4 - Added DuplicationHint parameter & default to duplicating the shorter chain
	 *  1.3 - Short CPs are now discarded
	 *  1.2 - now supports check AlignmentTools.isSequentialAlignment. XML protocol
	 *  1.1 - skipped, (trying to avoid confusion with jfatcat in all vs. all comparisons)
	 *  1.0 - initial release
	 */
	public static final String version = "1.5";

	public CeCPMain(){
		super();
		params = new CECPParameters();
	}

	@Override
	public String getAlgorithmName() {
		return CeCPMain.algorithmName;
	}

	@Override
	public String getVersion() {
		return CeCPMain.version;
	}

	public static void main(String[] args) throws ConfigurationException {
		CeCPUserArgumentProcessor processor = new CeCPUserArgumentProcessor(); //Responsible for creating a CeCPMain instance
		processor.process(args);
	}


	/**
	 * Aligns ca1 and ca2 using a heuristic to check for CPs.
	 * <p>
	 * Aligns ca1 against a doubled ca2, then cleans up the alignment.
	 * @param ca1
	 * @param ca2
	 * @param param
	 * @return the alignment, possibly containing a CP.
	 * @throws StructureException
	 */
	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object param) throws StructureException{
		if ( ! (param instanceof CECPParameters))
			throw new IllegalArgumentException("CE algorithm needs an object of call CeParameters as argument.");

		CECPParameters cpparams = (CECPParameters) param;
		this.params = cpparams;

		boolean duplicateRight;

		switch( cpparams.getDuplicationHint() ) {
		case LEFT:
			duplicateRight = false;
			break;
		case RIGHT:
			duplicateRight = true;
			break;
		case SHORTER:
			duplicateRight = ca1.length >= ca2.length;
			break;
		default:
			duplicateRight = true;
		}


		if( duplicateRight ) {
			return alignRight(ca1, ca2, cpparams);
		} else {
			if(debug) {
				System.out.println("Swapping alignment order.");
			}
			AFPChain afpChain = this.alignRight(ca2, ca1, cpparams);
			return invertAlignment(afpChain);
		}
	}

	/**
	 * Aligns the structures, duplicating ca2 regardless of
	 * {@link CECPParameters.getDuplicationHint() param.getDuplicationHint}.
	 * @param ca1
	 * @param ca2
	 * @param cpparams
	 * @return
	 * @throws StructureException
	 */
	private AFPChain alignRight(Atom[] ca1, Atom[] ca2, CECPParameters cpparams)
			throws StructureException {
		long startTime = System.currentTimeMillis();

		Atom[] ca2m = StructureTools.duplicateCA2(ca2);

		if(debug) {
			System.out.format("Duplicating ca2 took %s ms\n",System.currentTimeMillis()-startTime);
			startTime = System.currentTimeMillis();
		}

		// Do alignment
		AFPChain afpChain = super.align(ca1, ca2m,params);

		// since the process of creating ca2m strips the name info away, set it explicitely
		try {
			afpChain.setName2(ca2[0].getGroup().getChain().getStructure().getName());
		} catch( Exception e) {}

		if(debug) {
			System.out.format("Running %dx2*%d alignment took %s ms\n",ca1.length,ca2.length,System.currentTimeMillis()-startTime);
			startTime = System.currentTimeMillis();
		}
		afpChain = postProcessAlignment(afpChain, ca1, ca2m, calculator, cpparams);

		if(debug) {
			System.out.format("Finding CP point took %s ms\n",System.currentTimeMillis()-startTime);
			startTime = System.currentTimeMillis();
		}

		return afpChain;
	}



	/** Circular permutation specific code to be run after the standard CE alignment
	 *
	 * @param afpChain The finished alignment
	 * @param ca1 CA atoms of the first protein
	 * @param ca2m A duplicated copy of the second protein
	 * @param calculator The CECalculator used to create afpChain
	 * @throws StructureException
	 */
	public static AFPChain postProcessAlignment(AFPChain afpChain, Atom[] ca1, Atom[] ca2m,CECalculator calculator ) throws StructureException{
		return postProcessAlignment(afpChain, ca1, ca2m, calculator, null);
	}

	/** Circular permutation specific code to be run after the standard CE alignment
	 *
	 * @param afpChain The finished alignment
	 * @param ca1 CA atoms of the first protein
	 * @param ca2m A duplicated copy of the second protein
	 * @param calculator The CECalculator used to create afpChain
	 * @param param Parameters
	 * @throws StructureException
	 */
	public static AFPChain postProcessAlignment(AFPChain afpChain, Atom[] ca1, Atom[] ca2m,CECalculator calculator, CECPParameters param ) throws StructureException{

		// remove bottom half of the matrix
		Matrix doubledMatrix = afpChain.getDistanceMatrix();

		// the matrix can be null if the alignment is too short.
		if ( doubledMatrix != null ) {
			assert(doubledMatrix.getRowDimension() == ca1.length);
			assert(doubledMatrix.getColumnDimension() == ca2m.length);

			Matrix singleMatrix = doubledMatrix.getMatrix(0, ca1.length-1, 0, (ca2m.length/2)-1);
			assert(singleMatrix.getRowDimension() == ca1.length);
			assert(singleMatrix.getColumnDimension() == (ca2m.length/2));

			afpChain.setDistanceMatrix(singleMatrix);
		}
		// Check for circular permutations
		int alignLen = afpChain.getOptLength();
		if ( alignLen > 0) {
			afpChain = filterDuplicateAFPs(afpChain,calculator,ca1,ca2m,param);
		}
		return afpChain;
	}

	/**
	 * Swaps the order of structures in an AFPChain
	 * @param a
	 * @return
	 */
	public AFPChain invertAlignment(AFPChain a) {
		String name1 = a.getName1();
		String name2 = a.getName2();
		a.setName1(name2);
		a.setName2(name1);

		int len1 = a.getCa1Length();
		a.setCa1Length( a.getCa2Length() );
		a.setCa2Length( len1 );

		int beg1 = a.getAlnbeg1();
		a.setAlnbeg1(a.getAlnbeg2());
		a.setAlnbeg2(beg1);

		char[] alnseq1 = a.getAlnseq1();
		a.setAlnseq1(a.getAlnseq2());
		a.setAlnseq2(alnseq1);

		Matrix distab1 = a.getDisTable1();
		a.setDisTable1(a.getDisTable2());
		a.setDisTable2(distab1);

		int[] focusRes1 = a.getFocusRes1();
		a.setFocusRes1(a.getFocusRes2());
		a.setFocusRes2(focusRes1);

		//What are aftIndex and befIndex used for? How are they indexed?
		//a.getAfpAftIndex()


		String[][][] pdbAln = a.getPdbAln();
		if( pdbAln != null) {
			for(int block = 0; block < a.getBlockNum(); block++) {
				String[] paln1 = pdbAln[block][0];
				pdbAln[block][0] = pdbAln[block][1];
				pdbAln[block][1] = paln1;
			}
		}

		int[][][] optAln = a.getOptAln();
		if( optAln != null ) {
			for(int block = 0; block < a.getBlockNum(); block++) {
				int[] aln1 = optAln[block][0];
				optAln[block][0] = optAln[block][1];
				optAln[block][1] = aln1;
			}
		}
		a.setOptAln(optAln); // triggers invalidate()

		Matrix distmat = a.getDistanceMatrix();
		if(distmat != null)
			a.setDistanceMatrix(distmat.transpose());


		// invert the rotation matrices
		Matrix[] blockRotMat = a.getBlockRotationMatrix();
		Atom[] shiftVec = a.getBlockShiftVector();
		if( blockRotMat != null) {
			for(int block = 0; block < a.getBlockNum(); block++) {
				if(blockRotMat[block] != null) {
					// if y=x*A+b, then x=y*inv(A)-b*inv(A)
					blockRotMat[block] = blockRotMat[block].inverse();

					Calc.rotate(shiftVec[block],blockRotMat[block]);
					shiftVec[block] = Calc.invert(shiftVec[block]);
				}
			}
		}

		return a;
	}

	/**
	 * Takes as input an AFPChain where ca2 has been artificially duplicated.
	 * This raises the possibility that some residues of ca2 will appear in
	 * multiple AFPs. This method filters out duplicates and makes sure that
	 * all AFPs are numbered relative to the original ca2.
	 *
	 * <p>The current version chooses a CP site such that the length of the
	 * alignment is maximized.
	 *
	 * <p>This method does <i>not</i> update scores to reflect the filtered alignment.
	 * It <i>does</i> update the RMSD and superposition.
	 *
	 * @param afpChain The alignment between ca1 and ca2-ca2. Blindly assumes
	 *  that ca2 has been duplicated.
	 * @return A new AFPChain consisting of ca1 to ca2, with each residue in
	 *  at most 1 AFP.
	 * @throws StructureException
	 */
	public static AFPChain filterDuplicateAFPs(AFPChain afpChain, CECalculator ceCalc, Atom[] ca1, Atom[] ca2duplicated) throws StructureException {
		return filterDuplicateAFPs(afpChain, ceCalc, ca1, ca2duplicated, null);
	}
	public static AFPChain filterDuplicateAFPs(AFPChain afpChain, CECalculator ceCalc,
			Atom[] ca1, Atom[] ca2duplicated, CECPParameters params) throws StructureException {
		AFPChain newAFPChain = new AFPChain(afpChain);

		if(params == null)
			params = new CECPParameters();

		final int minCPlength = params.getMinCPLength();


		int ca2len = afpChain.getCa2Length()/2;
		newAFPChain.setCa2Length(ca2len);

		// Fix optimal alignment
		int[][][] optAln = afpChain.getOptAln();
		int[] optLen = afpChain.getOptLen();
		int alignLen = afpChain.getOptLength();
		if (alignLen < 1) return newAFPChain;

		assert(afpChain.getBlockNum() == 1); // Assume that CE returns just one block


		// Determine the region where ca2 and ca2' overlap

		// The bounds of the alignment wrt ca2-ca2'
		int nStart = optAln[0][1][0]; //alignment N-terminal
		int cEnd = optAln[0][1][alignLen-1]; // alignment C-terminal
		// overlap is between nStart and cEnd


		int firstRes = nStart; // start res number after trimming
		int lastRes = nStart+ca2len;  // last res number after trimming
		if(nStart >= ca2len || cEnd < ca2len) { // no circular permutation
			firstRes=nStart;
			lastRes=cEnd;
		} else {
			// Rule: maximize the length of the alignment

			int overlapLength = cEnd+1 - nStart - ca2len;
			if(overlapLength <= 0) {
				// no overlap!

				CPRange minCP = calculateMinCP(optAln[0][1], alignLen, ca2len, minCPlength);

				firstRes=nStart;
				lastRes=cEnd;

				// Remove short blocks
				if(firstRes > minCP.n) {
					firstRes = ca2len;

					if(debug) {
						System.out.format("Discarding n-terminal block as too " +
								"short (%d residues, needs %d)\n",
								minCP.mid, minCPlength);
					}
				}

				if(lastRes < minCP.c) {
					lastRes = ca2len-1;

					if(debug) {
						System.out.format("Discarding c-terminal block as too " +
								"short (%d residues, needs %d)\n",
								optLen[0] - minCP.mid, minCPlength);
					}
				}

			}
			else {
				// overlap!

				CutPoint cp = calculateCutPoint(optAln[0][1], nStart, cEnd,
						overlapLength, alignLen, minCPlength, ca2len, firstRes);

				// Adjust alignment length for trimming
				//alignLen -= cp.numResiduesCut; //TODO inaccurate

				firstRes = cp.firstRes;
				lastRes = cp.lastRes;

				//TODO Now have CP site, and could do a NxM alignment for further optimization.
				// For now, does not appear to be worth the 50% increase in time

				//TODO Bug: scores need to be recalculated
			}
		}



		// Fix numbering:
		// First, split up the atoms into left and right blocks
		List< ResiduePair > left = new ArrayList<ResiduePair>(); // residues from left of duplication
		List< ResiduePair > right = new ArrayList<ResiduePair>(); // residues from right of duplication

		for(int i=0;i<optLen[0];i++) {
			if( optAln[0][1][i] >= firstRes && optAln[0][1][i] <= lastRes ) { // not trimmed
				if(optAln[0][1][i] < ca2len) { // in first half of ca2
					left.add(new ResiduePair(optAln[0][0][i],optAln[0][1][i]));
				}
				else {
					right.add(new ResiduePair(optAln[0][0][i],optAln[0][1][i]-ca2len));
				}
			}
		}
		//assert(left.size()+right.size() == alignLen);
		alignLen = 0;

		// Now we don't care about left/right, so just call them "blocks"
		List<List<ResiduePair>> blocks = new ArrayList<List<ResiduePair>>(2);
		if( !left.isEmpty() ) {
			blocks.add(left);
			alignLen += left.size();
		}
		if( !right.isEmpty()) {
			blocks.add(right);
			alignLen += right.size();
		}
		left=null; right = null;

		// Put the blocks back together into arrays for the AFPChain
		int[][][] newAlign = new int[blocks.size()][][];
		int[] blockLengths = new int[blocks.size()];
		for(int blockNum = 0; blockNum < blocks.size(); blockNum++) {
			//Alignment
			List<ResiduePair> block = blocks.get(blockNum);
			newAlign[blockNum] = new int[2][block.size()];
			for(int i=0;i<block.size();i++) {
				ResiduePair pair = block.get(i);
				newAlign[blockNum][0][i] = pair.a;
				newAlign[blockNum][1][i] = pair.b;
			}

			// Block lengths
			blockLengths[blockNum] = block.size();
		}
		// Set Alignment
		newAFPChain.setOptAln(newAlign);
		newAFPChain.setOptLen(blockLengths );
		newAFPChain.setOptLength(alignLen);
		newAFPChain.setBlockNum(blocks.size());
		newAFPChain.setBlockResSize(blockLengths.clone());
		newAFPChain.setSequentialAlignment(blocks.size() == 1);

		// TODO make the AFPSet consistent
		// TODO lots more block properties & old AFP properties

		// Recalculate superposition
		Atom[] atoms1 = new Atom[alignLen];
		Atom[] atoms2 = new Atom[alignLen];

		int pos=0;
		for(List<ResiduePair> block:blocks ) {
			for(ResiduePair pair:block) {
				atoms1[pos] = ca1[pair.a];
				// Clone residue to allow modification
				Atom atom2 = ca2duplicated[pair.b];
				Group g = (Group) atom2.getGroup().clone();
				atoms2[pos] = g.getAtom( atom2.getName() );
				pos++;
			}
		}
		assert(pos == alignLen);

		// Sets the rotation matrix in ceCalc to the proper value
		double rmsd = -1;
		double tmScore = 0.;
		double[] blockRMSDs = new double[blocks.size()];
		Matrix[] blockRotationMatrices = new Matrix[blocks.size()];
		Atom[] blockShifts = new Atom[blocks.size()];

		if(alignLen>0) {
			
			// superimpose
			Matrix4d trans = SuperPositions.superpose(Calc.atomsToPoints(atoms1), 
					Calc.atomsToPoints(atoms2));

			Matrix matrix = Matrices.getRotationJAMA(trans);
			Atom shift = Calc.getTranslationVector(trans);

			for( Atom a : atoms2 ) 
				Calc.transform(a.getGroup(), trans);

			//and get overall rmsd
			rmsd = Calc.rmsd(atoms1, atoms2);
			tmScore = Calc.getTMScore(atoms1, atoms2, ca1.length, ca2len);

			// set all block rotations to the overall rotation
			// It's not well documented if this is the expected behavior, but
			// it seems to work.
			blockRotationMatrices[0] = matrix;
			blockShifts[0] = shift;
			blockRMSDs[0] = -1;

			for(int i=1;i<blocks.size();i++) {
				blockRMSDs[i] = -1; //TODO Recalculate for the FATCAT text format
				blockRotationMatrices[i] = (Matrix) blockRotationMatrices[0].clone();
				blockShifts[i] = (Atom) blockShifts[0].clone();
			}

		}
		newAFPChain.setOptRmsd(blockRMSDs);
		newAFPChain.setBlockRmsd(blockRMSDs);
		newAFPChain.setBlockRotationMatrix(blockRotationMatrices);
		newAFPChain.setBlockShiftVector(blockShifts);
		newAFPChain.setTotalRmsdOpt(rmsd);
		newAFPChain.setTMScore( tmScore );

		// Clean up remaining properties using the FatCat helper method
		Atom[] ca2 = new Atom[ca2len];
		System.arraycopy(ca2duplicated, 0, ca2, 0, ca2len);
		AFPAlignmentDisplay.getAlign(newAFPChain, ca1, ca2duplicated);

		return newAFPChain;
	}

	private static int[] countCtermResidues(int[] block, int blockLen,
			int cEnd, int overlapLength) {
		int[] cTermResCount = new int[overlapLength+1]; // # res at or to the right of i within overlap
		cTermResCount[overlapLength] = 0;
		int alignPos = blockLen - 1;
		for(int i=overlapLength-1;i>=0;i--) { // i starts at the c-term and increases to the left
			if(block[alignPos] == cEnd - overlapLength+1 + i) { // matches the aligned pair
				// the c-term contains the -ith overlapping residue
				cTermResCount[i] = cTermResCount[i+1]+1;
				alignPos--;
			} else {
				cTermResCount[i] = cTermResCount[i+1];
			}
		}
		return cTermResCount;
	}

	private static int[] countNtermResidues(int[] block, int nStart,
			int overlapLength) {
		int[] nTermResCount = new int[overlapLength+1]; // increases monotonically
		nTermResCount[0] = 0;
		int alignPos = 0; // index of the next aligned pair

		for(int i=1;i<=overlapLength;i++) {
			if(block[alignPos] == nStart + i-1 ) { // matches the aligned pair
				// the n-term contains the ith overlapping residue
				nTermResCount[i] = nTermResCount[i-1]+1;
				alignPos++;
			} else {
				nTermResCount[i] = nTermResCount[i-1];
			}
		}
		return nTermResCount;
	}


	/**
	 * A light class to store an alignment between two residues.
	 * @author Spencer Bliven
	 * @see #filterDuplicateAFPs()
	 */
	private static class ResiduePair {
		public int a;
		public int b;
		public ResiduePair(int a, int b) {
			this.a=a;
			this.b=b;
		}
		@Override
		public String toString() {
			return a+":"+b;
		}
	}


	/**
	 * Tiny wrapper for the disallowed regions of an alignment.
	 * @see CeCPMain#calculateMinCP(int[], int, int, int)
	 * @author Spencer Bliven
	 *
	 */
	protected static class CPRange {
		/**
		 * last allowed n-term
		 */
		public int n;
		/**
		 * midpoint of the alignment
		 */
		public int mid;
		/**
		 * first allowed c-term
		 */
		public int c;
	}

	/**
	 * Finds the alignment index of the residues minCPlength before and after
	 * the duplication.
	 *
	 * @param block The permuted block being considered, generally optAln[0][1]
	 * @param blockLen The length of the block (in case extra memory was allocated in block)
	 * @param ca2len The length, in residues, of the protein specified by block
	 * @param minCPlength The minimum number of residues allowed for a CP
	 * @return a CPRange with the following components:
	 *  <dl><dt>n</dt><dd>Index into <code>block</code> of the residue such that
	 *  	<code>minCPlength</code> residues remain to the end of <code>ca2len</code>,
	 *  	or -1 if no residue fits that criterium.</dd>
	 *  <dt>mid</dt><dd>Index of the first residue higher than <code>ca2len</code>.</dd>
	 *  <dt>c</dt><dd>Index of <code>minCPlength</code>-th residue after ca2len,
	 *  	or ca2len*2 if no residue fits that criterium.</dd>
	 *  </dl>
	 */
	protected static CPRange calculateMinCP(int[] block, int blockLen, int ca2len, int minCPlength) {
		CPRange range = new CPRange();

		// Find the cut point within the alignment.
		// Either returns the index i of the alignment such that block[i] == ca2len,
		// or else returns -i-1 where block[i] is the first element > ca2len.
		int middle = Arrays.binarySearch(block, ca2len);
		if(middle < 0) {
			middle = -middle -1;
		}
		// Middle is now the first res after the duplication
		range.mid = middle;

		int minCPntermIndex = middle-minCPlength;
		if(minCPntermIndex >= 0) {
			range.n = block[minCPntermIndex];
		} else {
			range.n = -1;
		}

		int minCPctermIndex = middle+minCPlength-1;
		if(minCPctermIndex < blockLen) {
			range.c = block[minCPctermIndex];
		} else {
			range.c = ca2len*2;
		}

		// Stub:
		// Best-case: assume all residues in the termini are aligned
		//range.n = ca2len - minCPlength;
		//range.c = ca2len + minCPlength-1;

		return range;
	}


	private static class CutPoint {
		public int numResiduesCut;
		public int firstRes;
		public int lastRes;
	}

	private static CutPoint calculateCutPoint(int[] block,int nStart, int cEnd,
			int overlapLength, int alignLen, int minCPlength, int ca2len, int firstRes) {

		// We require at least minCPlength residues in a block.
		//TODO calculate these explicitely based on the alignment

		// The last valid n-term
		CPRange minCP = calculateMinCP(block, alignLen, ca2len, minCPlength);
		int minCPnterm = minCP.n;
		// The first valid c-term
		int minCPcterm = minCP.c;


		// # res at or to the left of i within the overlap
		int[] nTermResCount = countNtermResidues(block, nStart,
				overlapLength);

		// Determine the position with the largest sum of lengths
		int[] cTermResCount = countCtermResidues(block, alignLen,
				cEnd, overlapLength);

		// Alignment length for a cut at the left of the overlap
		int maxResCount=-1;
		for(int i=0;i<=overlapLength;i++) { // i is the position of the CP within the overlap region
			// Calculate number of residues which remain after the CP
			int nRemain,cRemain;
			if(nStart+i <= minCPnterm) {
				nRemain = nTermResCount[overlapLength]-nTermResCount[i];
			} else {
				nRemain = 0;
			}
			if(cEnd-overlapLength+i >= minCPcterm) {
				cRemain = cTermResCount[0] - cTermResCount[i];
			} else {
				cRemain = 0;
			}

			// Look for the cut point which cuts off the minimum number of res
			if(nRemain + cRemain > maxResCount ) { // '>' biases towards keeping the n-term
				maxResCount = nRemain + cRemain;
				firstRes = nStart+ i;
			}
		}

		// Calculate the number of residues cut within the overlap
		int numResiduesCut = nTermResCount[overlapLength]+cTermResCount[0]-maxResCount;

		// Remove short blocks
		if(firstRes > minCPnterm) {
			// update number of residues cut for those outside the overlap
			numResiduesCut += 0; //TODO

			firstRes = ca2len;
		}

		int lastRes = firstRes+ca2len-1;
		if(lastRes < minCPcterm) {
			// update number of residues cut for those outside the overlap
			numResiduesCut += 0; //TODO

			lastRes = ca2len-1;
		}


		CutPoint cp = new CutPoint();
		cp.firstRes=firstRes;
		cp.numResiduesCut = numResiduesCut;
		cp.lastRes = lastRes;

		if(debug) {
			System.out.format("Found a CP at residue %d. Trimming %d aligned residues from %d-%d of block 0 and %d-%d of block 1.\n",
					firstRes,cp.numResiduesCut,nStart,firstRes-1,firstRes, cEnd-ca2len);
		}

		return cp;
	}


	// try showing a GUI
	// requires additional dependencies biojava-structure-gui and JmolApplet
	// TODO dmyersturnbull: This should probably be in structure-gui
	@SuppressWarnings("unused")
	private static void displayAlignment(AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException, StructureException {
		Atom[] ca1clone = StructureTools.cloneAtomArray(ca1);
		Atom[] ca2clone = StructureTools.cloneAtomArray(ca2);
		if (! GuiWrapper.isGuiModuleInstalled()) {
			System.err.println("The biojava-structure-gui and/or JmolApplet modules are not installed. Please install!");
			// display alignment in console
			System.out.println(afpChain.toCE(ca1clone, ca2clone));
		} else {
			Object jmol = GuiWrapper.display(afpChain,ca1clone,ca2clone);
			GuiWrapper.showAlignmentImage(afpChain, ca1clone,ca2clone,jmol);
		}
	}
}
