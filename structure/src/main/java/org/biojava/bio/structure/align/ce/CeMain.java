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
 * Created on Sep 15, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.ce;


import java.util.ArrayList;
import java.util.List;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.AbstractStructureAlignment;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.bio.structure.jama.Matrix;

/** The main class of the Java implementation of the Combinatorial Extension Algorithm (CE).
 * 
 * The original CE paper is available from here: <a href="http://peds.oxfordjournals.org/cgi/content/short/11/9/739">http://peds.oxfordjournals.org/cgi/content/short/11/9/739</a>
 * 
 * For a demo of how to use this algorithm, visit the BioJava web site:
 * <a href="">CE usage example</a>.
 * 
 * @author Andreas Prlic.
 *
 */
public class CeMain extends AbstractStructureAlignment implements StructureAlignment {

	public static final String algorithmName = "jCE";

	public static final String version = "1.0";

	protected CeParameters params;
	protected CECalculator calculator;
	protected Atom[] ca2clone;

	public CeMain(){
		super();
		params = new CeParameters();
		calculator = new CECalculator(params);
	}


	public static void main(String[] args){

		CeMain ce = new CeMain(); //Used only for printing help
		if (args.length  == 0 ) {			
			System.out.println(ce.printHelp());
			return;			
		}

		if ( args.length == 1){
			if (args[0].equalsIgnoreCase("-h") || args[0].equalsIgnoreCase("-help")|| args[0].equalsIgnoreCase("--help")){
				System.out.println(ce.printHelp());								
				return;
			}

		}

		CeUserArgumentProcessor processor = new CeUserArgumentProcessor(); //Responsible for creating a CeMain instance
		processor.process(args);

	}

	/**
	 * Align ca2 onto ca1.
	 */
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object param) throws StructureException{

		if ( ! (param instanceof CeParameters))
			throw new IllegalArgumentException("CE algorithm needs an object of call CeParameters as argument.");

		params = (CeParameters) param;

		int ca2length = params.getCheckCircular()? ca2.length*2 : ca2.length; 


		// we don't want to rotate input atoms, do we?
		ca2clone = new Atom[ca2length];

		int pos = 0;
		for (Atom a : ca2){
			Group g = (Group)a.getParent().clone(); // works because each group has only a CA atom

			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos++;
		}

		if(params.getCheckCircular()) {
			//System.out.println("Checking Circular permutations");

			// Duplicate ca2
			for (Atom a : ca2){
				Group g = (Group)a.getParent().clone();

				ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

				pos++;
			}
		}

		calculator = new CECalculator(params);

		//Build alignment ca1 to ca2-ca2
		AFPChain afpChain = new AFPChain();
		afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);
		calculator.traceFragmentMatrix( afpChain,ca1, ca2clone);
		calculator.nextStep( afpChain,ca1, ca2clone);

		afpChain.setAlgorithmName(getAlgorithmName());
		afpChain.setVersion(version);
		
		
		if ( afpChain.getNrEQR() == 0)
		   return afpChain;

		// Set the distance matrix
		int winSize = params.getWinSize();
		int winSizeComb1 = (winSize-1)*(winSize-2)/2;	
		double[][] m = calculator.initSumOfDistances(ca1.length, ca2length, winSize, winSizeComb1, ca1, ca2clone);

		// Draw a little green line across the center of the distance matrix
		/* TODO do this more elegantly. If anyone uses afpChain.getDistanceMatrix()
		 * for anything besides displaying a dotplot, they will be confused by
		 * the row of -1 through the middle
		 */
		if(params.getCheckCircular()) {

			//System.out.println("Running in Circular Permutation mode.");
			for(int i=0;i<m.length;i++) {
				m[i][ca2.length] = -1;
			}
		}

		afpChain.setDistanceMatrix(new Matrix(m));


		if(params.getCheckCircular()) {
			// Check for circular permutations
			afpChain = filterDuplicateAFPs(afpChain,calculator,ca1,ca2clone);
		}


		return afpChain;
	}



	/**
	 * Takes as input an AFPChain where ca2 has been artificially duplicated.
	 * This raises the possibility that some residues of ca2 will appear in 
	 * multiple AFPs. This method filters out duplicates and makes sure that
	 * all AFPs are numbered relative to the original ca2.
	 * 
	 * The current version chooses a CP site such that the length of the
	 * alignment is maximized.
	 * 
	 * @param afpChain The alignment between ca1 and ca2-ca2. Blindly assumes 
	 *  that ca2 has been duplicated.
	 * @return A new AFPChain consisting of ca1 to ca2, with each residue in
	 *  at most 1 AFP.
	 * @throws StructureException 
	 */
	private static AFPChain filterDuplicateAFPs(AFPChain afpChain, CECalculator ceCalc, Atom[] ca1, Atom[] ca2duplicated) throws StructureException {
		AFPChain newAFPChain = new AFPChain(afpChain);
		newAFPChain.setAlgorithmName(afpChain.getAlgorithmName());
		newAFPChain.setVersion(afpChain.getVersion());
		newAFPChain.setName1(afpChain.getName1());
		newAFPChain.setName2(afpChain.getName2());

		int ca2len = afpChain.getCa2Length()/2;

		// Fix optimal alignment		
		int[][][] align = afpChain.getOptAln();
		int alignLen = afpChain.getOptLength();
		assert(align.length == 1); // Assume that CE returns just one block

		// Determine the region where ca2 and ca2' overlap
		int nStart = align[0][1][0]; //alignment N-terminal
		int cEnd = align[0][1][alignLen-1]; // alignment C-terminal 
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
				firstRes=nStart;
				lastRes=cEnd;
			}
			else {
				// overlap!
				
				// # res at or to the left of i within the overlap
				int[] nTermResCount = new int[overlapLength]; // increases monotonically
				nTermResCount[0] = 1;

				int alignPos = 1; // index of the next aligned pair
				
				for(int i=1;i<overlapLength;i++) {
					if(align[0][1][alignPos] == nStart + i ) { // matches the aligned pair
						// the n-term contains the ith overlapping residue
						nTermResCount[i] = nTermResCount[i-1]+1;
						alignPos++;
					}
					else {
						nTermResCount[i] = nTermResCount[i-1];
					}
				}

				// Determine the position with the largest sum of lengths
				int cTermResCount = 0;
				alignPos = alignLen - 1;
				int minResCount = overlapLength+1;
				for(int i=0;i<overlapLength;i++) {

					if(nTermResCount[overlapLength-1-i] + cTermResCount <= minResCount ) {
						minResCount=nTermResCount[overlapLength-1-i] + cTermResCount;
						firstRes = nStart+overlapLength - i;
						lastRes = cEnd - i;
						assert(lastRes == firstRes+ca2len-1);
					}

					if(align[0][1][alignPos] == cEnd - i) {
						cTermResCount++;
						alignPos--;
					}
				}
			
				//Adjust alignment length for trimming
				alignLen -= minResCount;
				

				System.out.format("Found a CP at residue %d. Trimming %d residues (%d-%d,%d-%d).\n",
						firstRes,minResCount,nStart,firstRes-1,firstRes+ca2len, cEnd);
				//TODO Now have CP site, and could do a nxm alignment for further optimization.
				// For now, does not appear to be worth the 50% increase in time
			}
		}



		// Fix numbering:
		// First, split up the atoms into left and right blocks
		List< ResiduePair > left = new ArrayList<ResiduePair>(); // residues from left of duplication
		List< ResiduePair > right = new ArrayList<ResiduePair>(); // residues from right of duplication

		for(int i=0;i<afpChain.getOptLength();i++) {
			if( align[0][1][i] >= firstRes && align[0][1][i] <= lastRes ) { // not trimmed
				if(align[0][1][i] < ca2len) { // in first half of ca2
					left.add(new ResiduePair(align[0][0][i],align[0][1][i]));
				}
				else {
					right.add(new ResiduePair(align[0][0][i],align[0][1][i]-ca2len));
				}
			}
		}
		assert(left.size()+right.size() == alignLen);


		// Now we don't care about left/right, so just call them "blocks"
		List<List<ResiduePair>> blocks = new ArrayList<List<ResiduePair>>(2);
		if( !left.isEmpty() ) {
			blocks.add(left);
		}
		if( !right.isEmpty()) {
			blocks.add(right);
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
				atoms2[pos] = ca2duplicated[pair.b];
				pos++;
			}
		}
		assert(pos == alignLen);

		// Sets the rotation matrix in ceCalc to the proper value
		double rmsd = ceCalc.calc_rmsd(atoms1, atoms2, alignLen, true, false);

		double[] blockRMSDs = new double[blocks.size()];
		Matrix[] blockRotationMatrices = new Matrix[blocks.size()];
		Atom[] blockShifts = new Atom[blocks.size()];

		blockRMSDs[0] = rmsd;
		blockRotationMatrices[0] = ceCalc.getRotationMatrix();
		blockShifts[0] = ceCalc.getShift();
		for(int i=1;i<blocks.size();i++) {
			blockRMSDs[i] = rmsd;

			// Don't move blocks relative to the first block
			/*Matrix identity = new Matrix(3,3);
			for(int j=0;j<3;j++)
				identity.set(j, j, 1.);
			blockRotationMatrices[i] = identity;

			Atom zero = new AtomImpl();
			zero.setX(0.); zero.setY(0.); zero.setZ(0.);
			blockShifts[i] = zero;
			 */
			blockRotationMatrices[i] = (Matrix) blockRotationMatrices[0].clone();
			blockShifts[i] = (Atom) blockShifts[0].clone();
		}
		newAFPChain.setOptRmsd(blockRMSDs);
		newAFPChain.setBlockRmsd(blockRMSDs);
		newAFPChain.setBlockRotationMatrix(blockRotationMatrices);
		newAFPChain.setBlockShiftVector(blockShifts);

		// Clean up remaining properties using the FatCat helper method
		Atom[] ca2 = new Atom[ca2len];
		for(int i=0;i<ca2len;i++) {
			ca2[i]=ca2duplicated[i];
		}
		AFPAlignmentDisplay.getAlign(newAFPChain, ca1, ca2duplicated);
//		return afpChain;

		return newAFPChain;
	}


	public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException {

		if (params == null)
			params = new CeParameters();

		return align(ca1,ca2,params);
	}

	public String getAlgorithmName() {

		return CeMain.algorithmName;
	}

	public ConfigStrucAligParams getParameters() {

		return params;
	}

	public void setParameters(ConfigStrucAligParams params){
		if (! (params instanceof CeParameters )){
			throw new IllegalArgumentException("provided parameter object is not of type CeParameter");
		}
		this.params = (CeParameters) params;
	}

	public String getVersion() {
		return version;
	}

	public CECalculator getCECalculator() {
		return calculator;
	}

	/**
	 * A light class to store an alignment between two residues
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
	}
}
