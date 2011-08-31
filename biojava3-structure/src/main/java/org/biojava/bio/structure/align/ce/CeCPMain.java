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

package org.biojava.bio.structure.align.ce;


import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.jama.Matrix;

/** 
 * A wrapper for {@link CeMain} which sets default parameters to be appropriate for finding
 * circular permutations
 * 
 * @author Spencer Bliven.
 *
 */
public class CeCPMain extends CeMain {
	private static boolean debug = false;

	public static final String algorithmName = "jCE Circular Permutation";

	public static final String version = "1.0";

	public CeCPMain(){
		super();
		this.params.setMaxGapSize(0);
	}

	public String getAlgorithmName() {

		return algorithmName;
	}

	public String getVersion() {
		return version;
	}


	public static void main(String[] args){
		try {
			String name1, name2;

			//Concanavalin
			name1 = "2pel.A";
			name2 = "3cna";

			//small case
			//name1 = "1QDM.A";
			//name2 = "1NKL";

			CeMain ce = (CeMain) StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
			CeParameters params = (CeParameters) ce.getParameters();
			ce.setParameters(params);

			AtomCache cache = new AtomCache();

			Atom[] ca1 = cache.getAtoms(name1);
			Atom[] ca2 = cache.getAtoms(name2);
			ca1[0].getGroup().getChain().getParent().setName(name1);
			ca2[0].getGroup().getChain().getParent().setName(name2);


			System.out.format("Aligning %s to %s\n",
					ca1[1].getGroup().getChain().getParent().getName(),
					ca2[1].getGroup().getChain().getParent().getName());
			AFPChain afpChain = ce.align(ca1, ca2);


			// try showing a GUI
			// requires additional dependancies biojava3-structure-gui and JmolApplet
			if (! GuiWrapper.isGuiModuleInstalled()) {
				System.err.println("The biojava-structure-gui and/or JmolApplet modules are not installed. Please install!");
				// display alignment in console
				System.out.println(afpChain.toCE(ca1, ca2));
			} else {
				Object jmol = GuiWrapper.display(afpChain,ca1,ca2);
				GuiWrapper.showAlignmentImage(afpChain, ca1,ca2,jmol);
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object param) throws StructureException{
		long startTime = System.currentTimeMillis();

		// Duplicate ca2
		Atom[] ca2m = new Atom[ca2.length*2];

		int pos = 0;
		for (Atom a : ca2){
			Group g = (Group)a.getGroup().clone(); // works because each group has only a CA atom
			ca2m[pos] = g.getAtom(StructureTools.caAtomName);
			pos++;
		}
		for (Atom a : ca2){
			Group g = (Group)a.getGroup().clone();
			ca2m[pos] = g.getAtom(StructureTools.caAtomName);
			pos++;
		}

		if(debug) {
			System.out.format("Duplicating ca2 took %s ms\n",System.currentTimeMillis()-startTime);
			startTime = System.currentTimeMillis();
		}

		// Do alignment
		AFPChain afpChain = super.align(ca1, ca2m,params);

		if(debug) {
			System.out.format("Running Nx2N alignment took %s ms\n",System.currentTimeMillis()-startTime);
			startTime = System.currentTimeMillis();
		}

		// Draw a little green line across the center of the distance matrix
		/* TODO do this more elegantly. If anyone uses afpChain.getDistanceMatrix()
		 * for anything besides displaying a dotplot, they will be confused by
		 * the row of -1 through the middle
		 */

		Matrix doubledMatrix = afpChain.getDistanceMatrix();

		assert(doubledMatrix.getRowDimension() == ca1.length);
		assert(doubledMatrix.getColumnDimension() == ca2.length*2);
		/*double[][] doubledArray = doubledMatrix.getArray();

		double[][] singleArray = new double[ca1.length][ca2.length];
		// copy the first half of the matrix, discarding the rest
		for(int i=0;i<ca1.length;i++) {
			for(int j=0;j<ca2.length;j++) {
				singleArray[i][j] = doubledArray[i][j];
			}
		afpChain.setDistanceMatrix(new Matrix(singleArray));
		}*/
		Matrix singleMatrix = doubledMatrix.getMatrix(0, ca1.length-1, 0, ca2.length-1);
		assert(singleMatrix.getRowDimension() == ca1.length);
		assert(singleMatrix.getColumnDimension() == ca2.length);

		afpChain.setDistanceMatrix(singleMatrix);

		// Check for circular permutations
		afpChain = filterDuplicateAFPs(afpChain,calculator,ca1,ca2m);

		if(debug) {
			System.out.format("Finding CP point took %s ms\n",System.currentTimeMillis()-startTime);
			startTime = System.currentTimeMillis();
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
	public static AFPChain filterDuplicateAFPs(AFPChain afpChain, CECalculator ceCalc, Atom[] ca1, Atom[] ca2duplicated) throws StructureException {
		AFPChain newAFPChain = new AFPChain(afpChain);
		newAFPChain.setAlgorithmName(afpChain.getAlgorithmName());
		newAFPChain.setVersion(afpChain.getVersion());
		newAFPChain.setName1(afpChain.getName1());
		newAFPChain.setName2(afpChain.getName2());
		newAFPChain.setTMScore(afpChain.getTMScore());
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


				//System.out.format("Found a CP at residue %d. Trimming %d residues (%d-%d,%d-%d).\n",
				//		firstRes,minResCount,nStart,firstRes-1,firstRes+ca2len, cEnd);
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
