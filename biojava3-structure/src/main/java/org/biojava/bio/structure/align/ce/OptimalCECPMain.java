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


import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.jama.Matrix;

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
public class OptimalCECPMain extends CeMain {
	private static boolean debug = true;


	public static final String algorithmName = "jCE Optimal Circular Permutation";

	public static final String version = "1.0";

	protected OptimalCECPParameters params;
	
	/**
	 * 
	 */
	public OptimalCECPMain() {
		super();
		params = new OptimalCECPParameters();
	}
	
	@Override
	public String getAlgorithmName() {
		return OptimalCECPMain.algorithmName;
	}

	@Override
	public String getVersion() {
		return OptimalCECPMain.version;
	}
	
	/**
	 * @return an {@link OptimalCECPParameters} object 
	 */
	@Override
	public ConfigStrucAligParams getParameters() {
		return params;
	}

	/**
	 * @param params Should be an {@link OptimalCECPParameters} object specifying alignment options
	 */
	@Override
	public void setParameters(ConfigStrucAligParams params){
		if (! (params instanceof OptimalCECPParameters )){
			throw new IllegalArgumentException("provided parameter object is not of type CeParameter");
		}
		this.params = (OptimalCECPParameters) params;
	}
	
	/**
	 * Circularly permutes arr in place.
	 * 
	 * <p>Similar to {@link Collections#rotate(List, int)} but with reversed
	 * direction. Perhaps it would be more efficient to use the Collections version?
	 * @param <T>
	 * @param arr The array to be permuted
	 * @param cp The number of residues to shift leftward, or equivalently, the index of
	 *  the first element after the permutation point.
	 */
	private static <T> void permuteArray(T[] arr, int cp) {
		// Allow negative cp points for convenience.
		if(cp == 0) {
			return;
		}
		if(cp < 0) {
			cp = arr.length+cp;
		}
		if(cp < 0 || cp >= arr.length) {
			throw new ArrayIndexOutOfBoundsException(
					"Permutation point ("+cp+") must be between -ca2.length and ca2.length-1" );
		}
		
		List<T> temp = new ArrayList<T>(cp);

		// shift residues left
		for(int i=0;i<cp;i++) {
			temp.add(arr[i]);
		}
		for(int j=cp;j<arr.length;j++) {
			arr[j-cp]=arr[j];
		}
		for(int i=0;i<cp;i++) {
			arr[arr.length-cp+i] = temp.get(i);
		}
	}
	
	/**
	 * Circularly permutes arr in place.
	 * 
	 * <p>Similar to {@link Collections#rotate(List, int)} but with reversed
	 * direction. Perhaps it would be more efficient to use the Collections version?
	 * @param <T>
	 * @param arr The array to be permuted
	 * @param cp The number of residues to shift leftward, or equivalently, the index of
	 *  the first element after the permutation point.
	 *
	private static void permuteArray(int[] arr, int cp) {
		// Allow negative cp points for convenience.
		if(cp == 0) {
			return;
		}
		if(cp < 0) {
			cp = arr.length+cp;
		}
		if(cp < 0 || cp >= arr.length) {
			throw new ArrayIndexOutOfBoundsException(
					"Permutation point ("+cp+") must be between -ca2.length and ca2.length-1" );
		}
		
		List<Integer> temp = new ArrayList<Integer>(cp);

		// shift residues left
		for(int i=0;i<cp;i++) {
			temp.add(arr[i]);
		}
		for(int j=cp;j<arr.length;j++) {
			arr[j-cp]=arr[j];
		}
		for(int i=0;i<cp;i++) {
			arr[arr.length-cp+i] = temp.get(i);
		}
	}
	*/

	/**
	 * Aligns ca1 with ca2 permuted by <i>cp</i> residues.
	 * <p><strong>WARNING:</strong> Modifies ca2 during the permutation. Be sure
	 * to make a copy before calling this method.
	 * 
	 * @param ca1
	 * @param ca2
	 * @param param
	 * @param cp
	 * @return
	 * @throws StructureException 
	 */
	public AFPChain alignPermuted(Atom[] ca1, Atom[] ca2, Object param, int cp) throws StructureException {
		// initial permutation
		permuteArray(ca2,cp);
				
		// perform alignment
		AFPChain afpChain = super.align(ca1, ca2, param);
		
		// un-permute alignment
		permuteAFPChain(afpChain, -cp);
		
		if(afpChain.getName2() != null) {
			afpChain.setName2(afpChain.getName2()+" CP="+cp);
		}
		
		// Specify the permuted
		return afpChain;
	}

	/**
	 * Permute the second protein of afpChain by the specified number of residues.
	 * @param afpChain Input alignment
	 * @param cp Amount leftwards (or rightward, if negative) to shift the 
	 * @return A new alignment equivalent to afpChain after the permutations
	 */
	private static void permuteAFPChain(AFPChain afpChain, int cp) {
		
		int ca2len = afpChain.getCa2Length();
				
		//fix up cp to be positive
		if(cp == 0) {
			return;
		}
		if(cp < 0) {
			cp = ca2len+cp;
		}
		if(cp < 0 || cp >= ca2len) {
			throw new ArrayIndexOutOfBoundsException(
					"Permutation point ("+cp+") must be between -ca2.length and ca2.length-1" );
		}
		
		// Fix up optAln
		permuteOptAln(afpChain,cp);
		
		if(afpChain.getBlockNum() > 1)
			afpChain.setSequentialAlignment(false);
		// fix up matrices
		// ca1 corresponds to row indices, while ca2 corresponds to column indices.
		afpChain.setDistanceMatrix(permuteMatrix(afpChain.getDistanceMatrix(),0,-cp));
		// this is square, so permute both
		afpChain.setDisTable2(permuteMatrix(afpChain.getDisTable2(),-cp,-cp));

		//TODO fix up other AFP parameters?
		
	}
	
	/**
	 * Permutes <i>mat</i> by moving the rows of the matrix upwards by <i>cp</i>
	 * rows.
	 * @param mat The original matrix
	 * @param cpRows Number of rows upward to move entries
	 * @param cpCols Number of columns leftward to move entries
	 * @return The permuted matrix
	 */
	private static Matrix permuteMatrix(Matrix mat, int cpRows, int cpCols) {
		//fix up cp to be positive
		if(cpRows == 0 && cpCols == 0) {
			return mat.copy();
		}
		if(cpRows < 0) {
			cpRows = mat.getRowDimension()+cpRows;
		}
		if(cpRows < 0 || cpRows >= mat.getRowDimension()) {
			throw new ArrayIndexOutOfBoundsException( String.format(
					"Can't permute rows by %d: only %d rows.", 
					cpRows, mat.getRowDimension() )
			);
		}

		if(cpCols < 0) {
			cpCols = mat.getColumnDimension()+cpCols;
		}
		if(cpCols < 0 || cpCols >= mat.getColumnDimension()) {
			throw new ArrayIndexOutOfBoundsException( String.format(
					"Can't permute cols by %d: only %d rows.", 
					cpCols, mat.getColumnDimension() )
			);
		}
		
		int[] rows = new int[mat.getRowDimension()];
		for(int i=0;i<rows.length;i++) {
			rows[i] = (i+cpRows)%rows.length;
		}
		int[] cols = new int[mat.getColumnDimension()];
		for(int i=0;i<cols.length;i++) {
			cols[i] = (i+cpCols)%cols.length;
		}
		
		Matrix newMat = mat.getMatrix(rows, cols);
		assert(newMat.getRowDimension() == mat.getRowDimension());
		assert(newMat.getColumnDimension() == mat.getColumnDimension());
		assert(newMat.get(0, 0) ==
			mat.get(cpRows%mat.getRowDimension(), cpCols%mat.getColumnDimension()));
		
		
		return newMat;
	}

	/**
	 * Modifies the {@link AFPChain#setOptAln(int[][][]) optAln} of an AFPChain
	 * by permuting the second protein.
	 * 
	 * Sets residue numbers in the second protein to <i>(i-cp)%len</i>
	 * 
	 * @param afpChain
	 * @param cp Amount leftwards (or rightward, if negative) to shift the
	 */
	private static void permuteOptAln(AFPChain afpChain, int cp)
	{
		int ca2len = afpChain.getCa2Length();

		if( ca2len <= 0) {
			throw new IllegalArgumentException("No Ca2Length specified in "+afpChain);
		}

		// Allow negative cp points for convenience.
		if(cp == 0) {
			return;
		}
		if(cp <= -ca2len || cp >= ca2len) {
			// could just take cp%ca2len, but probably its a bug if abs(cp)>=ca2len
			throw new ArrayIndexOutOfBoundsException( String.format(
					"Permutation point %d must be between %d and %d for %s",
					cp, 1-ca2len,ca2len-1, afpChain.getName2() ) );
		}
		if(cp < 0) {
			cp = cp + ca2len;
		}

		// the unprocessed alignment
		int[][][] optAln = afpChain.getOptAln();
		int[] optLen = afpChain.getOptLen();

		// the processed alignment
		List<List<List<Integer>>> blocks = new ArrayList<List<List<Integer>>>(afpChain.getBlockNum()*2);
		
		//Update residue indices
		// newi = (oldi-cp) % N
		for(int block = 0; block < afpChain.getBlockNum(); block++) {
			if(optLen[block]<1)
				continue;
			
			// set up storage for the current block
			List<List<Integer>> currBlock = new ArrayList<List<Integer>>(2);
			currBlock.add( new ArrayList<Integer>());
			currBlock.add( new ArrayList<Integer>());
			blocks.add(currBlock);

			// pos = 0 case
			currBlock.get(0).add( optAln[block][0][0] );
			currBlock.get(1).add( (optAln[block][1][0]+cp ) % ca2len);
			
			for(int pos = 1; pos < optLen[block]; pos++) {
				//check if we need to start a new block
				//this happens when the new alignment crosses the protein terminus
				if( optAln[block][1][pos-1]+cp<ca2len && 
						optAln[block][1][pos]+cp >= ca2len) {
					currBlock = new ArrayList<List<Integer>>(2);
					currBlock.add( new ArrayList<Integer>());
					currBlock.add( new ArrayList<Integer>());
					blocks.add(currBlock);
				}
				currBlock.get(0).add( optAln[block][0][pos] );
				currBlock.get(1).add( (optAln[block][1][pos]+cp ) % ca2len);
			}
		}
		
		// save permuted blocks to afpChain
		assignOptAln(afpChain,blocks);
	}

	/**
	 * Sometimes it's convenient to store an alignment using java collections,
	 * where <tt>blocks.get(blockNum).get(0).get(pos)</tt> specifies the aligned
	 * residue at position <i>pos</i> of block <i>blockNum</i> of the first
	 * protein.
	 * 
	 * This method takes such a collection and stores it into the afpChain's
	 * {@link AFPChain#setOptAln(int[][][]) optAln}, setting the associated
	 * length variables as well.
	 * 
	 * @param afpChain
	 * @param blocks
	 */
	private static void assignOptAln(AFPChain afpChain, List<List<List<Integer>>> blocks)
	{
		
		int[][][] optAln = new int[blocks.size()][][];
		int[] optLen = new int[blocks.size()];
		int optLength = 0;
		int numBlocks = blocks.size();
		
		for(int block = 0; block < numBlocks; block++) {
			// block should be 2xN rectangular
			assert(blocks.get(block).size() == 2);
			assert( blocks.get(block).get(0).size() == blocks.get(block).get(1).size());

			optLen[block] = blocks.get(block).get(0).size();
			optLength+=optLen[block];
			
			optAln[block] = new int[][] {
					new int[optLen[block]],
					new int[optLen[block]]
			};
			for(int pos = 0; pos < optLen[block]; pos++) {
				optAln[block][0][pos] = blocks.get(block).get(0).get(pos);
				optAln[block][1][pos] = blocks.get(block).get(1).get(pos);
			}
		}
		
		
		afpChain.setBlockNum(numBlocks);
		afpChain.setOptAln(optAln);
		afpChain.setOptLen(optLen);
		afpChain.setOptLength(optLength);
		
		// TODO I don't know what these do. Should they be set?
		//afpChain.setBlockSize(blockSize);
		//afpChain.setBlockResList(blockResList);
		//afpChain.setChainLen(chainLen);
		
	}

	/**
	 * Finds the optimal alignment between two proteins allowing for a circular
	 * permutation (CP).
	 * 
	 * The precise algorithm is controlled by the 
	 * {@link OptimalCECPParameters parameters}. If the parameter
	 * {@link OptimalCECPParameters#isTryAllCPs() tryAllCPs} is true, all possible
	 * CP sites are tried and the optimal site is returned. Otherwise, the
	 * {@link OptimalCECPParameters#getCPPoint() cpPoint} parameter is used to
	 * determine the CP point, greatly reducing the computation required.
	 * 
	 * @param ca1 CA atoms of the first protein
	 * @param ca2 CA atoms of the second protein
	 * @param param {@link CeParameters} object
	 * @return The best-scoring alignment
	 * @throws StructureException
	 * 
	 * @see #alignOptimal(Atom[], Atom[], Object, AFPChain[])
	 */
	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object param)
	throws StructureException
	{
		if(params.isTryAllCPs()) {
			return alignOptimal(ca1,ca2,param,null);
		} else {
			int cpPoint = params.getCPPoint();
			return alignPermuted(ca1, ca2, param, cpPoint);
		}
	}

	/**
	 * Finds the optimal alignment between two proteins allowing for a circular
	 * permutation (CP).
	 * 
	 * This algorithm performs a CE alignment for each possible CP site. This is
	 * quite slow. Use {@link #alignHeuristic(Atom[], Atom[], Object)} for a
	 * faster algorithm.
	 * 
	 * @param ca1 CA atoms of the first protein
	 * @param ca2 CA atoms of the second protein
	 * @param param {@link CeParameters} object
	 * @param alignments If not null, should be an empty array of the same length as
	 *  ca2. This will be filled with the alignments from permuting ca2 by 
	 *  0 to n-1 residues.
	 * @return The best-scoring alignment
	 * @throws StructureException
	 */
	public AFPChain alignOptimal(Atom[] ca1, Atom[] ca2, Object param, AFPChain[] alignments)
	throws StructureException
	{
		long startTime = System.currentTimeMillis();
		
		if(alignments.length != ca2.length) {
			throw new IllegalArgumentException("scores param should have same length as ca2");
		}
		
		AFPChain unaligned = super.align(ca1, ca2, param);
		AFPChain bestAlignment = unaligned;
		
		if(debug) {
			// print progress bar header
			System.out.print("|");
			for(int cp=1;cp<ca2.length-1;cp++) {
				System.out.print("=");
			}
			System.out.println("|");
			System.out.print(".");
		}

		if(alignments != null) {
			alignments[0] = unaligned;
		}
		
		for(int cp=1;cp<ca2.length;cp++) {
			// clone ca2 to prevent side effects from propegating
			Atom[] ca2p = StructureTools.cloneCAArray(ca2);

			//permute one each time. Alters ca2p as a side effect
			AFPChain currentAlignment = alignPermuted(ca1,ca2p,param,cp);
			
			// increment progress bar
			if(debug) System.out.print(".");
			
			// fix up names, since cloning ca2 wipes it
			try {
				currentAlignment.setName2(ca2[0].getGroup().getChain().getParent().getName()+" CP="+cp);
			} catch( Exception e) {
				//null pointers, empty arrays, etc.
			}
			double currentScore = currentAlignment.getAlignScore();
			
			if(alignments != null) {
				alignments[cp] = currentAlignment;
			}
			
			if(currentScore>bestAlignment.getAlignScore()) {
				bestAlignment = currentAlignment;
			}
		}
		if(debug) {
			long elapsedTime = System.currentTimeMillis()-startTime;
			System.out.println();
			System.out.format("%d alignments took %.4f s (%.1f ms avg)\n",
					ca2.length, elapsedTime/1000., (double)elapsedTime/ca2.length);
		}

		
		return bestAlignment;
		
	}



	
	public static void main(String[] args){
		try {
			String name1, name2;

			int[] cps= new int[] {};
			
			//Concanavalin
			name1 = "2pel.A";
			name2 = "3cna";
			cps = new int[] {122,0,3};

			//small case
			//name1 = "d1qdmA1";
			//name1 = "1QDM.A";
			//name2 = "d1nklA_";
			/*cps = new int[] {
					//41, // CECP optimum
					19,59, // unpermuted local minima in TM-score
					//39, // TM-score optimum
					0,
			};*/
			
			//1itb selfsymmetry
			//name1 = "1ITB.A";
			//name2 = "1ITB.A";
			//cps = new int[] {92};
			

			OptimalCECPMain ce = (OptimalCECPMain) StructureAlignmentFactory.getAlgorithm(OptimalCECPMain.algorithmName);
			CeParameters params = (CeParameters) ce.getParameters();
			ce.setParameters(params);
			
			AtomCache cache = new AtomCache();

			Atom[] ca1 = cache.getAtoms(name1);
			Atom[] ca2 = cache.getAtoms(name2);

			AFPChain afpChain;
			

			
			// find optimal solution			
			AFPChain[] alignments = new AFPChain[ca2.length];
			afpChain = ce.alignOptimal(ca1, ca2, params, alignments);
			System.out.format("Optimal Score: %.2f\n", afpChain.getAlignScore());
			
			System.out.println("Pos\tScore\tTMScore\tLen\tRMSD\tBlocks");
			for(int i = 0; i< alignments.length; i++) {
				double tm = AFPChainScorer.getTMScore(alignments[i], ca1, ca2);
				System.out.format("%d\t%.2f\t%.2f\t%d\t%.2f\t%d\n",
						i,
						alignments[i].getAlignScore(),
						tm,
						alignments[i].getOptLength(),
						alignments[i].getTotalRmsdOpt(),
						alignments[i].getBlockNum()
				);
			}	
			
			//displayAlignment(afpChain,ca1,ca2);
			
			// permuted alignment
			for(int cp : cps) {
				// new copy of ca2, since alignPermuted has side effects
				//Atom[] ca2clone = cache.getAtoms(name2);
				//afpChain = ce.alignPermuted(ca1, ca2clone, params, cp);
				//displayAlignment(afpChain, ca1, ca2);
				
				displayAlignment(alignments[cp],ca1,ca2);
			}
			
			// CECP alignment
			CeCPMain cecp = new CeCPMain();
			afpChain = cecp.align(ca1, ca2);
			displayAlignment(afpChain,ca1,ca2);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Try showing a the afpChain in a GUI.
	 * 
	 * <p>requires additional dependencies biojava3-structure-gui and JmolApplet
	 * 
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @throws ClassNotFoundException
	 * @throws NoSuchMethodException
	 * @throws InvocationTargetException
	 * @throws IllegalAccessException
	 * @throws StructureException
	 */
	private static void displayAlignment(AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws ClassNotFoundException, NoSuchMethodException, InvocationTargetException, IllegalAccessException, StructureException {
		Atom[] ca1clone = StructureTools.cloneCAArray(ca1);
		Atom[] ca2clone = StructureTools.cloneCAArray(ca2);
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
