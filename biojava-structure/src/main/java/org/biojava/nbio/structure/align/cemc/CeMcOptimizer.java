package org.biojava.nbio.structure.align.cemc;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.Block;
import org.biojava.nbio.structure.align.model.BlockImpl;
import org.biojava.nbio.structure.align.model.BlockSet;
import org.biojava.nbio.structure.align.model.BlockSetImpl;
import org.biojava.nbio.structure.align.model.MultipleAlignment;
import org.biojava.nbio.structure.align.model.Pose.PoseMethod;
import org.biojava.nbio.structure.align.model.StructureAlignmentException;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * This class takes a MultipleAlignment seed previously generated and runs a Monte Carlo optimization.
 * It implements Callable in size to be run in parallel along other optimization instances.
 * 
 * @author Aleix Lafita
 *
 */
public class CeMcOptimizer implements Callable<MultipleAlignment> {
	
	private static final boolean debug = false;  //Prints the optimization moves and saves a file of the history
	private Random rnd;
	
	//Optimization parameters
	private static final int AFPmin = 4; //Minimum block length of aligned residues
	private static final int Lmin = 8;   //Minimum alignment length
	private int iterFactor = 100; //Factor to control the max number of iterations of optimization
	private double C = 5; //Probability function constant (probability of acceptance for bad moves)
	
	//Score function parameters
	private static final double M = 20.0; //Maximum score of a match
	private static final double A = 10.0; //Penalty for alignment distances
	private double d0 = 10; //Maximum distance that is not penalized - chosen from the seed alignment
	
	//Alignment Information
	private MultipleAlignment msa;  //alignment - seed and return type
	private int size; 				//number of structures in the alignment
	private int length;				//total number of residues aligned
	
	//Multiple Alignment Residues
	private List<List<Integer>> block;     //List to store the residues aligned, in the block. Dimensions are: [size][length]
	private List<List<Integer>> freePool; 	//List to store the residues not aligned. Dimensions are: [size][residues in the pool]
	
	//Score information
	private double rmsd;     // Average RMSD of all structure superpositions
	private double tmScore;  // Average TM-score of all structure superpositions
	private double mcScore;  // Optimization score, calculated as the original CEMC algorithm
	
	//Original atom arrays of all the structures
	private List<Atom[]> atomArrays;
	private List<Integer> structureLengths;
	
	//Superposition information
	private double[] colDistances;  	  //Stores the average distance of the algined residues in a column. Entries: length.
	private double[] rowDistances;  	  //Stores the average distance from one structure to all the others. Similarity measure between structures. Entries: size.
	
	//Variables that store the history of the optimization, in size to be able to plot the evolution of the system.
	private List<Integer> lengthHistory;
	private List<Double> rmsdHistory;
	private List<Double> scoreHistory;
	
	/**
	 * Constructor.
	 * @param seedAln MultipleAlignment to be optimize.
	 * @throws StructureException 
	 * @throws StructureAlignmentException 
	 */
	public CeMcOptimizer(MultipleAlignment seedAln, long randomSeed) throws StructureException, StructureAlignmentException {
		this.msa = (MultipleAlignment) seedAln.clone();
		rnd = new Random(randomSeed);
		initialize(seedAln);
	}

	@Override
	public MultipleAlignment call() throws Exception {
		
		//The maximum number of iterations depends on the maximum possible alignment length
		optimizeMC(iterFactor*Collections.min(structureLengths));
		if (debug) saveHistory("/scratch/cemc/CeMcOptimizationHistory.csv");
		return msa;
	}
	
	/**
	 * Initialize all the variables of the optimization object.
	 * @param seed MultipleAlignment starting point.
	 * @throws StructureException
	 * @throws StructureAlignmentException 
	 */
	private void initialize(MultipleAlignment seed) throws StructureException, StructureAlignmentException {
		
		//Initialize member variables
		msa = seed;
		size = seed.size();
		length = seed.length();
		atomArrays = msa.getAtomArrays();
		structureLengths = new ArrayList<Integer>();
		for (int i=0; i<size; i++) structureLengths.add(atomArrays.get(i).length);
		
		//Initialize alignment variables
		block = new ArrayList<List<Integer>>();
		freePool = new ArrayList<List<Integer>>();
		
		//Generate the initial state of the system from the aligned blocks of the AFPChain
		for (int i=0; i<size; i++){
			ArrayList<Integer> residues = new ArrayList<Integer>();
			for (int bs=0; bs<msa.getBlockSetNum(); bs++){
				for (int b=0; b<msa.getBlockSets().get(bs).getBlockNum(); b++){
					for (int l=0; l<msa.getBlockSets().get(bs).getBlocks().get(b).length(); l++){
						Integer residue = msa.getBlockSets().get(bs).getBlocks().get(b).getAlignRes().get(i).get(l);
						residues.add(residue);
					}
				}
			}
			block.add(residues);
			freePool.add(new ArrayList<Integer>());
		}
		
		//Add any residue not aligned to the free pool for every structure
		for (int i=0; i<size; i++){
			for (int k=0; k<atomArrays.get(i).length; k++){
				if (!block.get(i).contains(k)) freePool.get(i).add(k);
			}
		}
		updateScore();
		calculatePenaltyDistance();
		updateScore();
	}
	/**
	 *  Optimization method based in a Monte-Carlo approach. 4 types of moves:
	 *  
	 *  	1- Shift Row: if there are enough freePool residues available.
	 *  	2- Expand Block: if there are enough freePool residues available.
	 *  	3- Shrink Block: move a block column to the freePool.
	 *  	4- Split and Shrink Block: split a block in the middle and shrink one column.
	 * 
	 */
	private void optimizeMC(int maxIter) throws StructureException, StructureAlignmentException{
		
		//Initialize the history variables
		lengthHistory = new ArrayList<Integer>();
		rmsdHistory = new ArrayList<Double>();
		scoreHistory = new ArrayList<Double>();
		
		int conv = 0;  //Number of steps without an alignment improvement
		int i = 1;
		
		while (i<maxIter && conv<(maxIter/20)){
			
			//Save the state of the system in case the modifications are not favorable
			List<List<Integer>> lastBlock = new ArrayList<List<Integer>>();
			List<List<Integer>> lastFreePool = new ArrayList<List<Integer>>();
			for (int k=0; k<size; k++){
				List<Integer> b = new ArrayList<Integer>();
				List<Integer> p = new ArrayList<Integer>();
				for (int l=0; l<length; l++) b.add(block.get(k).get(l));
				for (int l=0; l<freePool.get(k).size(); l++) p.add(freePool.get(k).get(l));
				lastBlock.add(b);
				lastFreePool.add(p);
			}
			double lastScore = mcScore;
			double lastRMSD = rmsd;
			double lastTMscore = tmScore;
			
			
			boolean moved = false;
			
			while (!moved){
				//Randomly select one of the steps to modify the alignment
				int move = rnd.nextInt(4);
				switch (move){
				case 0: moved = shiftRow();
						if (debug) System.out.println("did shift");
						break;
				case 1: moved = expandBlock();
						if (debug) System.out.println("did expand");
						break;
				case 2: moved = shrinkBlock();
						if (debug) System.out.println("did shrink");
						break;
				case 3: moved = splitBlock();
						if (debug) System.out.println("did split");
						break;
				}
			}
			
			//Get the properties of the new alignment
			updateScore();
			
			double AS = mcScore-lastScore;  //Change in the optimization Score
			double prob=1.0;
			
			if (AS<0){
				
				//Probability of accepting the new alignment given that produces a negative score change
				prob = probabilityFunction(AS,i,maxIter);
				double p = rnd.nextDouble();
				//Reject the move
				if (p>prob){
					block = lastBlock;
					freePool = lastFreePool;
					length = block.get(0).size();
					mcScore = lastScore;
					rmsd = lastRMSD;
					tmScore = lastTMscore;
					conv ++; //Increment the number of steps without a change in score
					
				} else conv = 0;
				
			} else conv=0;
			
			if (debug) 	System.out.println("Step: "+i+": --prob: "+prob+", --score: "+AS+", --conv: "+conv);
			
			if (i%100==1){
				lengthHistory.add(length);
				rmsdHistory.add(rmsd);
				scoreHistory.add(mcScore);
			}
			
			i++;
		}
		
		int[][][] newAlgn = new int[size][2][length];
		for (int su=0; su<size; su++){
			int[] chain1 = new int[length];
			int[] chain2 = new int[length];
			for (int k=0; k<length; k++){
				chain1[k] = block.get(su).get(k);
				chain2[k] = block.get((su+1)%size).get(k);
			}
			newAlgn[su][0] = chain1;
			newAlgn[su][1] = chain2;
		}
		
		//Override the MultipleAlignment with the optimized alignment to return
		msa.setBlockSets(new ArrayList<BlockSet>());
		BlockSet bs = new BlockSetImpl(msa);
		Block bk = new BlockImpl(bs);
		bk.setAlignRes(block);
		
		//Update Superposition (will be changed with the superposition calculated in the algorithm)
		msa.updateCache(PoseMethod.REFERENCE);
	}

	/**
	 *  Move all the block residues of one structure one position to the left or right and move the corresponding
	 *  boundary residues from the freePool to the block, and viceversa.
	 */
	private boolean shiftRow(){
		
		boolean moved = false;

		int su = rnd.nextInt(size); //Select randomly the structure that is going to be shifted
		int rl = rnd.nextInt(2);  //Select between moving right (0) or left (1)
		int res = rnd.nextInt(length); //Residue as a pivot to make the shift
		
		if (freePool.get(su).size()==0) return moved;  //If the freePool is empty the structure cannot be shifted
		
		switch(rl){
		case 0: //Move to the right
			//Check that there is at least one residue in the freePool to the left (smaller) than the pivot
			if (freePool.get(su).get(0)<block.get(su).get(res)){
				
				int leftGap = res-1;  //Find the nearest gap to the left of the res
				for (int i = res; i>=0; i--){
					if(leftGap < 0){
						break;
					} else if (block.get(su).get(i) > block.get(su).get(leftGap)+1){
						break;
					}
					leftGap--;
				}
				
				int rightGap = res+1;  //Find the nearest gap to the right of the res
				for (int i = res; i<=length; i++){
					if(rightGap == length){
						break;
					} else if (block.get(su).get(i)+1 < block.get(su).get(rightGap)){
						break;
					}
					rightGap++;
				}
				
				//Move the residue at the left of the block from the freePool to the block
				Integer residue = block.get(su).get(leftGap+1)-1;
				block.get(su).add(leftGap+1,residue);
				freePool.get(su).remove(residue);
				
				//Move the residue at the right of the block to the freePool
				if (rightGap == length){
					freePool.get(su).add(block.get(su).get(rightGap));
					block.get(su).remove(rightGap);
				}
				else{
					freePool.get(su).add(block.get(su).get(rightGap+1));
					block.get(su).remove(rightGap+1);
				}
				Collections.sort(freePool.get(su));
				
				moved = true;
			}
			break;
			
		case 1: //Move to the left
			//Check that there is at least one residue in the freePool to the right (bigger) than the pivot
			if (freePool.get(su).get(freePool.get(su).size()-1)>block.get(su).get(res)){
				
				int leftGap = res-1;  //Find the nearest gap to the left of the res
				for (int i = res; i>=0; i--){
					if(leftGap <= 0){
						leftGap = 0;
						break;
					} else if (block.get(su).get(i) > block.get(su).get(leftGap)+1){
						leftGap++;
						break;
					}
					leftGap--;
				}
				
				int rightGap = res+1;  //Find the nearest gap to the right of the res
				for (int i = res; i<=length; i++){
					if(rightGap >= length){
						rightGap = length;
						break;
					} else if (block.get(su).get(i)+1 < block.get(su).get(rightGap)){
						break;
					}
					rightGap++;
				}
				
				//Move the residue at the right of the block from the freePool to the block
				Integer residue = block.get(su).get(rightGap-1)+1;
				if (rightGap == length){
					block.get(su).add(residue);
					freePool.get(su).remove(residue);
				}
				else {
					block.get(su).add(rightGap,residue);
					freePool.get(su).remove(residue);
				}
				
				//Move the residue at the left of the block to the freePool
				freePool.get(su).add(block.get(su).get(leftGap));
				Collections.sort(freePool.get(su));
				block.get(su).remove(leftGap);
				
				moved = true;
			}
			break;
		}
		return moved;
	}
	
	/**
	 *  It extends at the beginning or end a group of consecutive aligned residues by moving the residues from the
	 *  freePool to the block.
	 */
	private boolean expandBlock(){
		
		boolean moved = false;
		
		//If any freePool is empty, the block cannot be expanded (one or more structures cannot)
		for (int su=0; su<size; su++){
			if (freePool.get(su).size()==0) return moved;
		}
			
		int rl = rnd.nextInt(2);  //Select between expanding right (0) or left (1)
		int res = rnd.nextInt(length); //Residue as a pivot to expand the structures
		
		switch (rl) {
		case 0:
			
			//Check that there is at least one residue in the freePool to the right (bigger) than the pivot in each structure
			for (int su=0; su<size; su++){
				if (freePool.get(su).get(freePool.get(su).size()-1)<block.get(su).get(res)) return moved;
			}
			
			//Find the next expandable group of residues from the pivot to the right (bigger)
			int rightRes = res+1;
			for (int i=res; i<length; i++){
				//Break if the end of the aligned residues has been found
				if (rightRes == length){
					break;
				}
				boolean groupFound = true;
				//Only a group is found when there is a gap in all the structures in that position
				for (int su=0; su<size; su++){
					if (block.get(su).get(rightRes)-1 == block.get(su).get(i)){
						groupFound = false;
						break;
					}
				}
				if (groupFound){
					break;
				}
				rightRes++;
			}
			
			//Special case: when the rightRes==length and there is no freePool residue higher, avoid adding a residue outside the structure
			for (int su=0; su<size; su++){
				if (rightRes==length && freePool.get(su).get(freePool.get(su).size()-1)<block.get(su).get(length-1)){
					return moved;
				}
			}			
			
			//Expand the block with the residues and delete them from the freePool
			for (int su=0; su<size; su++){
				Integer residue = block.get(su).get(rightRes-1)+1;
				if (rightRes == length){
					block.get(su).add(residue);
					freePool.get(su).remove(residue);
				} else{
					block.get(su).add(rightRes, residue);
					freePool.get(su).remove(residue);
				}
			}
			length++;
			moved = true;
			break;
			
		case 1:
			
			//Check that there is at least one residue in the freePool to the left (smaller) than the first block in each structure
			for (int su=0; su<size; su++){
				if (freePool.get(su).get(0)>block.get(su).get(res)) return moved;
			}
			
			//Find the next expandable group of residues from the pivot to the left (smaller)
			int leftRes = res-1;
			for (int i=res; i>0; i--){
				//Break if the start of the aligned residues has been found
				if (leftRes < 0){
					break;
				}
				boolean groupFound = true;
				//Only a group is found when there is a gap in all the structures in that position
				for (int su=0; su<size; su++){
					if (block.get(su).get(leftRes)+1 == block.get(su).get(i)){
						groupFound = false;
						break;
					}
				}
				if (groupFound){
					break;
				}
				leftRes--;
			}
			
			//Special case: when the leftRes==-1 and there is no freePool residue lower, avoid adding a residue outside the structure
			for (int su=0; su<size; su++){
				if (leftRes<0 && freePool.get(su).get(0)>block.get(su).get(0)){
					return moved;
				}
			}
			
			//Expand the block with the residues and delete them from the freePool
			for (int su=0; su<size; su++){
				Integer residue = block.get(su).get(leftRes+1)-1;
				block.get(su).add(leftRes+1,residue);
				freePool.get(su).remove(residue);
			}
			length++;
			moved = true;
			break;
		}
		return moved;
	}
	
	/**
	 *  It moves aligned residues from the start or end of a group of consecutive aligned residues 
	 *  from the block to the freePool.
	 */
	private boolean shrinkBlock(){
		
		boolean moved = false;
		if (length <= Lmin) return moved; //Do not let shrink moves if the structure is smaller than the minimum length
		
		int rl = rnd.nextInt(2);  //Select between shrink right (0) or left (1)
		int res = rnd.nextInt(length); //Residue as a pivot to shrink the structures
		
		switch (rl) {
		case 0:
			//Find the last aligned group of residues before a gap to the right (bigger)
			int rightRes = res+1;
			for (int i=res; i<length; i++){
				//Break if the end of the aligned residues has been found
				if (rightRes == length){
					break;
				}
				boolean groupFound = true;
				//Only a group is found when there is a gap in all the structures in the next position
				for (int su=0; su<size; su++){
					if (block.get(su).get(rightRes)-1 == block.get(su).get(i)){
						groupFound = false;
						break;
					}
				}
				if (groupFound){
					break;
				}
				rightRes++;
			}
			if ((rightRes-res) <= AFPmin) return moved;  //If the block (consecutive residues) is short don't shrink
			//Shrink the block and add the residues to the freePool
			for (int su=0; su<size; su++){
				Integer residue = block.get(su).get(rightRes-1);
				block.get(su).remove(rightRes-1);
				freePool.get(su).add(residue);
				Collections.sort(freePool.get(su));
			}
			length--;
			moved = true;
			break;
			
		case 1:
			//Find the first aligned group of residues after a gap to the left (smaller) than the pivot
			int leftRes = res-1;
			for (int i=res; i>0; i--){
				//Break if the start of the aligned residues has been found
				if (leftRes < 0){
					break;
				}
				boolean groupFound = true;
				//Only a group is found when there is a gap in all the structures in that position
				for (int su=0; su<size; su++){
					if (block.get(su).get(leftRes)+1 == block.get(su).get(i)){
						groupFound = false;
						break;
					}
				}
				if (groupFound){
					break;
				}
				leftRes--;
			}
			if ((res-leftRes) <= AFPmin) return moved;  //If the block (consecutive residues) is short don't shrink
			//Shrink the block and add the residues to the freePool
			for (int su=0; su<size; su++){
				Integer residue = block.get(su).get(leftRes+1);
				block.get(su).remove(leftRes+1);
				freePool.get(su).add(residue);
				Collections.sort(freePool.get(su));
			}
			length--;
			moved = true;
			break;
		}
		return moved;
	}
	
	private boolean splitBlock(){
		
		boolean moved = false;
		if (length <= Lmin) return moved; //Let split moves everywhere if the alignment is larger than the minimum length
		
		int res = rnd.nextInt(length); //Residue as a pivot to split the structures
		
		for (int su=0; su<size; su++){
			Integer residue = block.get(su).get(res);
			block.get(su).remove(res);
			freePool.get(su).add(residue);
			Collections.sort(freePool.get(su));
		}
		length--;
		moved = true;
		return moved;
	}
	
	/**
	 *  Calculates the average RMSD and Score of all structure superimpositions of the structure, corresponding to the
	 *  aligned residues in block. It also updates the Monte Carlo score, the optimized value, from the distances matrix.
	 *  It uses a different transformation for every structure pair (flexible).
	 */
	private void updateScore() throws StructureException{
		
		//The first index is the score and the second the RMSD
		double[] ScoreRMSD = {0.0,0.0};
		colDistances = new double[length];
		rowDistances = new double[size];
		
		//Reset variables
		rmsd = 0.0;
		tmScore = 0.0;
		mcScore = 0.0;
		
		//Construct all possible pairwise comparisons of the molecule (size)*(size-1)/2
		for (int i=0; i<size; i++){
			for (int j=i+1; j<size; j++){
				calculateScore(i, j, ScoreRMSD);
				rmsd += ScoreRMSD[1];
				tmScore += ScoreRMSD[0];
			}
		}
		//Divide the colDistances entries for the total number of comparisons to get the average
		int total = (size)*(size-1)/2;
		for (int i=0; i<length; i++) colDistances[i] /= total;
		for (int i=0; i<size; i++) rowDistances[i] /= (size-1)*length;
		
		//Assign the new values to the member variables
		rmsd /= total;
		tmScore /= total;
		mcScore = scoreFunctionMC();
	}
	
	/**
	 *  Raw implementation of the RMSD, TMscore and distances calculation for a better algorithm efficiency.
	 *  Superimpose two structures by their aligned residues and get the RMSD and TM-score.
	 */
	private void calculateScore(int str1, int str2, double[] ScoreRMSD) throws StructureException{
		
		Atom[] arr1 = new Atom[length];
		Atom[] arr2 = new Atom[length];
		int pos = 0;
		
		//Calculate the aligned atom arrays
		for (int k=0; k<length; k++){
			arr1[pos] = atomArrays.get(str1)[block.get(str1).get(k)];
			arr2[pos] = (Atom) atomArrays.get(str2)[block.get(str2).get(k)].clone();
			pos++;
		}
		
		//Superimpose the two structures in correspondence to the new alignment
		SVDSuperimposer svd = new SVDSuperimposer(arr1, arr2);
		Matrix matrix = svd.getRotation();
		Atom shift = svd.getTranslation();
		
		//Transform atoms of the second strcuture
		Calc.rotate(arr2, matrix);
		Calc.shift(arr2, shift);
		
		//Get the rmsd and score of the rotation
		ScoreRMSD[1] = SVDSuperimposer.getRMS(arr1, arr2);
		ScoreRMSD[0] = SVDSuperimposer.getTMScore(arr1, arr2, structureLengths.get(str1), structureLengths.get(str2));
		
		//Calculate the distances between C alpha atoms of the same column and store them in colDistances
		for (int k=0; k<arr1.length; k++){
			double distance = Math.abs(Calc.getDistance(arr1[k], arr2[k]));
			colDistances[k] += distance;
			rowDistances[str1] += distance;
			rowDistances[str2] += distance;
		}
	}
	
	/**
	 *  Calculates the optimization score from the column average distances.
	 *  
	 *  Function: sum(M/(d1/d0)^2)  and add the penalty A if (d1>d0).
	 */
	private double scoreFunctionMC(){
		
		double score = 0.0;
		
		//Loop through all the columns
		for (int col=0; col<length; col++){
			
			double d1 = colDistances[col];
			double colScore = M/Math.pow(1+d1/d0,2);
			
			if (d1>d0) colScore-=A;
			score += colScore;
		}
		return score;
	}
	
	/**
	 *  Calculates the probability of accepting a bad move given the iteration step and the score change.
	 *  
	 *  Function: p=(C-AS)/m^0.5   *from the CEMC algorithm.
	 *  Added a normalization factor so that the probability approaches 0 when the maxIter is reached.
	 */
	private double probabilityFunction(double AS, int m, int maxIter) {
		
		double prob = (C+AS)/Math.sqrt(m);
		double norm = (1-(m*1.0)/maxIter);  //Normalization factor
		return Math.min(Math.max(prob*norm,0.0),1.0);
	}
	
	/**
	 *  Calculate the maximum distance which is not penalized in the score function. Only used at the beginning.
	 *  Options are:
	 *    1- Pick the 90% of the distances range (softer condition, results in longer structures).
	 *    2- Pick the average value of the top 10% distances (can be softer than 1 depending on the range scale).
	 *    3- Pick the value at the boundary of the top 10% distances (hardest condition, restricts the structures to the core only).
	 *    4- A function of the RMSD of the seed alignment and the size (this is the softest condition of all, but longer core alignments are obtained).
	 *    		Justification: the more structures the higher the variability between the columns.
	 *    
	 *  A minimum distance of 5A is set always to avoid short alignments in the very good symmetric cases.
	 */
	private void calculatePenaltyDistance() {
		
		double[] distances = colDistances.clone();
		Arrays.sort(distances);
		
		//Option 1: 90% of distances range
		double range = distances[distances.length-1] - distances[0];
		double d1 = distances[0] + range*0.9;
		
		//Option 2: average of the top 10% distances
		int index10 = (distances.length-1) - distances.length/10;
		double d2 = 0.0;
		for (double dist:Arrays.copyOfRange(distances,index10,distances.length)) d2+=dist;
		d2 /= (distances.length-index10);
		
		//Option 3: boundary of top 10%
		double d3 = distances[index10];
		
		//Option 4: the highest distance of the seed alignment.
		double d4 = rmsd*(1+0.25*size);
		
		d0=Math.max(d4,5);  //The minimum d0 value is 5A
	}
	
	/**
	 * Save the evolution of the optimization process as a csv file.
	 */
	private void saveHistory(String filePath) throws IOException{
		
	    FileWriter writer = new FileWriter(filePath);
	    writer.append("Step,Length,RMSD,Score\n");
	    
	    for (int i=0; i<lengthHistory.size(); i++){
	    		writer.append(i*100+","+lengthHistory.get(i)+","+rmsdHistory.get(i)+","+scoreHistory.get(i)+"\n");
	    }
	    writer.flush();
	    writer.close();
	}
}
