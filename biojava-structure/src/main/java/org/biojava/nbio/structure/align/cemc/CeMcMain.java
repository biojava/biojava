package org.biojava.nbio.structure.align.cemc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.CallableStructureAlignment;
import org.biojava.nbio.structure.align.MultipleStructureAligner;
import org.biojava.nbio.structure.align.ce.CeCPMain;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsemble;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsembleImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.multiple.MultipleSuperimposer;
import org.biojava.nbio.structure.align.multiple.ReferenceSuperimposer;

/** 
 * The main class of the Java implementation of the Combinatorial Extension - Monte Carlo (CEMC) Algorithm,
 * as it was originally described by C.Guda, E.D.Scheeff, P.E. Bourne and I.N. Shindyalov (2001).
 * <p>
 * The original CEMC paper is available from <a href="http://psb.stanford.edu/psb-online/proceedings/psb01/guda.pdf">here</a>.
 * <p>
 * The usage follows the {@link MultipleStructureAligner} interface.
 * A Demo on how to use the algorithm can be found in {@link DemoCEMC}.
 * 
 * @author Aleix Lafita
 *
 */
public class CeMcMain implements MultipleStructureAligner {
	
	/**
	 *  Version history:<p>
	 *  1.0 - Initial code implementation from CEMC article without partial gaps.<p>
	 *  2.0 - Update to support CP and partial gaps.<p>
	 */
	public static final String version = "2.0";
	public static final String algorithmName = "jCEMC";
	
	private CeMcParameters params;
	private MultipleAlignmentEnsemble ensemble;
	int reference = 0;
	
	/**
	 * Default constructor. 
	 * Default parameters are used.
	 */
	public CeMcMain(){
		ensemble = null;
		params = new CeMcParameters();
	}

	/**
	 * Creates a MultipleAlignment seed for MC optimization from the representative Atoms
	 * of the structures. If there are N structures:
	 * <ul><li>Generate (N*(N-1))/2 all-to-all alignments in parallel using the Java API.
	 * <li>Choose the closest structure to all others as the reference.
	 * <li>Generate a MultipleAlignment by combining all the alignments to the reference.
	 * </ul>
	 * 
	 * @param atomArrays List of Atoms to align of the structures
	 * @return MultipleAlignment seed alignment
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 * @throws StructureException 
	 */
	private MultipleAlignment generateSeed(List<Atom[]> atomArrays) throws InterruptedException, ExecutionException, StructureException{
		
		int size = atomArrays.size();
		
		//List to store the all-to-all alignments. Contains the n^2 pairwise alignments as a 2D List
		List<List<AFPChain>> afpAlignments = new ArrayList<List<AFPChain>>();
		for (int i=0; i<size; i++){
			afpAlignments.add(new ArrayList<AFPChain>());
			for (int j=0; j<size; j++)
				afpAlignments.get(i).add(null);
		}
		
	    //Initialize the executor with a variable number of threads and the Future Object List
		ExecutorService executor = Executors.newCachedThreadPool();
	    List<Future<AFPChain>> afpFuture = new ArrayList<Future<AFPChain>>();
	    
	    //Create all the possible protein pairwise combinations (N*(N-1)/2)and call the CeCP pairwise alignment algorithm
	  	for (int i=0; i<size; i++){	  		
	  		for (int j=i+1; j<size; j++){
	    
	  			Callable<AFPChain> worker = new CallableStructureAlignment(atomArrays.get(i), atomArrays.get(j), CeCPMain.algorithmName);
	  			Future<AFPChain> submit = executor.submit(worker);
	  			afpFuture.add(submit);
	  		}
	  	}
	  	
	  	//Store the resulting AFPChains in the 2D List
	  	int index = 0;  //the alignment index
	  	for (int i=0; i<size; i++){
	  		for (int j=i; j<size; j++){
	  			if (i!=j){
		  			afpAlignments.get(i).add(j,afpFuture.get(index).get());
		  			afpAlignments.get(j).add(i,afpFuture.get(index).get());
		  			index++;
	  			}
	  		}
	  	}
	    executor.shutdown(); //Finish the executor because all the pairwise alignments have been calculated.
	    
	    reference = chooseReferenceRMSD(afpAlignments);
	    return combineReferenceAlignments(afpAlignments.get(reference), atomArrays, reference);
	}
	
	/**
	 * This method takes the all-to-all pairwise alignments Matrix (as a double List of AFPChain) and
	 * calculates the structure with the lowest average RMSD against all others. The index of this structure
	 * is returned.
	 * 
	 * @param alignments List double containing all-to-all pairwise alignments
	 * @return int reference index
	 */
	private static int chooseReferenceRMSD(List<List<AFPChain>> afpAlignments){
		
		int size = afpAlignments.size();
		
	    List<Double> RMSDs = new ArrayList<Double>();
	    for (int i=0; i<afpAlignments.size(); i++){
	    	double rmsd=0.0;
	    	for (int j=0; j<size; j++){
	    		if (i!=j) rmsd += afpAlignments.get(i).get(j).getTotalRmsdOpt();
	    	}
	    	RMSDs.add(rmsd);
	    }
	    int reference = 0;
	    for (int i=1; i<size; i++){
	    	if (RMSDs.get(i) < RMSDs.get(reference)) reference = i;
	    }
	    return reference;
	}
	
	/**
	 * This method takes a list of pairwise alignments to the reference structure and calculates a 
	 * MultipleAlignment resulting from combining their residue equivalencies.
	 * <p>
	 * It uses the blocks in AFPChain as {@link Block}s in the MultipleAlignment, so only CP are considered,
	 * and flexible parts are ignored.
	 * 
	 * @param afpList the list of pairwise alignments to the reference
	 * @param atomArrays List of Atoms of the structures
	 * @param ref index of the reference structure
	 * @return MultipleAlignment seed alignment
	 * @throws StructureException 
	 */
	private static MultipleAlignment combineReferenceAlignments(List<AFPChain> afpList, List<Atom[]> atomArrays, int ref) throws StructureException {
		
		int size = atomArrays.size();  //the number of structures
		int length = 0;  //the number of residues of the reference structure
		if (ref==0) length = afpList.get(1).getCa1Length(); //because the 0-0 alignment is null
		else length = afpList.get(0).getCa2Length();
		
		//Stores the alignment equivalencies of all the structures as a double List: equivalencies[str][res]
		List<List<Integer>> equivalencies = new ArrayList<List<Integer>>();
		for (int str=0; str<size; str++){
			equivalencies.add(new ArrayList<Integer>());
			for (int res=0; res<length; res++){
				if (str==ref) equivalencies.get(str).add(res);  //identity
				else equivalencies.get(str).add(null);
			}
		}
		
		//Now we parse the AFPChains adding the residue equivalencies
		for (int str=0; str<size; str++){
			if (str==ref) continue;  //avoid the self-comparison because it is null
			for (int bk=0; bk<afpList.get(str).getBlockNum(); bk++){
				for (int i=0; i<afpList.get(str).getOptLen()[bk]; i++){
					int res1 = 0;  //reference index
					int res2 = 0;
					//The low index is always in the first chain (0) and the higher in the second (1)
					if(str>ref){
						res1 = afpList.get(str).getOptAln()[bk][0][i];
						res2 = afpList.get(str).getOptAln()[bk][1][i];
					}
					else if (str<ref){
						res1 = afpList.get(str).getOptAln()[bk][1][i];
						res2 = afpList.get(str).getOptAln()[bk][0][i];
					}
					equivalencies.get(str).set(res1,res2);
				}
			}
		}
		
		//Now that we have translated the equivalencies we create the MultipleAlignment
		MultipleAlignment seed = new MultipleAlignmentImpl();
		seed.getEnsemble().setAtomArrays(atomArrays);
		BlockSet blockSet = new BlockSetImpl(seed);
		new BlockImpl(blockSet);
		
		//Store last positions in the block different than null to detect CP
		int[] lastResidues = new int[size];
		Arrays.fill(lastResidues, -1);
		
		//We loop through all the equivalencies checking for CP: start new Block if CP
		for (int pos=0; pos<length; pos++){
			boolean cp = false;
			for (int str=0; str<size; str++){
				if (equivalencies.get(str).get(pos) == null) continue;  //there is a gap, ignore position
				else if (equivalencies.get(str).get(pos) < lastResidues[str]){
					cp = true;  //Because the current residue is lower than the last added for at least one structure
					break;
				} else lastResidues[str] = equivalencies.get(str).get(pos);
			}
			if (cp){  //if there is a CP create a new Block, initialize AlignRes and clear lastResidues
				new BlockImpl(blockSet);
				Arrays.fill(lastResidues,-1);
			}
			
			//Now add the equivalent residues into the Block AlignRes variable
			for (int str=0; str<size; str++){
				if (blockSet.getBlocks().get(blockSet.getBlocks().size()-1).getAlignRes() == null){ //If it is empty initialize it
					List<List<Integer>> alnRes = new ArrayList<List<Integer>>(size);
					for (int k=0; k<size; k++) alnRes.add(new ArrayList<Integer>());
					blockSet.getBlocks().get(blockSet.getBlocks().size()-1).setAlignRes(alnRes);
				}
				blockSet.getBlocks().get(blockSet.getBlocks().size()-1).getAlignRes().get(str).add(equivalencies.get(str).get(pos));
			}
		}
		MultipleSuperimposer imposer= new ReferenceSuperimposer(ref);
		imposer.superimpose(seed);
		return seed;
	}

	@Override
	public MultipleAlignment align(List<Atom[]> atomArrays, Object params) throws StructureException {
		
		MultipleAlignment result = null;
		ensemble = new MultipleAlignmentEnsembleImpl();
		ensemble.setAtomArrays(atomArrays);
		
		//Generate the seed alignment from all-to-all pairwise alignments
		try {
			result = generateSeed(atomArrays);
			ExecutorService executor = Executors.newCachedThreadPool();
			List<Future<MultipleAlignment>> afpFuture = new ArrayList<Future<MultipleAlignment>>();
			int seed = 0;
			
			//Repeat the optimization 10 times in parallel, to obtain a more robust result.
			for (int i=0; i<1; i++){
				Callable<MultipleAlignment> worker = new CeMcOptimizer(result, seed+i,reference);
	  			Future<MultipleAlignment> submit = executor.submit(worker);
	  			afpFuture.add(submit);
			}

			double maxScore = Double.NEGATIVE_INFINITY;
			//When all the optimizations are finished take the one with the best result (best CEMC-Score)
			for (int i=0; i<afpFuture.size(); i++){
				MultipleAlignment align = afpFuture.get(i).get();
				double score = align.getScore(MultipleAlignmentScorer.CEMC_SCORE);
				if (score > maxScore){
					result = align;
					maxScore = score;
				}
			}
			
			ensemble.addMultipleAlignment(result);
			executor.shutdown();
			return result;
			
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}
		
		return result;
	}
	
	@Override
	public MultipleAlignment align(List<Atom[]> atomArrays) throws StructureException {
		CeMcParameters params = new CeMcParameters();
		return align(atomArrays,params);
	}

	@Override
	public ConfigStrucAligParams getParameters() {
		return params;
	}

	@Override
	public void setParameters(ConfigStrucAligParams parameters) {
		if (! (params instanceof CeMcParameters )){
			throw new IllegalArgumentException("Provided parameter object is not of type CeMcParameter");
		}
		this.params = (CeMcParameters) params;
	}

	@Override
	public String getAlgorithmName() {
		return algorithmName;
	}

	@Override
	public String getVersion() {
		return version;
	}
}
