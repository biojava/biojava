package org.biojava.nbio.structure.align.cemc;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.CallableStructureAlignment;
import org.biojava.nbio.structure.align.MultipleStructureAlignment;
import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.Block;
import org.biojava.nbio.structure.align.model.BlockImpl;
import org.biojava.nbio.structure.align.model.BlockSet;
import org.biojava.nbio.structure.align.model.BlockSetImpl;
import org.biojava.nbio.structure.align.model.MultipleAlignment;
import org.biojava.nbio.structure.align.model.MultipleAlignmentEnsemble;
import org.biojava.nbio.structure.align.model.MultipleAlignmentEnsembleImpl;
import org.biojava.nbio.structure.align.model.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.model.StructureAlignmentException;
import org.biojava.nbio.structure.align.superimpose.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.superimpose.MultipleSuperimposer;
import org.biojava.nbio.structure.align.superimpose.ReferenceSuperimposer;

/** 
 * The main class of the Java implementation of the Combinatorial Extension - Monte Carlo Algorithm (CEMC),
 * as has been originally developed by C.Guda, E.D.Scheeff, P.E. Bourne and I.N. Shindyalov (2001).
 * The original CEMC paper is available from <a href="http://psb.stanford.edu/psb-online/proceedings/psb01/guda.pdf">here</a>.
 * 
 * There is still no demo on how to use this algorithm.
 * 
 * @author Aleix Lafita
 *
 */
public class CeMcMain implements MultipleStructureAlignment{
	
	/**
	 *  version history:
	 *  1.0 - Initial code implementation from CEMC article.
	 */
	public static final String version = "1.0";
	public static final String algorithmName = "jCEMC";
	private CeMcParameters params;
	private MultipleAlignmentEnsemble ensemble;
	
	/**
	 * Default constructor. Instantiates an empty CeMcMain object with the default parameters.
	 */
	public CeMcMain(){
		ensemble = null;
		params = new CeMcParameters();
	}

	/**
	 * Creates the seed of the multiple structure alignment, before optimization, from all-to-all pairwise alignments
	 * of the structures. The alignments are generated in parallel using the Java API for concurrency management. 
	 * The closest structure to all others is chosen as the reference and all the alignments to it are taken to generate 
	 * an ungapped seed MultipleAlignment.
	 * This method is static because can be used outside this alignment class.
	 * @param atomArrays List of Atoms to align of the structures
	 * @return MultipleAlignment seed alignment
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 * @throws StructureException 
	 * @throws StructureAlignmentException 
	 */
	public MultipleAlignment generateSeed(List<Atom[]> atomArrays) throws InterruptedException, ExecutionException, StructureAlignmentException, StructureException{
		
		int size = atomArrays.size();
		
		//List to store the all-to-all alignments. Contains the n^2 pairwise alignments as a 2D matrix indices.
		List<List<AFPChain>> afpAlignments = new ArrayList<List<AFPChain>>();
		for (int i=0; i<size; i++){
			afpAlignments.add(new ArrayList<AFPChain>());
			for (int j=0; j<size; j++)
				afpAlignments.get(i).add(null);
		}
		
	    //Initialize the executor with a variable number of threads and the Future Object List
		ExecutorService executor = Executors.newCachedThreadPool();
	    List<Future<AFPChain>> afpFuture = new ArrayList<Future<AFPChain>>();
	    
	    //Create all the possible protein pairwise combinations and call the alignment
	  	for (int i=0; i<size; i++){	  		
	  		for (int j=i+1; j<size; j++){
	    
	  			Callable<AFPChain> worker = new CallableStructureAlignment(atomArrays.get(i), atomArrays.get(j), CeMain.algorithmName);
	  			Future<AFPChain> submit = executor.submit(worker);
	  			afpFuture.add(submit);
	  		}
	  	}
	  	
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
	    
	    //Choose the reference structure from the lowest average RMSD
	    List<Double> RMSDs = new ArrayList<Double>();
	    for (int i=0; i<size; i++){
	    	double rmsd=0.0;
	    	for (int j=0; j<size; j++){
	    		if (i!=j) rmsd += afpAlignments.get(i).get(j).getTotalRmsdOpt();
	    	}
	    	RMSDs.add(rmsd);
	    }
	    int ref = 0;
	    for (int i=1; i<size; i++){
	    	if (RMSDs.get(i) < RMSDs.get(ref)) ref = i;
	    }
	    
	    return seedFromReference(afpAlignments.get(ref), atomArrays, ref);
	}
	
	/**
	 * This method takes a list of pairwise alignments to the reference structure and calculates the 
	 * MultipleAlignment resulting from them. It ignores blocks in AFPChain (flexible parts) and builds 
	 * the Blocks as the definition of MultipleAlignment dictates {@link Block}. Gaps are not included.
	 * @param afpList the list of pairwise alignments to the reference
	 * @param atomArrays List of Atoms of the structures
	 * @param ref index of the reference structure
	 * @return MultipleAlignment seed alignment
	 * @throws StructureAlignmentException 
	 * @throws StructureException 
	 */
	private MultipleAlignment seedFromReference(List<AFPChain> afpList, List<Atom[]> atomArrays, int ref) throws StructureAlignmentException, StructureException {
		
		int size = atomArrays.size();  //the number of structures
		int length = 0;  //the number of residues of the reference structure
		if (ref==0) length = afpList.get(1).getCa1Length();
		else length = afpList.get(0).getCa2Length();
		
		//Stores the alignment equivalencies of all the structures as a double List (first index reference residue, second structure nr) - note this is the inverse of AlignRes in Block
		List<List<Integer>> equivalencies = new ArrayList<List<Integer>>();
		for (int i=0; i<length; i++){
			equivalencies.add(new ArrayList<Integer>());
			for (int j=0; j<size; j++){
				if (j==ref) equivalencies.get(i).add(i);
				else equivalencies.get(i).add(null);
			}
		}
		
		//Now we parse the AFPChains adding the residue equivalencies
		for (int j=0; j<size; j++){
			if (j==ref) continue;  //avoid the self-comparison because it is null
			for (int bk=0; bk<afpList.get(j).getBlockNum(); bk++){
				for (int i=0; i<afpList.get(j).getOptLen()[bk]; i++){
					int res1 = 0;
					int res2 = 0;
					//The low index is always in the first chain and the higher in the second
					if(j>ref){
						res1 = afpList.get(j).getOptAln()[bk][0][i];
						res2 = afpList.get(j).getOptAln()[bk][1][i];
					}
					else if (j<ref){
						res1 = afpList.get(j).getOptAln()[bk][1][i];
						res2 = afpList.get(j).getOptAln()[bk][0][i];
					}
					equivalencies.get(res1).set(j,res2);
				}
			}
		}
		
		//Now that we have the equivalencies we create the MultipleAlignment
		MultipleAlignment seed = new MultipleAlignmentImpl(ensemble);
		BlockSet blockSet = new BlockSetImpl(seed);
		new BlockImpl(blockSet);  //This automatically adds a Block to the Block list in BlockSet
		
		//We loop through all the equivalencies checking for gapped columns and CP
		//Case CP: we start a new block
		//Case gap: we do not include the equivalencies, the seed cannot contain gaps
		for (int i=0; i<length; i++){
			boolean gap = false;
			boolean cp = false;
			for (int j=0; j<size; j++){
				if (equivalencies.get(i).get(j) == null){  //there is a gap, do not consider the 
					gap = true;
					break;
				}
				if (blockSet.getBlocks().get(blockSet.getBlocks().size()-1).length() > 0){
					if (equivalencies.get(i).get(j) < blockSet.getBlocks().get(blockSet.getBlocks().size()-1)
							.getAlignRes().get(j).get( blockSet.getBlocks().get(blockSet.getBlocks().size()-1).length()-1)){
						cp = true;  //Because the current residue is lower than the last added for at least one structure
					}
				}
			}
			if (gap) continue;
			if (cp)  //if there is a CP create a new Block
				new BlockImpl(blockSet);
			
			//Now add the equivalent residues into the Block AlignRes variable
			for (int j=0; j<size; j++){
				if (blockSet.getBlocks().get(blockSet.getBlocks().size()-1).getAlignRes().size() == 0) //If it is empty initialize a list for every structure
					for (int k=0; k<size; k++) blockSet.getBlocks().get(blockSet.getBlocks().size()-1).getAlignRes().add(new ArrayList<Integer>());
				blockSet.getBlocks().get(blockSet.getBlocks().size()-1).getAlignRes().get(j).add(equivalencies.get(i).get(j));
			}
		}
		MultipleSuperimposer imposer= new ReferenceSuperimposer();
		imposer.superimpose(seed);
		return seed;
	}

	@Override
	public MultipleAlignment align(List<Atom[]> atomArrays, Object params) throws StructureException, StructureAlignmentException {
		
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
			for (int i=0; i<10; i++){
				Callable<MultipleAlignment> worker = new CeMcOptimizer(result, seed+i);
	  			Future<MultipleAlignment> submit = executor.submit(worker);
	  			afpFuture.add(submit);
			}

			Double resultScore = result.getScore(MultipleAlignmentScorer.SCORE_REF_TMSCORE);
			if(resultScore == null) {
				resultScore = MultipleAlignmentScorer.getRefTMScore(result,0);
				result.putScore(MultipleAlignmentScorer.SCORE_REF_TMSCORE,resultScore);
			}

			//When all the optimizations are finished take the one with the best result (best TM-score)
			for (int i=0; i<afpFuture.size(); i++){
				MultipleAlignment align = afpFuture.get(i).get();
				//TODO change TMScore? -SB
				Double tmScore = align.getScore(MultipleAlignmentScorer.SCORE_REF_TMSCORE);
				if(tmScore == null) {
					tmScore = MultipleAlignmentScorer.getRefTMScore(align,0);
					align.putScore(MultipleAlignmentScorer.SCORE_REF_TMSCORE,tmScore);
				}
				
				if (tmScore > resultScore){
					result = align;
					resultScore = tmScore;
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
	public MultipleAlignment align(List<Atom[]> atomArrays) throws StructureException, StructureAlignmentException {
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
