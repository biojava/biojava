package org.biojava.nbio.structure.align.multiple.mc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.CallableStructureAlignment;
import org.biojava.nbio.structure.align.MultipleStructureAligner;
import org.biojava.nbio.structure.align.StructureAlignment;
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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/** 
 * Main class of the Java implementation of the Combinatorial Extension - 
 * Monte Carlo (CEMC) Algorithm,
 * as it was originally described by C.Guda, E.D.Scheeff, P.E. Bourne and 
 * I.N. Shindyalov (2001).
 * The original CEMC paper is available from 
 * <a href="http://psb.stanford.edu/psb-online/proceedings/psb01/guda.pdf">
 * here</a>.
 * <p>
 * This implementation is a generalized version that allows any pairwise 
 * structure alignment algorithm as input, thus supporting any non-topological
 * or flexible alignment. The seed can also directly be the input for the 
 * optimization. For that, look at {@link MultipleMcOptimizer}.
 * <p>
 * A Demo on how to use the algorithm can be found in the demo package.
 * 
 * @author Aleix Lafita
 * @since 4.1.0
 *
 */
public class MultipleMcMain implements MultipleStructureAligner {

	private static final Logger logger = 
			LoggerFactory.getLogger(MultipleMcMain.class);

	/**
	 *  Version history:<p>
	 *  1.0 - Initial code implementation from CEMC article.<p>
	 *  1.1 - Support CP, non-topological and flexible seed alignments.<p>
	 */
	public static final String version = "1.1";
	public static final String algorithmName = "jMultipleMC";

	private MultipleMcParameters params;
	private MultipleAlignmentEnsemble ensemble;
	private StructureAlignment pairwise;
	private int reference = 0;

	/**
	 * Default constructor. 
	 * Default parameters are used.
	 * @param pairwise the pairwise structure alignment used to generate the
	 * 			multiple alignment seed.
	 */
	public MultipleMcMain(StructureAlignment pairwiseAlgo){
		ensemble = null;
		params = new MultipleMcParameters();
		pairwise = pairwiseAlgo;
		if (pairwise == null) pairwise = new CeCPMain();
	}

	/**
	 * Creates a MultipleAlignment seed for MC optimization from the 
	 * representative Atoms of the structures. If there are N structures:
	 * <ul><li>Generate (N*(N-1))/2 all-to-all alignments in parallel using 
	 * the Java API.
	 * <li>Choose the closest structure to all others as the reference.
	 * <li>Generate a MultipleAlignment by combining all the alignments to 
	 * the reference, that will be used as a seed for the MC optimization.
	 * </ul>
	 * 
	 * @param atomArrays List of Atoms to align of the structures
	 * @return MultipleAlignment seed alignment
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 * @throws StructureException 
	 */
	private MultipleAlignment generateSeed(List<Atom[]> atomArrays) 
			throws InterruptedException, 
			ExecutionException, StructureException {

		int size = atomArrays.size();

		//List to store the all-to-all alignments
		List<List<AFPChain>> afpAlignments = new ArrayList<List<AFPChain>>();
		for (int i=0; i<size; i++){
			afpAlignments.add(new ArrayList<AFPChain>());
			for (int j=0; j<size; j++)
				afpAlignments.get(i).add(null);
		}

		int threads = params.getNrThreads();
		ExecutorService executor = Executors.newFixedThreadPool(threads);
		List<Future<AFPChain>> afpFuture = new ArrayList<Future<AFPChain>>();

		//Create all the possible protein pairwise combinations 
		//(N*(N-1)/2) and call the pairwise alignment algorithm
		for (int i=0; i<size; i++){	  		
			for (int j=i+1; j<size; j++){

				Callable<AFPChain> worker = new CallableStructureAlignment(
						atomArrays.get(i), atomArrays.get(j), 
						pairwise.getAlgorithmName(), pairwise.getParameters());

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
		executor.shutdown();

		reference = chooseReferenceRMSD(afpAlignments);
		boolean flexible = false;
		if (pairwise.getAlgorithmName().contains("flexible"))
			flexible = true;

		return combineReferenceAlignments(afpAlignments.get(reference),
				atomArrays, reference, flexible);
	}

	/**
	 * This method takes the all-to-all pairwise alignments Matrix (as a 
	 * double List of AFPChain) and calculates the structure with the 
	 * lowest average RMSD against all others. 
	 * The index of this structure is returned.
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
				if (i!=j) 
					rmsd += afpAlignments.get(i).get(j).getTotalRmsdOpt();
			}
			RMSDs.add(rmsd);
		}
		int reference = 0;
		for (int i=1; i<size; i++){
			if (RMSDs.get(i) < RMSDs.get(reference)) reference = i;
		}
		logger.info("Reference structure is "+reference);
		return reference;
	}

	/**
	 * This method takes a list of pairwise alignments to the reference 
	 * structure and calculates a MultipleAlignment resulting from combining 
	 * their residue equivalencies.
	 * <p>
	 * It uses the blocks in AFPChain as {@link Block}s in the 
	 * MultipleAlignment, so considers non-topological
	 * alignments, if the alignment is rigid. If the alignment is flexible, 
	 * it considers the blocks as {@link BlockSets}.
	 * 
	 * @param afpList the list of pairwise alignments to the reference
	 * @param atomArrays List of Atoms of the structures
	 * @param ref index of the reference structure
	 * @param flexible uses BlockSets if true, uses Blocks otherwise
	 * @return MultipleAlignment seed alignment
	 * @throws StructureException 
	 */
	private static MultipleAlignment combineReferenceAlignments(
			List<AFPChain> afpList, List<Atom[]> atomArrays, 
			int ref, boolean flexible) throws StructureException {

		int size = atomArrays.size();
		int length = 0;  //the number of residues of the reference structure
		if (ref==0) length = afpList.get(1).getCa1Length();
		else length = afpList.get(0).getCa2Length();
		SortedSet<Integer> flexibleBoundaries = new TreeSet<Integer>();

		//Stores the equivalencies of all the structures as a double List
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
			if (str==ref) continue;  //avoid self-comparison
			for (int bk=0; bk<afpList.get(str).getBlockNum(); bk++){
				for (int i=0; i<afpList.get(str).getOptLen()[bk]; i++){
					int res1 = 0;  //reference index
					int res2 = 0;
					//The low index is always in the first chain (0)
					if(str>ref){
						res1 = afpList.get(str).getOptAln()[bk][0][i];
						res2 = afpList.get(str).getOptAln()[bk][1][i];
					}
					else if (str<ref){
						res1 = afpList.get(str).getOptAln()[bk][1][i];
						res2 = afpList.get(str).getOptAln()[bk][0][i];
					}
					equivalencies.get(str).set(res1,res2);

					//Add the flexible boundaries if flexible
					if(flexible && i==0) flexibleBoundaries.add(res1);
				}
			}
		}

		//We have translated the equivalencies, we create the MultipleAlignment
		MultipleAlignment seed = new MultipleAlignmentImpl();
		seed.getEnsemble().setAtomArrays(atomArrays);
		BlockSet blockSet = new BlockSetImpl(seed);
		new BlockImpl(blockSet);

		//Store last positions in the block different than null to detect CP
		int[] lastResidues = new int[size];
		Arrays.fill(lastResidues, -1);

		//We loop through all the equivalencies checking for CP
		for (int pos=0; pos<length; pos++){
			//Start a new BlockSet if the position means a boundary
			if (flexibleBoundaries.contains(pos) && 
					blockSet.getBlocks().get(blockSet.getBlocks().size()-1).
					getAlignRes() != null){

				blockSet = new BlockSetImpl(seed);
				new BlockImpl(blockSet);
			}

			boolean cp = false;
			for (int str=0; str<size; str++){
				if (equivalencies.get(str).get(pos) == null){
					continue;  //there is a gap, ignore position
				} else if (equivalencies.get(str).get(pos)<lastResidues[str]){
					cp = true;  //current residue is lower than the last
					break;
				} else lastResidues[str] = equivalencies.get(str).get(pos);
			}
			if (cp){  //if there is a CP create a new Block
				new BlockImpl(blockSet);
				Arrays.fill(lastResidues,-1);
			}

			//Now add the equivalent residues into the Block AlignRes variable
			for (int str=0; str<size; str++){
				Block lastB = blockSet.getBlocks().get(
						blockSet.getBlocks().size()-1);

				if (lastB.getAlignRes() == null){
					//Initialize the aligned residues list
					List<List<Integer>> alnRes = 
							new ArrayList<List<Integer>>(size);

					for (int k=0; k<size; k++) {
						alnRes.add(new ArrayList<Integer>());
					}
					lastB.setAlignRes(alnRes);
				}
				lastB.getAlignRes().get(str).add(
						equivalencies.get(str).get(pos));
			}
		}
		logger.info("Seed alignment has "+seed.getBlocks()+" Blocks.");
		return seed;
	}

	@Override
	public MultipleAlignment align(List<Atom[]> atomArrays, Object parameters)
			throws StructureException {

		MultipleAlignment result = null;
		ensemble = new MultipleAlignmentEnsembleImpl();
		ensemble.setAtomArrays(atomArrays);
		ensemble.setAlgorithmName(algorithmName);
		ensemble.setVersion(version);
		ensemble.setIoTime(System.currentTimeMillis());
		setParameters((ConfigStrucAligParams) parameters);

		//Generate the seed alignment and optimize it
		try {
			result = generateSeed(atomArrays);
		} catch (InterruptedException e) {
			logger.warn("Seed generation failed.",e);
		} catch (ExecutionException e) {
			logger.warn("Seed generation failed.",e);
		}

		//Repeat the optimization in parallel - DISALLOWED
		/*int threads = params.getNrThreads();
			ExecutorService executor = Executors.newFixedThreadPool(threads);
			List<Future<MultipleAlignment>> msaFuture = 
					new ArrayList<Future<MultipleAlignment>>();

			for (int i=0; i<params.getNrThreads(); i++){
				//Change the random seed for each parallelization
				MultipleMcParameters paramsMC = (MultipleMcParameters) params;
				paramsMC.setRandomSeed(paramsMC.getRandomSeed()+i);

				Callable<MultipleAlignment> worker = 
						new MultipleMcOptimizer(
								result, paramsMC, reference);

	  			Future<MultipleAlignment> submit = executor.submit(worker);
	  			msaFuture.add(submit);
			}

			double maxScore = Double.NEGATIVE_INFINITY;
			//Take the one with the best result (best MC-Score)
			for (int i=0; i<msaFuture.size(); i++){
				MultipleAlignment align = msaFuture.get(i).get();
				double s = align.getScore(MultipleAlignmentScorer.MC_SCORE);
				if (s > maxScore){
					result = align;
					maxScore = s;
				}
			}
			Long runtime = System.currentTimeMillis()-ensemble.getIoTime();
			ensemble.setCalculationTime(runtime);

			result.setEnsemble(ensemble);
			ensemble.addMultipleAlignment(result);
			executor.shutdown();*/

		MultipleMcOptimizer optimizer = new MultipleMcOptimizer(
				result, params, reference);

		Long runtime = System.currentTimeMillis()-ensemble.getIoTime();
		ensemble.setCalculationTime(runtime);

		result = optimizer.optimize();
		result.setEnsemble(ensemble);
		ensemble.addMultipleAlignment(result);

		return result;

	}

	@Override
	public MultipleAlignment align(List<Atom[]> atomArrays) 
			throws StructureException {

		if (params == null) {
			logger.info("Using DEFAULT MultipleMc Parameters");
			params = new MultipleMcParameters();
		}
		return align(atomArrays,params);
	}

	@Override
	public ConfigStrucAligParams getParameters() {
		return params;
	}

	@Override
	public void setParameters(ConfigStrucAligParams parameters) {
		if (!(parameters instanceof MultipleMcParameters)){
			throw new IllegalArgumentException(
					"Provided parameter object is not of type MultipleMC");
		}
		this.params = (MultipleMcParameters) parameters;
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
