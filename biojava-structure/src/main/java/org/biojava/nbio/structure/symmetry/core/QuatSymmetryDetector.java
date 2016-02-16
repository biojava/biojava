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
 * Created on 2013-05-23
 *
 */
package org.biojava.nbio.structure.symmetry.core;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.symmetry.utils.CombinationGenerator;
import org.jgrapht.UndirectedGraph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.UndirectedSubgraph;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Point3d;

import java.math.BigInteger;
import java.util.*;

/**
 * Detects global and local quaternary protein structure symmetry in a structure.
 * 
 * The QuatSymmetryParameter settings affect the calculated results.
 * 
 * @author Peter Rose
 *
 */
public class QuatSymmetryDetector {
	
	private static final Logger logger = LoggerFactory
			.getLogger(QuatSymmetryDetector.class);
	
	private Structure structure = null;
	private QuatSymmetryParameters parameters = null;

	private List<QuatSymmetryResults> globalSymmetry = new ArrayList<QuatSymmetryResults>();
	private List<List<QuatSymmetryResults>> localSymmetries = new ArrayList<List<QuatSymmetryResults>>();
	private int proteinChainCount = 0;
	private boolean complete = false;

	public QuatSymmetryDetector(Structure structure, QuatSymmetryParameters parameters) {
		this.structure = structure;
		this.parameters = parameters;
	}
	
	/**
	 * Returns true if structure contains protein subunits. The other methods
	 * in this class will only return data if protein subunits are present.
	 * Always use this method first, before retrieving global or local symmetry 
	 * results.
	 * 
	 * @return true if protein subunits are present
	 */
	public boolean hasProteinSubunits() {
		run();
		return proteinChainCount > 0;
	}
	
	/**
	 * Returns list of global quaternary structure symmetry results
	 * 
	 * @return list of global quaternary structure symmetry results
	 */
	public List<QuatSymmetryResults> getGlobalSymmetry() {
		run();
		return globalSymmetry;
	}
	
	/**
	 * Returns a list of lists of local quaternary structure symmetry results
	 * 
	 * @return list of lists of local quaternary structure symmetry results
	 */
	public List<List<QuatSymmetryResults>> getLocalSymmetries() {
		run();
		return localSymmetries;
	}
	
	private void run() {
		if (complete) {
			return;
		}
		complete = true;
		//Cluster chains by sequence
		ClusterProteinChains clusterer = new ClusterProteinChains(structure, parameters);
		proteinChainCount = clusterer.getProteinChainCount();
		
		if (! hasProteinSubunits()) {
			return;
		}
		
		int nucleicAcidChainCount = clusterer.getNucleicAcidChainCount();
		
		// sort seq. identity thresholds from smallest to largest. This reduces the total number of calculations necessary.
		double[] thresholds = parameters.getSequenceIdentityThresholds().clone();
		Arrays.sort(thresholds);
		
		for (int index = 0; index < thresholds.length; index++) {
			// Map structure info to the sequence cluster
			ChainClusterer chainClusterer = new ChainClusterer(clusterer.getSequenceAlignmentClusters(thresholds[index]));
			
			// determine global symmetry
			Subunits globalSubunits = createGlobalSubunits(chainClusterer, nucleicAcidChainCount);
			QuatSymmetryResults gSymmetry = calcQuatSymmetry(globalSubunits);
			gSymmetry.setSequenceIdentityThreshold(thresholds[index]);		
			globalSymmetry.add(gSymmetry);
//			SymmetryDeviation sd = new SymmetryDeviation(globalSubunits, gSymmetry.getRotationGroup());
	

			// determine local symmetry if global structure is 
			// (1) asymmetric (C1)
			// (2) heteromeric (belongs to more than 1 sequence cluster)
			// (3) more than 2 chains (heteromers with just 2 chains cannot have local symmetry)
			
			// TODO example 2PT7: global C2, but local C6 symm., should that be included here ...?
			// i.e., include all heteromers here, for example if higher symmetry is possible by stoichiometry? A6B2 -> local A6  can have higher symmetry
			if (parameters.isLocalSymmetry() && globalSubunits.getSubunitCount() <= parameters.getMaximumLocalSubunits()) {
				if (gSymmetry.getSymmetry().equals("C1") && proteinChainCount > 2) {
					List<QuatSymmetryResults> lSymmetry = new ArrayList<QuatSymmetryResults>();
					
					long start = System.nanoTime();

					for (Subunits subunits: createLocalSubunits(chainClusterer)) {
						QuatSymmetryResults result = calcQuatSymmetry(subunits);
						addToLocalSymmetry(result, lSymmetry);
						
						double time = (System.nanoTime()- start)/1000000000;
						if (time > parameters.getLocalTimeLimit()) {
							logger.warn("Exceeded time limit for local symmetry calculations: " + time +
									" seconds. Quat symmetry results may be incomplete");
							break;
						}
					}
					localSymmetries.add(lSymmetry);
				}
			}
			
			if (! gSymmetry.getSubunits().isPseudoStoichiometric()) {
				break;
			}
		}
		

		trimGlobalSymmetryResults();
		trimLocalSymmetryResults();
		setPseudoSymmetry();
		setPreferredResults();
	}
	
	
	/**
	 * trims asymmetric global symmetry results that are C1 and have pseudoStoichiometry
	 */
	private void trimGlobalSymmetryResults() {
		for (Iterator<QuatSymmetryResults> iter = globalSymmetry.iterator(); iter.hasNext();) {
			QuatSymmetryResults result = iter.next();
			if (result.getSymmetry().equals("C1") && result.getSubunits().isPseudoStoichiometric()) {
				iter.remove();
			}
		}
	}
	
	/**
	 * trims local symmetry results if any global symmetry is found. This only happens in special cases.
	 */
	private void trimLocalSymmetryResults() {
		boolean hasGlobalSymmetry = false;
		for (QuatSymmetryResults result: globalSymmetry) {
			if (! result.getSymmetry().equals("C1")) {
				hasGlobalSymmetry = true;
				break;
			}
		} 
		
		if (hasGlobalSymmetry) {
			localSymmetries.clear();
		}
	}
	
	/**
	 * Sets pseudosymmetry flag for results that have pseudosymmetry
	 * @param thresholds sequence identity thresholds
	 */
	private void setPseudoSymmetry() {
		setPseudoSymmetry(globalSymmetry);
		for (List<QuatSymmetryResults> localSymmetry: localSymmetries) {
			setPseudoSymmetry(localSymmetry);
		}
	}
	
	/**
	 * Sets pseudosymmetry based on the analysis of all symmetry results
	 */
	private void setPseudoSymmetry(List<QuatSymmetryResults> results) {
		// find maximum symmetry order for cases of non-pseudostoichiometry
		int maxOrder = 0;
		for (QuatSymmetryResults result: results) {
			if (result.getRotationGroup() != null && !result.getSubunits().isPseudoStoichiometric()) {
				if (result.getRotationGroup().getOrder() > maxOrder) {
					maxOrder = result.getRotationGroup().getOrder();
				}
			} 
		}

		// if the order of symmetry for a pseudstoichiometric case is higher 
		// than the order for a non-pseudostoichiometry, set the pseudoSymmetry flag to true (i.e. 4HHB)
		for (QuatSymmetryResults result: results) {
			if (result.getRotationGroup() != null) {
				if (result.getRotationGroup().getOrder() > maxOrder) {
					result.getSubunits().setPseudoSymmetric(true);
				}
			}
		}
	}
	
	/**
	 * Sets preferred results flag for symmetry result that should be shown by default in visualization programs
	 * @param thresholds sequence identity thresholds
	 */
	private void setPreferredResults() {
		int[] score = new int[globalSymmetry.size()];
		
		// check global symmetry
		for (int i = 0; i < globalSymmetry.size(); i++) {
			QuatSymmetryResults result = globalSymmetry.get(i);
			if (! result.getSymmetry().equals("C1")) {
				score[i] += 2;
			}
			if (! result.getSubunits().isPseudoStoichiometric()) {
				score[i]++;
			}
		}

		int bestGlobal = 0;
		int bestScore = 0;
		for (int i = 0; i < score.length; i++) {
			if (score[i] > bestScore) {
				bestScore = score[i];
				bestGlobal = i;
			}
		}
		if (bestScore >= 2) {
			QuatSymmetryResults g = globalSymmetry.get(bestGlobal);
			g.setPreferredResult(true);
			return;
		}

		// check local symmetry
		score = new int[localSymmetries.size()];

		for (int i = 0; i < localSymmetries.size(); i++) {
			List<QuatSymmetryResults> results = localSymmetries.get(i);
			for (QuatSymmetryResults result: results) {
				if (! result.getSymmetry().equals("C1")) {
					score[i] += 2;
				}
				if (! result.getSubunits().isPseudoStoichiometric()) {
					score[i]++;
				}
			}
		}
	
		int bestLocal = 0;
		bestScore = 0;
		for (int i = 0; i < score.length; i++) {
			if (score[i] > bestScore) {
				bestScore = score[i];
				bestLocal = i;
			}
		}
		if (bestScore > 0) {
			List<QuatSymmetryResults> results = localSymmetries.get(bestLocal);
			for (QuatSymmetryResults result: results) {
				result.setPreferredResult(true);
			}
		} else {
			QuatSymmetryResults g = globalSymmetry.get(bestGlobal);
			g.setPreferredResult(true);
		}
	}
	
	private void addToLocalSymmetry(QuatSymmetryResults testResults, List<QuatSymmetryResults> localSymmetry) {
		if (testResults.getSymmetry().equals("C1")) {
			return;
		}

		for (QuatSymmetryResults results: localSymmetry) {
			if (results.getSubunits().overlaps(testResults.getSubunits())) {
				return;
			}
		}
		testResults.setLocal(true);
		localSymmetry.add(testResults);
	}

	private Subunits createGlobalSubunits(ChainClusterer chainClusterer, int nucleicAcidChainCount) {
		Subunits subunits = new Subunits(chainClusterer.getCalphaCoordinates(), 
				chainClusterer.getSequenceClusterIds(),
				chainClusterer.getPseudoStoichiometry(),
				chainClusterer.getMinSequenceIdentity(),
				chainClusterer.getMaxSequenceIdentity(),
				chainClusterer.getFolds(),
				chainClusterer.getChainIds(),
				chainClusterer.getModelNumbers());
		subunits.setNucleicAcidChainCount(nucleicAcidChainCount);
		return subunits;
	}
	
	private List<Subunits> createLocalSubunits(ChainClusterer chainClusterer) {
		List<Subunits> subunits = new ArrayList<Subunits>();
		List<List<Integer>> subClusters = decomposeClusters(chainClusterer.getCalphaCoordinates(), chainClusterer.getSequenceClusterIds());
		for (List<Integer> subCluster: subClusters) {
			subunits.add(createLocalSubunit(subCluster, chainClusterer));
		}
		return subunits;
	}
	
	private Subunits createLocalSubunit(List<Integer> subCluster, ChainClusterer chainClusterer) {
	      List<Point3d[]> subCalphaCoordinates = new ArrayList<Point3d[]>(subCluster.size());   
	      List<Integer> subSequenceIds = new ArrayList<Integer>(subCluster.size());
	      List<Boolean> subPseudoStoichiometry = new ArrayList<Boolean>(subCluster.size());
	      List<Double> subMinSequenceIdentity = new ArrayList<Double>(subCluster.size());
	      List<Double> subMaxSequenceIdentity = new ArrayList<Double>(subCluster.size());
	      List<String> subChainIds = new ArrayList<String>(subCluster.size());
	      List<Integer> subModelNumbers = new ArrayList<Integer>(subCluster.size());

	      for (int index: subCluster) {
	    	  subCalphaCoordinates.add(chainClusterer.getCalphaCoordinates().get(index));
	    	  subSequenceIds.add(chainClusterer.getSequenceClusterIds().get(index));
	    	  subPseudoStoichiometry.add(chainClusterer.getPseudoStoichiometry().get(index));
	    	  subMinSequenceIdentity.add(chainClusterer.getMinSequenceIdentity().get(index));
	    	  subMaxSequenceIdentity.add(chainClusterer.getMaxSequenceIdentity().get(index));
	    	  subChainIds.add(chainClusterer.getChainIds().get(index));
	    	  subModelNumbers.add(chainClusterer.getModelNumbers().get(index));
	      }

	      standardizeSequenceIds(subSequenceIds);
	      
	      Integer[] array = subSequenceIds.toArray(new Integer[subSequenceIds.size()]);
	      List<Integer> subFolds = getFolds(array, subSequenceIds.size());
	      Subunits subunits = new Subunits(subCalphaCoordinates, 
					subSequenceIds,
					subPseudoStoichiometry,
					subMinSequenceIdentity,
					subMaxSequenceIdentity,
			        subFolds,
					subChainIds,
					subModelNumbers);
			return subunits;
	}

	/**
	 * Resets list of arbitrary sequence ids into integer order: 0, 1, ...
	 * @param subSequenceIds
	 */
	private static void standardizeSequenceIds(List<Integer> subSequenceIds) {
		int count = 0;
	      int current = subSequenceIds.get(0);
	      for (int i = 0; i < subSequenceIds.size(); i++) {
	    	  if (subSequenceIds.get(i) > current) {
	    		  current = subSequenceIds.get(i);
	    		  count++;
	    	  }
	    	  subSequenceIds.set(i, count);
	      }
	}
	
	private List<List<Integer>> decomposeClusters(List<Point3d[]> caCoords, List<Integer> clusterIds) {
		List<List<Integer>> subClusters = new ArrayList<List<Integer>>();

		int last = getLastMultiSubunit(clusterIds);
		List<Point3d[]> subList = caCoords;
		if (last < caCoords.size()) {
			subList = caCoords.subList(0, last);
		} else {
			last = caCoords.size();
		}

		SubunitGraph subunitGraph = new SubunitGraph(subList);
		UndirectedGraph<Integer, DefaultEdge> graph = subunitGraph.getProteinGraph();
		logger.debug("Graph: {}", graph.toString());

		for (int i = last; i > 1; i--) {
			CombinationGenerator generator = new CombinationGenerator(last, i);
			int[] indices = null;
			Integer[] subCluster = new Integer[i];
			
			// avoid combinatorial explosion, i.e. for 1FNT
			BigInteger maxCombinations = BigInteger.valueOf(parameters.getMaximumLocalCombinations());
			logger.debug("Number of combinations: {}", generator.getTotal());
		    if (generator.getTotal().compareTo(maxCombinations) > 0) {
		    	logger.warn("Number of combinations exceeds limit for biounit with {} subunits in groups of {} subunits. Will not check local symmetry for them", last, i);
		    	continue;
		    }
			
			while (generator.hasNext()) {
				indices = generator.getNext();	
				
				
				// only consider sub clusters that can have rotational symmetry based on the number of subunits
				// TODO this however may exclude a few cases of helical symmetry. Checking all combinations for
				// helical symmetry would be prohibitively slow.
				for (int j = 0; j < indices.length; j++) {
					subCluster[j] = clusterIds.get(indices[j]);
				}	
				List<Integer> folds = getFolds(subCluster, last);

				if (folds.size() < 2) {
					continue;
				}
				
				List<Integer> subSet = new ArrayList<Integer>(indices.length);
				for (int index: indices) {
					subSet.add(index);
				}

				// check if this subset of subunits interact with each other
				UndirectedGraph<Integer, DefaultEdge> subGraph = new UndirectedSubgraph<Integer, DefaultEdge>(graph, new HashSet<Integer>(subSet), null);
				if (isConnectedGraph(subGraph)) {
					subClusters.add(subSet);
					if (subClusters.size() > parameters.getMaximumLocalResults()) {
						return subClusters;
					}
				}
			}
		}

		return subClusters;
	}
	
	private static int getLastMultiSubunit(List<Integer> clusterIds) {
		for (int i = 0, n = clusterIds.size(); i < n; i++) {
			if (i < n-2) {
				if (clusterIds.get(i)!=clusterIds.get(i+1) && 
						clusterIds.get(i+1) != clusterIds.get(i+2)) {
					return i+1;
				}
			}
			if (i == n-2) {
				if (clusterIds.get(i)!=clusterIds.get(i+1)) {
					return i+1;
				}
			}
		}
		return clusterIds.size();
	}
	
	private static boolean isConnectedGraph(UndirectedGraph<Integer, DefaultEdge> graph) {
		ConnectivityInspector<Integer, DefaultEdge> inspector = new ConnectivityInspector<Integer, DefaultEdge>(graph);
		return inspector.isGraphConnected();
	}
	
	private static List<Integer> getFolds(Integer[] subCluster, int size) {
		List<Integer> denominators = new ArrayList<Integer>();
		int[] counts = new int[size];
		for (int element: subCluster) {
			counts[element]++;
		}

		for (int d = 1; d <= subCluster.length; d++) {
			int count = 0;
			for (int i = 0; i < size; i++) {
				if (counts[i] > 0 && (counts[i] % d == 0)) {
					count += counts[i];
				}
			}
			if (count == subCluster.length) {
				denominators.add(d);
			}
		}
		
		Collections.sort(denominators);
		return denominators;
	}
	
	private QuatSymmetryResults calcQuatSymmetry(Subunits subunits){
		return calcQuatSymmetry(subunits, parameters);
	}
	
	public static QuatSymmetryResults calcQuatSymmetry(Subunits subunits,
			QuatSymmetryParameters parameters) {
		if (subunits.getSubunitCount() == 0) {
			return null;
		}
		
		RotationGroup rotationGroup = null;
		String method = null;
		if (subunits.getFolds().size() == 1) {			
			// no symmetry possible, create empty ("C1") rotation group
			method = "norotation";
			rotationGroup =  new RotationGroup();
			rotationGroup.setC1(subunits.getSubunitCount());
		} else if (subunits.getSubunitCount() == 2 && subunits.getFolds().contains(2)) {
			method = "C2rotation";
			QuatSymmetrySolver solver = new C2RotationSolver(subunits, parameters);
			rotationGroup = solver.getSymmetryOperations();
		} else {
			method = "rotation";
			QuatSymmetrySolver solver = new RotationSolver(subunits, parameters);
			rotationGroup = solver.getSymmetryOperations();
		}
		
		QuatSymmetryResults results = new QuatSymmetryResults(subunits, rotationGroup, method);
		
		// asymmetric structures cannot be pseudosymmetric
		String symmetry = results.getSymmetry();
		if (symmetry.equals("C1")) {
			subunits.setPseudoSymmetric(false);
		}
		
		// Check structures with Cn symmetry (n = 1, ...) for helical symmetry		
		if (symmetry.startsWith("C")) {			
			HelixSolver hc = new HelixSolver(subunits, rotationGroup.getOrder(), parameters);
			HelixLayers helixLayers = hc.getSymmetryOperations();

			if (helixLayers.size() > 0) {		
				// if the complex has no symmetry (c1) or has Cn (n>=2) symmetry and the 
				// helical symmetry has a lower RMSD than the cyclic symmetry, set helical symmetry
				// If the RMSD for helical and cyclic symmetry is similar, a slight preference is
				// given to the helical symmetry by the helixRmsdThreshold parameter.
				double cRmsd = rotationGroup.getScores().getRmsd();
				double hRmsd = helixLayers.getScores().getRmsd();
//				System.out.println("cRMSD: " + cRmsd + " hRMSD: " + hRmsd);
				double deltaRmsd = hRmsd - cRmsd;
				if (symmetry.equals("C1") || 
						(!symmetry.equals("C1") && deltaRmsd <= parameters.getHelixRmsdThreshold())) {
					method = "rottranslation";
					results = new QuatSymmetryResults(subunits, helixLayers, method);
				}
			}
		}

		return results;
	}
}
