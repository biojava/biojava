package org.biojava.bio.structure.contact;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.asa.AsaCalculator;
import org.biojava3.core.util.SingleLinkageClusterer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * A list of interfaces between 2 molecules (2 sets of atoms) 
 * 
 * @author duarte_j
 *
 */
public class StructureInterfaceList implements Serializable, Iterable<StructureInterface> {
	
	private static final Logger logger = LoggerFactory.getLogger(StructureInterfaceList.class);
	
	/**
	 * Default minimum area for a contact between two chains to be considered a
	 * valid interface.
	 * @see #removeInterfacesBelowArea(double);
	 */
	public static final double DEFAULT_MINIMUM_INTERFACE_AREA = 35.0;
	/**
	 * Default number of points to use when calculating ASAs
	 * @see #calcAsas(int, int, int)
	 */
	public static final int DEFAULT_ASA_SPHERE_POINTS = 3000;
	/**
	 * Default minimum size of cofactor molecule (non-chain HET atoms) that will be used
	 * @see #calcAsas(int, int, int)
	 */
	public static final int DEFAULT_MIN_COFACTOR_SIZE = 40;

	/**
	 * Any 2 interfaces with contact overlap score larger than this value
	 * will be considered to be clustered
	 */
	public static final double CONTACT_OVERLAP_SCORE_CLUSTER_CUTOFF = 0.2;
	
	private static final long serialVersionUID = 1L;

	private List<StructureInterface> list;

	private List<StructureInterfaceCluster> clusters;


	public StructureInterfaceList() {
		this.list = new ArrayList<StructureInterface>();
	}

	public void add(StructureInterface interf) {
		this.list.add(interf);
	}

	public int size() {
		return this.list.size();
	}

	/**
	 * Gets the interface corresponding to given id.
	 * The ids go from 1 to n
	 * If {@link #sort()} was called then the order is descendent by area.
	 * @param id
	 * @return
	 */
	public StructureInterface get(int id) {
		return list.get(id-1);
	}

	/**
	 * Calculates ASAs for all interfaces in list, both for the unbound 
	 * chains and for the complex of the two chains together. 
	 * Also sorts the interfaces based on calculated BSA areas (descending).
	 * 
	 * <p>Uses default parameters
	 */
	public void calcAsas() {
		calcAsas( DEFAULT_ASA_SPHERE_POINTS,
				Runtime.getRuntime().availableProcessors(),
				DEFAULT_MIN_COFACTOR_SIZE );
	}
	/**
	 * Calculates ASAs for all interfaces in list, both for the unbound 
	 * chains and for the complex of the two chains together. 
	 * Also sorts the interfaces based on calculated BSA areas (descending) 
	 * @param nSpherePoints
	 * @param nThreads
	 * @param cofactorSizeToUse the minimum size of cofactor molecule (non-chain HET atoms) that will be used
	 */
	public void calcAsas(int nSpherePoints, int nThreads, int cofactorSizeToUse) {

		// asa/bsa calculation 
		// NOTE in principle it is more efficient to calculate asas only once per unique chain
		// BUT! the rolling ball algorithm gives slightly different values for same molecule in different 
		// rotations (due to sampling depending on orientation of axes grid). 
		// Both NACCESS and our own implementation behave like that.
		// That's why we calculate ASAs for each rotation-unique molecule, otherwise 
		// we get discrepancies (not very big but annoying) which lead to things like negative (small) bsa values


		Map<String, Atom[]> uniqAsaChains = new TreeMap<String, Atom[]>();
		Map<String, double[]> chainAsas = new TreeMap<String, double[]>();

		// first we gather rotation-unique chains (in terms of AU id and transform id)
		for (StructureInterface interf:list) {
			String molecId1 = interf.getMoleculeIds().getFirst()+interf.getTransforms().getFirst().getTransformId();
			String molecId2 = interf.getMoleculeIds().getSecond()+interf.getTransforms().getSecond().getTransformId();

			uniqAsaChains.put(molecId1, interf.getFirstAtomsForAsa(cofactorSizeToUse)); 
			uniqAsaChains.put(molecId2, interf.getSecondAtomsForAsa(cofactorSizeToUse));
		}

		long start = System.currentTimeMillis();

		// we only need to calculate ASA for that subset (any translation of those will have same values)
		for (String molecId:uniqAsaChains.keySet()) {

			AsaCalculator asaCalc = new AsaCalculator(uniqAsaChains.get(molecId), 
					AsaCalculator.DEFAULT_PROBE_SIZE, nSpherePoints, nThreads);

			double[] atomAsas = asaCalc.calculateAsas();			

			chainAsas.put(molecId, atomAsas);

		}
		long end = System.currentTimeMillis();

		logger.debug("Calculated uncomplexed ASA for "+uniqAsaChains.size()+" orientation-unique chains. "
					+ "Time: "+((end-start)/1000.0)+" s");

		start = System.currentTimeMillis();

		// now we calculate the ASAs for each of the complexes 
		for (StructureInterface interf:list) {

			String molecId1 = interf.getMoleculeIds().getFirst()+interf.getTransforms().getFirst().getTransformId();
			String molecId2 = interf.getMoleculeIds().getSecond()+interf.getTransforms().getSecond().getTransformId();

			interf.setAsas(chainAsas.get(molecId1), chainAsas.get(molecId2), nSpherePoints, nThreads, cofactorSizeToUse);

		}
		end = System.currentTimeMillis();

		logger.debug("Calculated complexes ASA for "+list.size()+" pairwise complexes. "
					+ "Time: "+((end-start)/1000.0)+" s");


		// finally we sort based on the ChainInterface.comparable() (based in interfaceArea)
		sort();
	}

	/**
	 * Sorts the interface list and reassigns ids based on new sorting
	 */
	public void sort() {
		Collections.sort(list);
		int i=1;
		for (StructureInterface interf:list) {
			interf.setId(i);
			i++;
		}
	}

	/**
	 * Removes from this interface list all interfaces with areas
	 * below the default cutoff area
	 * @see #DEFAULT_MINIMUM_INTERFACE_AREA
	 */
	public void removeInterfacesBelowArea() {
		removeInterfacesBelowArea(DEFAULT_MINIMUM_INTERFACE_AREA);
	}

	/**
	 * Removes from this interface list all interfaces with areas
	 * below the given cutoff area
	 * @param area
	 */
	public void removeInterfacesBelowArea(double area) {
		Iterator<StructureInterface> it = iterator();
		while (it.hasNext()) {
			StructureInterface interf = it.next();
			if (interf.getTotalArea()<area) {
				it.remove();
			}
		}
	}
	
	/**
	 * Calculate the interface clusters for this StructureInterfaceList 
	 * using a contact overlap score to measure the similarity of interfaces.
	 * Subsequent calls will use the cached value without recomputing the clusters.
	 * @return
	 */
	public List<StructureInterfaceCluster> getClusters() {
		if (clusters!=null) {
			return clusters;
		}
		
		clusters = new ArrayList<StructureInterfaceCluster>();
		
		double[][] matrix = new double[list.size()][list.size()];
		
		for (int i=0;i<list.size();i++) {
			for (int j=i+1;j<list.size();j++) {
				StructureInterface iInterf = list.get(i);
				StructureInterface jInterf = list.get(j);
				
				double scoreDirect = iInterf.getContactOverlapScore(jInterf, false);
				double scoreInvert = iInterf.getContactOverlapScore(jInterf, true);
				
				double maxScore = Math.max(scoreDirect, scoreInvert);
				
				matrix[i][j] = maxScore;
			}
			
		}
		
		SingleLinkageClusterer slc = new SingleLinkageClusterer(matrix, true);
		
		Map<Integer,Set<Integer>> clusteredIndices = slc.getClusters(CONTACT_OVERLAP_SCORE_CLUSTER_CUTOFF);
		for (int clusterIdx:clusteredIndices.keySet()) {
			List<StructureInterface> members = new ArrayList<StructureInterface>();
			for (int idx:clusteredIndices.get(clusterIdx)) {
				members.add(list.get(idx));
			}
			StructureInterfaceCluster cluster = new StructureInterfaceCluster();			
			cluster.setMembers(members);
			clusters.add(cluster);
		}
		
		// finally we have to set the back-references in each StructureInterface
		for (StructureInterfaceCluster cluster:clusters) {
			for (StructureInterface interf:cluster.getMembers()) {
				interf.setCluster(cluster);
			}
		}		
		
		return clusters;
	}

	@Override
	public Iterator<StructureInterface> iterator() {
		return list.iterator();
	}
	
	@Override
	public String toString() {
		return list.toString();
	}

}
