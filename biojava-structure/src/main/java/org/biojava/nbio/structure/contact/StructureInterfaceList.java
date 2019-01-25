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
 */
package org.biojava.nbio.structure.contact;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;

import org.biojava.nbio.core.util.SingleLinkageClusterer;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.asa.AsaCalculator;
import org.biojava.nbio.structure.xtal.CrystalBuilder;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * A list of interfaces between 2 molecules (2 sets of atoms)
 *
 * @author Jose Duarte
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
	public static final double DEFAULT_CONTACT_OVERLAP_SCORE_CLUSTER_CUTOFF = 0.2;

	private static final long serialVersionUID = 1L;

	private List<StructureInterface> list;

	private List<StructureInterfaceCluster> clusters = null;
	private List<StructureInterfaceCluster> clustersNcs = null;

	private Map<String, String> chainOrigNamesMap;

	public StructureInterfaceList() {
		this.list = new ArrayList<>();
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


		Map<String, Atom[]> uniqAsaChains = new TreeMap<>();
		Map<String, double[]> chainAsas = new TreeMap<>();

		List<StructureInterface> redundancyReducedList;
		if (clustersNcs != null) {
			redundancyReducedList = new ArrayList<>();
			for (StructureInterfaceCluster ncsCluster : clustersNcs) {
				// we use the first one in list as the only one for which we calculate ASAs
				redundancyReducedList.add(ncsCluster.getMembers().get(0));
			}

		} else {
			redundancyReducedList = list;
		}

		// first we gather rotation-unique chains (in terms of AU id and transform id)
		for (StructureInterface interf:redundancyReducedList) {
			String molecId1 = interf.getMoleculeIds().getFirst()+interf.getTransforms().getFirst().getTransformId();
			String molecId2 = interf.getMoleculeIds().getSecond()+interf.getTransforms().getSecond().getTransformId();

			uniqAsaChains.put(molecId1, interf.getFirstAtomsForAsa(cofactorSizeToUse));
			uniqAsaChains.put(molecId2, interf.getSecondAtomsForAsa(cofactorSizeToUse));
		}

		logger.debug("Will calculate uncomplexed ASA for {} orientation-unique chains.", uniqAsaChains.size());

		long start = System.currentTimeMillis();

		// we only need to calculate ASA for that subset (any translation of those will have same values)
		for (String molecId:uniqAsaChains.keySet()) {

			logger.debug("Calculating uncomplexed ASA for molecId {}, with {} atoms", molecId, uniqAsaChains.get(molecId).length);

			AsaCalculator asaCalc = new AsaCalculator(uniqAsaChains.get(molecId),
					AsaCalculator.DEFAULT_PROBE_SIZE, nSpherePoints, nThreads);

			double[] atomAsas = asaCalc.calculateAsas();

			chainAsas.put(molecId, atomAsas);

		}
		long end = System.currentTimeMillis();

		logger.debug("Calculated uncomplexed ASA for {} orientation-unique chains. Time: {} s", uniqAsaChains.size(), ((end-start)/1000.0));

		logger.debug ("Will calculate complexed ASA for {} pairwise complexes.", redundancyReducedList.size());

		start = System.currentTimeMillis();

		// now we calculate the ASAs for each of the complexes
		for (StructureInterface interf:redundancyReducedList) {

			String molecId1 = interf.getMoleculeIds().getFirst()+interf.getTransforms().getFirst().getTransformId();
			String molecId2 = interf.getMoleculeIds().getSecond()+interf.getTransforms().getSecond().getTransformId();

			logger.debug("Calculating complexed ASAs for interface {} between molecules {} and {}", interf.getId(), molecId1, molecId2);

			interf.setAsas(chainAsas.get(molecId1), chainAsas.get(molecId2), nSpherePoints, nThreads, cofactorSizeToUse);

		}
		end = System.currentTimeMillis();

		logger.debug("Calculated complexes ASA for {} pairwise complexes. Time: {} s", redundancyReducedList.size(), ((end-start)/1000.0));

		// now let's populate the interface area value for the NCS-redundant ones from the reference interface (first one in list)
		if (clustersNcs!=null) {
			if (chainOrigNamesMap==null)  {
				logger.warn("No chainOrigNamesMap is set. Considering NCS interfaces in same order as reference. This is likely a bug.");
			}
			for (StructureInterfaceCluster ncsCluster : clustersNcs) {
				StructureInterface refInterf = ncsCluster.getMembers().get(0);
				String refMolecId1 = refInterf.getMoleculeIds().getFirst();
				for (int i=1;i<ncsCluster.getMembers().size();i++) {
					StructureInterface member = ncsCluster.getMembers().get(i);
					member.setTotalArea(refInterf.getTotalArea());
					String molecId1 = member.getMoleculeIds().getFirst();
					if (areMolecIdsSameOrder(refMolecId1, molecId1)) {
						// we add the reference interface GroupAsas as the GroupAsas for all other members, like that
						// ResidueNumbers won't match in their chain ids, but otherwise all info is there without using a lot of memory
						member.setFirstGroupAsas(refInterf.getFirstGroupAsas());
						member.setSecondGroupAsas(refInterf.getSecondGroupAsas());
					} else {
						member.setFirstGroupAsas(refInterf.getSecondGroupAsas());
						member.setSecondGroupAsas(refInterf.getFirstGroupAsas());
					}
				}
			}
		}

		// finally we sort based on the ChainInterface.comparable() (based in interfaceArea)
		sort();
	}

	private boolean areMolecIdsSameOrder(String refMolecId, String molecId) {

		if (chainOrigNamesMap==null) {
			// we've already warned above
			return true;
		}

		String refMolecIdOrig = chainOrigNamesMap.get(refMolecId);
		String molecIdOrig = chainOrigNamesMap.get(molecId);

		return (refMolecIdOrig.equals(molecIdOrig));
	}

	/**
	 * Sorts the interface list and reassigns ids based on new sorting
	 *
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
	 * Get the interface clusters for this StructureInterfaceList grouped by NCS-equivalence.
	 * This means that for any two interfaces in the same cluster:
	 * 1. The chains forming the first interface are NCS-copies of the chains forming the second interface, in any order.
	 * 2. Relative orientation of the chains is preserved, i.e. the contacts are identical.
	 * @return list of {@link StructureInterfaceCluster} objects.
	 * @since 5.0.0
	 */
	public List<StructureInterfaceCluster> getClustersNcs() {
		return clustersNcs;
	}

	/**
	 * Add an interface to the list, possibly defining it as NCS-equivalent to an interface already in the list.
	 * Used to build up the NCS clustering.
	 * @param interfaceNew
	 *          an interface to be added to the list.
	 * @param interfaceRef
	 *          interfaceNew will be added to the cluster which contains interfaceRef.
	 *          If interfaceRef is null, new cluster will be created for interfaceNew.
	 * @since 5.0.0
	 */
	public void addNcsEquivalent(StructureInterface interfaceNew, StructureInterface interfaceRef) {

		this.add(interfaceNew);
		if (clustersNcs == null) {
			clustersNcs = new ArrayList<>();
		}

		if (interfaceRef == null) {
			StructureInterfaceCluster newCluster = new StructureInterfaceCluster();
			newCluster.addMember(interfaceNew);
			clustersNcs.add(newCluster);
			return;
		}

		Optional<StructureInterfaceCluster> clusterRef =
				clustersNcs.stream().
					filter(r->r.getMembers().stream().
						anyMatch(c -> c.equals(interfaceRef))).
						findFirst();

		if (clusterRef.isPresent()) {
			clusterRef.get().addMember(interfaceNew);
			return;
		}

		logger.warn("The specified reference interface, if not null, should have been added to this set previously. " +
				"Creating new cluster and adding both interfaces. This is likely a bug.");
		this.add(interfaceRef);
		StructureInterfaceCluster newCluster = new StructureInterfaceCluster();
		newCluster.addMember(interfaceRef);
		newCluster.addMember(interfaceNew);
		clustersNcs.add(newCluster);
	}

	/**
	 * Sets a map with mapping from NCS chain names to original chain names.
	 * Necessary when {@link #addNcsEquivalent(StructureInterface, StructureInterface)} is used and NCS equivalent
	 * interfaces exist in this list and their names need mapping when setting ASAs.
	 * @param chainOrigNamesMap a map of NCS chain name to original chain name
	 */
	public void setChainOrigNamesMap(Map<String, String> chainOrigNamesMap) {
		this.chainOrigNamesMap = chainOrigNamesMap;
	}

	/**
	 * Removes from this interface list all interfaces with areas
	 * below the default cutoff area.
     * Note that this must be called after {@link #calcAsas(int, int, int)}, otherwise all areas would
     * be 0 and thus all removed.
	 * @see #DEFAULT_MINIMUM_INTERFACE_AREA
	 */
	public void removeInterfacesBelowArea() {
		removeInterfacesBelowArea(DEFAULT_MINIMUM_INTERFACE_AREA);
	}

	/**
	 * Removes from this interface list all interfaces with areas
	 * below the given cutoff area.
     * Note that this must be called after {@link #calcAsas(int, int, int)}, otherwise all areas would
     * be 0 and thus all removed.
	 * @param area the minimum interface buried surface area to keep. Interfaces below this value will be removed.
	 */
	public void removeInterfacesBelowArea(double area) {

		list.removeIf(interf -> interf.getTotalArea() < area);

	    if (clustersNcs != null) {
	    	clustersNcs.removeIf(ncsCluster -> ncsCluster.getMembers().get(0).getTotalArea() < area);
        }
	}

	/**
	 * Calculate the interface clusters for this StructureInterfaceList
	 * using a contact overlap score to measure the similarity of interfaces.
	 * Subsequent calls will use the cached value without recomputing the clusters.
	 * The contact overlap score cutoff to consider a pair in the same cluster is
	 * the value {@link #DEFAULT_CONTACT_OVERLAP_SCORE_CLUSTER_CUTOFF}
	 * @return
	 */
	public List<StructureInterfaceCluster> getClusters() {
		return getClusters(DEFAULT_CONTACT_OVERLAP_SCORE_CLUSTER_CUTOFF);
	}

	/**
	 * Calculate the interface clusters for this StructureInterfaceList
	 * using a contact overlap score to measure the similarity of interfaces.
	 * Subsequent calls will use the cached value without recomputing the clusters.
	 * The clusters will be assigned ids by sorting descending by {@link StructureInterfaceCluster#getTotalArea()}
	 * @param contactOverlapScoreClusterCutoff the contact overlap score above which a pair will be
	 * clustered
	 * @return
	 */
	public List<StructureInterfaceCluster> getClusters(double contactOverlapScoreClusterCutoff) {
		if (clusters!=null) {
			return clusters;
		}

		clusters = new ArrayList<>();

		// nothing to do if we have no interfaces
		if (list.size()==0) return clusters;

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

		Map<Integer, Set<Integer>> clusteredIndices = slc.getClusters(contactOverlapScoreClusterCutoff);
		for (int clusterIdx:clusteredIndices.keySet()) {
			List<StructureInterface> members = new ArrayList<>();
			for (int idx:clusteredIndices.get(clusterIdx)) {
				members.add(list.get(idx));
			}
			StructureInterfaceCluster cluster = new StructureInterfaceCluster();
			cluster.setMembers(members);
			double averageScore = 0.0;
			int countPairs = 0;
			for (int i=0;i<members.size();i++) {
                int iIdx = list.indexOf(members.get(i));
				for (int j=i+1;j<members.size();j++) {
					averageScore += matrix[iIdx][list.indexOf(members.get(j))];
					countPairs++;
				}
			}
			if (countPairs>0) {
				averageScore = averageScore/countPairs;
			} else {
				// if only one interface in cluster we set the score to the maximum
				averageScore = 1.0;
			}
			cluster.setAverageScore(averageScore);
			clusters.add(cluster);
		}

		// finally we have to set the back-references in each StructureInterface
		for (StructureInterfaceCluster cluster:clusters) {
			for (StructureInterface interf:cluster.getMembers()) {
				interf.setCluster(cluster);
			}
		}

		// now we sort by areas (descending) and assign ids based on that sorting
		clusters.sort((o1, o2) -> {
			return Double.compare(o2.getTotalArea(), o1.getTotalArea()); //note we invert so that sorting is descending
		});
		int id = 1;
		for (StructureInterfaceCluster cluster:clusters) {
			cluster.setId(id);
			id++;
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

	/**
	 * Calculates the interfaces for a structure using default parameters
	 * @param struc
	 * @return
	 */
	public static StructureInterfaceList calculateInterfaces(Structure struc) {
		CrystalBuilder builder = new CrystalBuilder(struc);
		StructureInterfaceList interfaces = builder.getUniqueInterfaces();
		logger.debug("Calculating ASA for "+interfaces.size()+" potential interfaces");
		interfaces.calcAsas(StructureInterfaceList.DEFAULT_ASA_SPHERE_POINTS, //fewer for performance
				Runtime.getRuntime().availableProcessors(),
				StructureInterfaceList.DEFAULT_MIN_COFACTOR_SIZE);
		interfaces.removeInterfacesBelowArea();
		interfaces.getClusters();
		logger.debug("Found "+interfaces.size()+" interfaces");
		return interfaces;
	}

}
