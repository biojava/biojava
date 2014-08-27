package org.biojava.bio.structure.contact;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.asa.AsaCalculator;


/**
 * A list of interfaces between 2 molecules (2 sets of atoms) 
 * 
 * @author duarte_j
 *
 */
public class StructureInterfaceList implements Serializable, Iterable<StructureInterface> {

	private static final long serialVersionUID = 1L;
	
	private List<StructureInterface> list;
	
	private boolean debug;
	

	public StructureInterfaceList() {
		this.list = new ArrayList<StructureInterface>();
		this.debug = false;
	}
	
	public void add(StructureInterface interf) {
		this.list.add(interf);
	}
	
	public int size() {
		return this.list.size();
	}
	
	public void setDebug(boolean debug) {
		// TODO do the debugging output through logging once we have slf4j in BJ
		this.debug = debug;
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
		
		if (debug) 
			System.out.println("Calculated uncomplexed ASA for "+uniqAsaChains.size()+" orientation-unique chains. "
					+ "Time: "+((end-start)/1000.0)+" s");
		
		start = System.currentTimeMillis();
		
		// now we calculate the ASAs for each of the complexes 
		for (StructureInterface interf:list) {
			
			String molecId1 = interf.getMoleculeIds().getFirst()+interf.getTransforms().getFirst().getTransformId();
			String molecId2 = interf.getMoleculeIds().getSecond()+interf.getTransforms().getSecond().getTransformId();
			
			interf.setAsas(chainAsas.get(molecId1), chainAsas.get(molecId2), nSpherePoints, nThreads, cofactorSizeToUse);
			
		}
		end = System.currentTimeMillis();
		
		if (debug) 
			System.out.println("Calculated complexes ASA for "+list.size()+" pairwise complexes. "
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

	@Override
	public Iterator<StructureInterface> iterator() {
		return list.iterator();
	}
}
