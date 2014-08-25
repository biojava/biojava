package org.biojava.bio.structure;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.biojava.bio.structure.asa.AsaCalculator;



public class ChainInterfaceList implements Serializable, Iterable<ChainInterface> {

	private static final long serialVersionUID = 1L;
	
	private List<ChainInterface> list;
	
	private boolean debug;
	

	public ChainInterfaceList() {
		this.list = new ArrayList<ChainInterface>();
		this.debug = false;
	}
	
	public void add(ChainInterface interf) {
		this.list.add(interf);
	}
	
	public int size() {
		return this.list.size();
	}
	
	public void setDebug(boolean debug) {
		this.debug = debug;
	}
	
	/**
	 * Gets the interface corresponding to given id.
	 * The ids go from 1 to n
	 * If {@link #sort()} was called then the order is descendent by area.
	 * @param id
	 * @return
	 */
	public ChainInterface get(int id) {
		return list.get(id-1);
	}
	
	/**
	 * Calculates ASAs for all interfaces in list, both for the uncomplexed 
	 * chains and for the complex of the two chains together. 
	 * Also sorts the interfaces based on calculated BSA areas (descending) 
	 * @param nSpherePoints
	 * @param nThreads
	 * @param hetAtoms whether to use any HET atoms that are part of the chain
	 * @param cofactorSizeToUse the minimum size of cofactor molecule (non-chain HET atoms) that will be used
	 */
	public void calcAsas(int nSpherePoints, int nThreads, boolean hetAtoms, int cofactorSizeToUse) {
				
		// asa/bsa calculation 
		// NOTE in principle it is more efficient to calculate asas only once per unique chain
		// BUT! surprisingly the rolling ball algorithm gives slightly different values for same molecule in different 
		// rotations (due to sampling depending on orientation of axes grid). 
		// Both NACCESS and our own implementation behave like that.
		// That's why we calculate ASAs for each rotation-unique molecule, otherwise 
		// we get discrepancies (not very big but annoying) which lead to things like negative (small) bsa values
		
		
		Map<String, Atom[]> uniqAsaChains = new TreeMap<String, Atom[]>();
		Map<String, double[]> chainAsas = new TreeMap<String, double[]>();
		
		// first we gather rotation-unique chains (in terms of AU id and transform id)
		for (ChainInterface interf:list) {
			String molecId1 = interf.getChains().getFirst().getChainID()+interf.getTransforms().getFirst().getTransformId();
			String molecId2 = interf.getChains().getSecond().getChainID()+interf.getTransforms().getSecond().getTransformId();
			
			uniqAsaChains.put(molecId1, interf.getFirstAtoms(hetAtoms, cofactorSizeToUse)); 
			uniqAsaChains.put(molecId2, interf.getSecondAtoms(hetAtoms, cofactorSizeToUse));
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
		for (ChainInterface interf:list) {
			
			String molecId1 = interf.getChains().getFirst().getChainID()+interf.getTransforms().getFirst().getTransformId();
			String molecId2 = interf.getChains().getSecond().getChainID()+interf.getTransforms().getSecond().getTransformId();
			
			interf.setAsas(chainAsas.get(molecId1), chainAsas.get(molecId2), nSpherePoints, nThreads, hetAtoms, cofactorSizeToUse);
			
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
		for (ChainInterface interf:list) {
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
		Iterator<ChainInterface> it = iterator();
		while (it.hasNext()) {
			ChainInterface interf = it.next();
			if (interf.getInterfaceArea()<area) {
				it.remove();
			}
		}
	}

	@Override
	public Iterator<ChainInterface> iterator() {
		return list.iterator();
	}
}
