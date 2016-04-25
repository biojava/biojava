package org.biojava.nbio.structure.io.mmtf;

import java.util.List;
import java.util.Map;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;

/**
 * Class to store the summary data for a given structure
 * @author Anthony Bradley
 *
 */
public class MmtfSummaryDataBean {

	private Map<String, Integer> chainIdToIndexMap;
	private List<Chain> allChains;
	private List<Atom> allAtoms;
	private int numBonds;
	
	/**
	 * @return the allChains
	 */
	public List<Chain> getAllChains() {
		return allChains;
	}
	/**
	 * @param allChains the allChains to set
	 */
	public void setAllChains(List<Chain> allChains) {
		this.allChains = allChains;
	}
	/**
	 * @return the allAtoms
	 */
	public List<Atom> getAllAtoms() {
		return allAtoms;
	}
	/**
	 * @param allAtoms the allAtoms to set
	 */
	public void setAllAtoms(List<Atom> allAtoms) {
		this.allAtoms = allAtoms;
	}
	/**
	 * @return the numBonds
	 */
	public int getNumBonds() {
		return numBonds;
	}
	/**
	 * @param numBonds the numBonds to set
	 */
	public void setNumBonds(int numBonds) {
		this.numBonds = numBonds;
	}
	/**
	 * @return the chainIdToIndexMap
	 */
	public Map<String, Integer> getChainIdToIndexMap() {
		return chainIdToIndexMap;
	}
	/**
	 * @param chainIdToIndexMap the chainIdToIndexMap to set
	 */
	public void setChainIdToIndexMap(Map<String, Integer> chainIdToIndexMap) {
		this.chainIdToIndexMap = chainIdToIndexMap;
	}
}
