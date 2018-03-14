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
package org.biojava.nbio.structure.io.mmtf;

import java.util.List;
import java.util.Map;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;

/**
 * Class to store the summary data for a given structure.
 * @author Anthony Bradley
 *
 */
public class MmtfSummaryDataBean {

	private Map<String, Integer> chainIdToIndexMap;
	private List<Chain> allChains;
	private List<Atom> allAtoms;
	private int numBonds;
	
	/**
	 * @return the list of chains (in all models) in the structure
	 */
	public List<Chain> getAllChains() {
		return allChains;
	}
	/**
	 * @param allChains the list of chains (in all models) in the structure
	 */
	public void setAllChains(List<Chain> allChains) {
		this.allChains = allChains;
	}
	/**
	 * @return the list of atoms (in all models) in the structure
	 */
	public List<Atom> getAllAtoms() {
		return allAtoms;
	}
	/**
	 * @param allAtoms the list of atoms (in all models) in the structure
	 */
	public void setAllAtoms(List<Atom> allAtoms) {
		this.allAtoms = allAtoms;
	}
	/**
	 * @return the number of covalent bonds in the structure
	 */
	public int getNumBonds() {
		return numBonds;
	}
	/**
	 * @param numBonds the number of covalent bonds in the structure
	 */
	public void setNumBonds(int numBonds) {
		this.numBonds = numBonds;
	}
	/**
	 * @return the map of chain ids (strings asymId) to the index of that chain in the allChains list. 
	 * This only applies for the first model in the structure.
	 */
	public Map<String, Integer> getChainIdToIndexMap() {
		return chainIdToIndexMap;
	}
	/**
	 * @param chainIdToIndexMap the map of chain ids (strings asymId) to the index of that chain in the allChains list. 
	 * This only applies for the first model in the structure.
	 */
	public void setChainIdToIndexMap(Map<String, Integer> chainIdToIndexMap) {
		this.chainIdToIndexMap = chainIdToIndexMap;
	}
}
