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
package org.biojava.bio.structure.io;

import org.biojava.bio.structure.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

/**
 * Heuristically finding of Compounds (called Entities in mmCIF dictionary)
 * in a given Structure. Compounds are the groups of sequence identical NCS-related polymer chains
 * in the Structure.
 * 
 * @author duarte_j
 */
public class CompoundFinder {

	private Structure s;
	
	private static final Logger logger = LoggerFactory.getLogger(CompoundFinder.class);
	
	/**
	 * Above this ratio of mismatching residue types for same residue numbers we 
	 * consider the 2 chains to be a mismatch and don't put them together in the same entity
	 */
	public static final double RATIO_GAPS_FOR_MISMATCH = 0.50;
	
	
	public CompoundFinder(Structure s) {
		this.s = s;
	}
	
	/**
	 * Utility method that employs some heuristics to find the Compounds
	 * for this Structure in case the information is missing in PDB/mmCIF file
	 * @return
	 */
	public List<Compound> findCompounds() {
		
		TreeMap<String,Compound> chainIds2entities = findCompoundsFromAtomSequences();

		return findUniqueCompounds(chainIds2entities);
	}
	
	/**
	 * Utility method to obtain a list of unique entities from the chainIds2entities map
	 * @return
	 */
	private static List<Compound> findUniqueCompounds(TreeMap<String,Compound> chainIds2entities) {
		
		List<Compound> list = new ArrayList<Compound>();
		
		for (Compound cluster:chainIds2entities.values()) {
			boolean present = false;
			for (Compound cl:list) {
				if (cl==cluster) {
					present = true;
					break;
				}
			}
			if (!present) list.add(cluster);
		}
		return list;
	} 
	
	private TreeMap<String,Compound> findCompoundsFromAtomSequences() {		

		TreeMap<String, Compound> chainIds2compounds = new TreeMap<String,Compound>();

		int molId = 1;
		outer:
			for (int i=0;i<s.getChains().size();i++) {
				for (int j=i+1;j<s.getChains().size();j++) {

					Chain c1 = s.getChain(i);
					Chain c2 = s.getChain(j);
					
					if (areResNumbersAligned(c1, c2)) {

						if (	!chainIds2compounds.containsKey(c1.getChainID()) &&
								!chainIds2compounds.containsKey(c2.getChainID())) {

							logger.debug("Creating Compound with chains {},{}",c1.getChainID(),c2.getChainID());

							Compound ent = new Compound();
							ent.addChain(c1);
							ent.addChain(c2);
							ent.setMolId(molId++);
							chainIds2compounds.put(c1.getChainID(), ent);
							chainIds2compounds.put(c2.getChainID(), ent);

						} else {
							Compound ent = chainIds2compounds.get(c1.getChainID());

							if (ent==null) {
								logger.debug("Adding chain {} to Compound {}",c1.getChainID(),c2.getChainID());
								ent = chainIds2compounds.get(c2.getChainID());
								ent.addChain(c1);
								chainIds2compounds.put(c1.getChainID(), ent);

							} else {
								logger.debug("Adding chain {} to Compound {}",c2.getChainID(),c1.getChainID());
								ent.addChain(c2);
								chainIds2compounds.put(c2.getChainID(), ent);								
							}
						}
					} 


					if (chainIds2compounds.size()==s.getChains().size()) // we've got all chains in Compounds
						break outer;
				}
			}

		// anything not in a Compound will be its own Compound
		for (Chain c:s.getChains()) {
			if (!chainIds2compounds.containsKey(c.getChainID())) {
				logger.debug("Creating a 1-member Compound for chain {}",c.getChainID());
				
				Compound ent = new Compound();
				ent.addChain(c);
				ent.setMolId(molId++);
				
				chainIds2compounds.put(c.getChainID(),ent);
			}
		}


		return chainIds2compounds;
	}
	
	
	private static boolean areResNumbersAligned(Chain c1, Chain c2) {

		boolean isC1prot = StructureTools.isProtein(c1);
		boolean isC2prot = StructureTools.isProtein(c2);
		
		// different kind of chain: we won't try to align them
		if (isC1prot != isC2prot ) return false;
		
		List<Group> c1AtomGroups = null;
		if (isC1prot) {
			c1AtomGroups = c1.getAtomGroups(GroupType.AMINOACID);
		}
		else {
			c1AtomGroups = c1.getAtomGroups(GroupType.NUCLEOTIDE);
		}
		
		int countGaps = 0;
		
		for (Group g1:c1AtomGroups) {
			try {
				Group g2 = c2.getGroupByPDB(g1.getResidueNumber());
				if (!g2.getPDBName().equals(g1.getPDBName())) {
					logger.debug("Mismatch of residues between chains {},{} for residue number {}: {} {}",
							c1.getChainID(),c2.getChainID(),g1.getResidueNumber(), g1.getPDBName(), g2.getPDBName());
					return false;
				}
			} catch (StructureException e) {
				// the group doesn't exist (no density) in the chain, go on
				countGaps++;
				continue;
			}
		}
		
		if ((double)countGaps/(double)c1AtomGroups.size() > RATIO_GAPS_FOR_MISMATCH) {
			logger.info("More than {} of the residues ({} out of {}) are gaps in chain {} when compared to chain {}. Will not group these 2 chains in an entity", 
					RATIO_GAPS_FOR_MISMATCH, countGaps, c1AtomGroups.size(), c2.getChainID(), c1.getChainID());
			return false;
		}

		return true;
	}
	


}
 