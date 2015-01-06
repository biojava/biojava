package org.biojava.bio.structure.io;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Compound;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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
	 * Below this ratio of aminoacid/nucleotide residues to the sequence total,
	 * we use simple majority of aminoacid/nucleotide residues to decide the character 
	 * of the chain (protein/nucleotide) 
	 */
	public static final double RATIO_RESIDUES_TO_TOTAL = 0.95;
	
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
	
	@SuppressWarnings("unused")
	private TreeMap<String,Compound> findCompoundsFromSeqresSequences() {

		TreeMap<String,Compound> chainIds2entities = new TreeMap<String,Compound>();
				
		// map of sequences to list of chain identifiers
		Map<String, List<String>> uniqSequences = new HashMap<String, List<String>>();
		// finding the entities (groups of identical chains)
		for (Chain chain:s.getChains()) {
			
			// TODO chain.getSeqResSequence() below will get 'XXXXX' sequences for nucleotides, we need to change that to get the right sequence
			// e.g. with this procedure 1g1n results in 2 entities, where there are actually 3
			// in the meanwhile we warn about that:
			if (!CompoundFinder.isProtein(chain)) {
				logger.warn("Chain {} looks like a nucleotide chain, entity finding will not be accurate for  it",chain.getChainID());
			}

			String seq = chain.getSeqResSequence();
				
			if (uniqSequences.containsKey(seq)) {
				uniqSequences.get(seq).add(chain.getChainID());
			} else {
				List<String> list = new ArrayList<String>();
				list.add(chain.getChainID());
				uniqSequences.put(seq, list);
			}		

		}

		for (List<String> chainIds:uniqSequences.values()) {
			// sorting ids in alphabetic order
			Collections.sort(chainIds);
			List<Chain> chains = new ArrayList<Chain>();
			for (String chainId:chainIds) {
				// chains will be sorted in ids' alphabetic order
				try {
					chains.add(s.getChainByPDB(chainId));
				} catch (StructureException e) {
					// this basically can't happen, if it does it is some kind of bug
					logger.error("Unexpected exception!",e);
				}
			}
			// the representative will be the one with first chain id in alphabetic order 
			Compound entity = new Compound();
			entity.setChains(chains);
			for (Chain member:chains) {
				chainIds2entities.put(member.getChainID(), entity);
			}
		}
		
		return chainIds2entities;
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

		boolean isC1prot = isProtein(c1);
		boolean isC2prot = isProtein(c2);
		
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
	
	/**
	 * Tell whether given chain is a protein chain
	 * @param c
	 * @return true if protein, false if nucleotide or ligand
	 */
	public static boolean isProtein(Chain c) {
		return getPredominantGroupType(c) == GroupType.AMINOACID;
	}
	/**
	 * Tell whether given chain is DNA or RNA
	 * @param c
	 * @return true if nucleic acid, false if protein or ligand
	 */
	public static boolean isNucleicAcid(Chain c) {
		return getPredominantGroupType(c) == GroupType.NUCLEOTIDE;
	}
	/**
	 * Gets the predominant GroupType for a given Chain
	 * @param c
	 * @return
	 */
	public static GroupType getPredominantGroupType(Chain c) {
		int sizeAminos = c.getAtomGroups(GroupType.AMINOACID).size();
		int sizeNucleotides = c.getAtomGroups(GroupType.NUCLEOTIDE).size();
		List<Group> hetAtoms = c.getAtomGroups(GroupType.HETATM);
		int sizeHetatoms = hetAtoms.size();
		int sizeWaters = 0;
		for (Group g:hetAtoms) {
			if (g.isWater()) sizeWaters++;
		}
		
		int fullSize = sizeAminos + sizeNucleotides + sizeHetatoms - sizeWaters;
		
		if ((double)sizeAminos/(double)fullSize>RATIO_RESIDUES_TO_TOTAL) return GroupType.AMINOACID;
		
		if ((double)sizeNucleotides/(double)fullSize>RATIO_RESIDUES_TO_TOTAL) return GroupType.NUCLEOTIDE;
		
		if ((double)(sizeHetatoms-sizeWaters)/(double)fullSize > RATIO_RESIDUES_TO_TOTAL) return GroupType.HETATM;
		
		// finally if neither condition works, we try based on majority, but log it
		GroupType max;
		if(sizeNucleotides > sizeAminos) {
			if(sizeNucleotides > sizeHetatoms) {
				max = GroupType.NUCLEOTIDE;
			} else {
				max = GroupType.HETATM;
			}
		} else {
			if(sizeAminos > sizeHetatoms) {
				max = GroupType.AMINOACID;
			} else {
				max = GroupType.HETATM;
			}
		}
		logger.debug("Ratio of residues to total for chain {} is below {}. Assuming it is a {} chain. "
				+ "Counts: # aa residues: {}, # nuc residues: {}, # het residues: {}, # waters: {}, "
				+ "ratio aa/total: {}, ratio nuc/total: {}",
				c.getChainID(), RATIO_RESIDUES_TO_TOTAL, max,
				sizeAminos, sizeNucleotides, sizeHetatoms, sizeWaters,
				(double)sizeAminos/(double)fullSize,(double)sizeNucleotides/(double)fullSize) ;

		return max;
	}


}
 