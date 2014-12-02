package org.biojava.bio.structure.io;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeMap;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Entity;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class EntityFinder {

	private Structure s;
	
	private static final Logger logger = LoggerFactory.getLogger(EntityFinder.class);
	
	/**
	 * Identity value for 2 chains to be considered part of same entity
	 */
	public static final double IDENTITY_THRESHOLD = 0.99999;
	
	/**
	 * Gap coverage value (num gaps over length of sequence) for each chain of the match: 
	 * 2 chains with more gap coverage than this value will not be considered part of the same entity
	 */
	public static final double GAP_COVERAGE_THRESHOLD = 0.3;
	
	public EntityFinder(Structure s) {
		this.s = s;
	}
	
	public TreeMap<String,Entity> findEntities() {		

		TreeMap<String, Entity> chainIds2entities = new TreeMap<String,Entity>();
		List<String> chainIds = new ArrayList<String>();
		
		for (Chain c:s.getChains()) {
			chainIds.add(c.getChainID());			
		}
		
		Collections.sort(chainIds); //making sure they are in alphabetical order

		outer:
			for (int i=0;i<chainIds.size();i++) {
				for (int j=i+1;j<chainIds.size();j++) {
					Chain c1 = null;
					Chain c2 = null;
					try {
						c1 = s.getChainByPDB(chainIds.get(i));
						c2 = s.getChainByPDB(chainIds.get(j));
					} catch (StructureException e) {
						logger.error("Unexpected exception!",e);
					}

					if (areResNumbersAligned(c1, c2)) {

						if (	!chainIds2entities.containsKey(c1.getChainID()) &&
								!chainIds2entities.containsKey(c2.getChainID())) {

							logger.debug("Creating entity with chains {},{}",c1.getChainID(),c2.getChainID());

							Entity ent = new Entity();
							ent.addMember(c1);
							ent.setRepresentative(c1); // this will be the first alphabetically since they are sorted
							ent.addMember(c2);
							chainIds2entities.put(c1.getChainID(), ent);
							chainIds2entities.put(c2.getChainID(), ent);

						} else {
							Entity ent = chainIds2entities.get(c1.getChainID());

							if (ent==null) {
								logger.debug("Adding chain {} to entity {}",c1.getChainID(),c2.getChainID());
								ent = chainIds2entities.get(c2.getChainID());
								ent.addMember(c1);
								chainIds2entities.put(c1.getChainID(), ent);

							} else {
								logger.debug("Adding chain {} to entity {}",c2.getChainID(),c1.getChainID());
								ent.addMember(c2);
								chainIds2entities.put(c2.getChainID(), ent);								
							}
						}
					} 


					if (chainIds2entities.size()==chainIds.size()) // we've got all chains in entities
						break outer;
				}
			}
		// anything not in an entity will be its own entity
		for (Chain c:s.getChains()) {
			if (!chainIds2entities.containsKey(c.getChainID())) {
				logger.debug("Creating a 1-member entity for chain {}",c.getChainID());
				chainIds2entities.put(c.getChainID(),getTrivialEntity(c));
			}
		}




		return chainIds2entities;
	}
	
	private static Entity getTrivialEntity(Chain c) {
		List<Chain> members = new ArrayList<Chain>();
		members.add(c);
		Entity ent = new Entity(c,members);
		return ent;
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
		
		if ((double)countGaps/(double)c1AtomGroups.size() > 0.50) {
			logger.info("More than half the residues ({} out of {}) are gaps in chain {} when compared to chain {}. Will not group these 2 chains in an entity", 
					countGaps, c1AtomGroups.size(), c2.getChainID(), c1.getChainID());
			return false;
		}

		return true;
	}
	
	/**
	 * Tell whether given chain is a protein chain (true) or a nucleotide chain (false)
	 * @param c
	 * @return true if protein, false if nucleotide
	 */
	private static boolean isProtein(Chain c) {
		int sizeAminos = c.getAtomGroups(GroupType.AMINOACID).size();
		int sizeNucleotides = c.getAtomGroups(GroupType.NUCLEOTIDE).size();
		List<Group> hetAtoms = c.getAtomGroups(GroupType.HETATM);
		int sizeHetatoms = hetAtoms.size();
		int sizeWaters = 0;
		for (Group g:hetAtoms) {
			if (g.getPDBName().equals("HOH")) sizeWaters++;
		}
		int fullSize = sizeAminos + sizeNucleotides + sizeHetatoms - sizeWaters;
		
		if ((double)sizeAminos/(double)fullSize>0.95) return true;
		
		if ((double)sizeNucleotides/(double)fullSize>0.95) return false;
		
		// finally if neither condition works, we try based on majority
		if (sizeAminos>sizeNucleotides) {
			logger.debug("Ratio of residues to total for chain {} is below 95%. Assuming it is a protein chain. Counts: # aa residues: {}, # nuc residues: {}, # het residues: {}, # waters: {}, ratio aa/total: {}, ratio nuc/total: {}",
				c.getChainID(), sizeAminos, sizeNucleotides, sizeHetatoms, sizeWaters,
				(double)sizeAminos/(double)fullSize,(double)sizeNucleotides/(double)fullSize) ;
			
			return true;
			
		} else {
			logger.debug("Ratio of residues to total for chain {} is below 95%. Assuming it is a nucleotide chain. Counts: # aa residues: {}, # nuc residues: {}, # het residues: {}, # waters: {}, ratio aa/total: {}, ratio nuc/total: {}",
					c.getChainID(), sizeAminos, sizeNucleotides, sizeHetatoms, sizeWaters,
					(double)sizeAminos/(double)fullSize,(double)sizeNucleotides/(double)fullSize) ;

			return false;
		}
		
	}
	

}
 