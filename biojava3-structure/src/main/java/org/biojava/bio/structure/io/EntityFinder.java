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
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.PairwiseSequenceAligner;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.exceptions.CompoundNotFoundException;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class EntityFinder {

	private Structure s;
	
	private static final Logger logger = LoggerFactory.getLogger(EntityFinder.class);
	
	public EntityFinder(Structure s) {
		this.s = s;
	}
	
	public TreeMap<String,Entity> findEntities() {
		
		
		TreeMap<String, Entity> chainIds2entities = new TreeMap<String,Entity>();
		List<String> chainIds = new ArrayList<String>();
		for (Chain c:s.getChains()) {
			if ( c.getAtomGroups(GroupType.AMINOACID).size() > 1 &&
					c.getAtomGroups(GroupType.NUCLEOTIDE).size() < 1) {	
				chainIds.add(c.getChainID());
			}
		}
		Collections.sort(chainIds); //making sure they are in alphabetical order
		try {
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
					ProteinSequence s1 = new ProteinSequence(c1.getAtomSequence());
					ProteinSequence s2 = new ProteinSequence(c2.getAtomSequence());

					SequencePair<ProteinSequence, AminoAcidCompound> pair = align(s1,s2);
					
					int nonGaps = pair.getLength()-getNumGaps(pair); 
					double identity = (double)pair.getNumIdenticals()/(double)nonGaps;
					logger.debug("Identity for chain pair {},{}: {}", c1.getChainID(), c2.getChainID(), identity);
					
					if (identity > 0.99999) {
						if (	!chainIds2entities.containsKey(c1.getChainID()) &&
								!chainIds2entities.containsKey(c2.getChainID())) {
							logger.debug("Creating entity with chains {},{}",c1.getChainID(),c2.getChainID());
							if (areResNumbersAligned(c1, c2)) { 
								Entity ent = new Entity();
								ent.addMember(c1);
								ent.setRepresentative(c1); // this will be the first alphabetically since they are sorted
								ent.addMember(c2);
								chainIds2entities.put(c1.getChainID(), ent);
								chainIds2entities.put(c2.getChainID(), ent);
							} else {
								logger.warn("Chains {},{} with 100% identity have residue numbers misaligned, can't include them in the same entity");
							}
						} else {
							Entity ent = chainIds2entities.get(c1.getChainID());
							
							if (ent==null) {
								logger.debug("Adding chain {} to entity {}",c1.getChainID(),c2.getChainID());
								ent = chainIds2entities.get(c2.getChainID());
								if (areResNumbersAligned(c1, c2)) { 
									ent.addMember(c1);
									chainIds2entities.put(c1.getChainID(), ent);
								} else {
									logger.warn("Chains {},{} with 100% identity have residue numbers misaligned, can't include them in the same entity");
								}
							} else {
								logger.debug("Adding chain {} to entity {}",c2.getChainID(),c1.getChainID());
								if (areResNumbersAligned(c1, c2)) {
									ent.addMember(c2);
									chainIds2entities.put(c2.getChainID(), ent);
								} else {
									logger.warn("Chains {},{} with 100% identity have residue numbers misaligned, can't include them in the same entity");
								}
							}
						}
					} else if (identity>0.95) {
						logger.info("Identity for chains {},{} above 0.95. Mismatch in {} out of {} non-gap-aligned residues (identity {})",
								c1.getChainID(),c2.getChainID(),nonGaps-pair.getNumIdenticals(),nonGaps,identity);
						logger.info("\n"+pair.toString(100)); 
					} else if (identity>1) {
						logger.warn("Identity for chains {},{} above 1. {} identicals out of {} non-gap-aligned residues (identity {})",
								c1.getChainID(),c2.getChainID(),pair.getNumIdenticals(),nonGaps,identity);
						logger.warn("\n"+pair.toString(100));
					}
					
					if (chainIds2entities.size()==chainIds.size()) // we've got all chains in entities
						break outer;
				}
			}
			// anything not in an entity will be its own entity
			for (Chain c:s.getChains()) {
				if (!chainIds2entities.containsKey(c.getChainID())) {
					logger.debug("Creating a 1-member entity for chain "+c.getChainID());
					chainIds2entities.put(c.getChainID(),getTrivialEntity(c));
				}
			}
		
		} catch (CompoundNotFoundException e) {
			logger.warn("Some unknown compounds in protein sequences, will use an entity per chain. Error: "+e.getMessage());
			return getTrivialEntities();
		}

		
		return chainIds2entities;
	}
	
	private SequencePair<ProteinSequence, AminoAcidCompound> align(ProteinSequence s1, ProteinSequence s2) {
		SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum65();

		GapPenalty penalty = new SimpleGapPenalty();

		short gop = 8;
		short extend = 1;
		penalty.setOpenPenalty(gop);
		penalty.setExtensionPenalty(extend);


		PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> smithWaterman =
				Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.LOCAL, penalty, matrix);

		return smithWaterman.getPair();
	}
	
	private TreeMap<String,Entity> getTrivialEntities() {
		TreeMap<String,Entity> list = new TreeMap<String,Entity>();
		
		for (Chain c:s.getChains()) {
			list.put(c.getChainID(),getTrivialEntity(c)); 
		}
		
		return list;
	}
	
	private static Entity getTrivialEntity(Chain c) {
		List<Chain> members = new ArrayList<Chain>();
		members.add(c);
		Entity ent = new Entity(c,members);
		return ent;
	}
	
	private static boolean areResNumbersAligned(Chain rep, Chain c) {

		for (Group repG:rep.getAtomGroups()) {
			try {
				Group g = c.getGroupByPDB(repG.getResidueNumber());
				if (!g.getPDBName().equals(repG.getPDBName())) {
					logger.info("Mismatch of residues between chains of same entity {},{} for residue number {}",
							rep.getChainID(),c.getChainID(),repG.getResidueNumber());
					return false;
				}
			} catch (StructureException e) {
				// the group doesn't exist (no density) in the chain, go on
				continue;
			}
		}

		return true;
	}
	
	private static int getNumGaps(SequencePair<ProteinSequence, AminoAcidCompound> pair) {
		int numGaps = 0;
		for (int alignmentIndex=1;alignmentIndex<=pair.getLength();alignmentIndex++) {
			if (pair.hasGap(alignmentIndex)) numGaps++;
		}
		return numGaps;
	}
}
 