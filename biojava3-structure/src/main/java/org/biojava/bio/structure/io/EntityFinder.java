package org.biojava.bio.structure.io;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Entity;
import org.biojava.bio.structure.Structure;
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

		try {
			outer:
			for (int i=0;i<s.getChains().size();i++) {
				for (int j=i+1;j<s.getChains().size();j++) {
					Chain c1 = s.getChain(i);
					Chain c2 = s.getChain(j);
					ProteinSequence s1 = new ProteinSequence(c1.getAtomSequence());
					ProteinSequence s2 = new ProteinSequence(c2.getAtomSequence());

					SequencePair<ProteinSequence, AminoAcidCompound> pair = align(s1,s2);
					
					double identity = (double)pair.getNumIdenticals()/(double)pair.getLength();
					logger.info("Identity for chain pair {},{}: {}", c1.getChainID(), c2.getChainID(), identity);
					
					if (identity > 0.9999) {
						if (	!chainIds2entities.containsKey(c1.getChainID()) &&
								!chainIds2entities.containsKey(c2.getChainID())) {
							logger.info("Creating entity with chains {},{}",c1.getChainID(),c2.getChainID());
							Entity ent = new Entity();
							ent.addMember(c1);
							ent.setRepresentative(c1); // TODO this should be the first alphabetically
							ent.addMember(c2);
							chainIds2entities.put(c1.getChainID(), ent);
							chainIds2entities.put(c2.getChainID(), ent);
						} else {
							Entity ent = chainIds2entities.get(c1.getChainID());
							
							if (ent==null) {
								logger.info("Adding chain {} to entity {}",c1.getChainID(),c2.getChainID());
								ent = chainIds2entities.get(c2.getChainID());
								ent.addMember(c1);
								chainIds2entities.put(c1.getChainID(), ent);
							} else {
								logger.info("Adding chain {} to entity {}",c2.getChainID(),c1.getChainID());
								ent.addMember(c2);
								chainIds2entities.put(c2.getChainID(), ent);
							}
						}
					}
					if (chainIds2entities.size()==s.getChains().size()) // we've got all chains in entities
						break outer;
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
			List<Chain> members = new ArrayList<Chain>();
			members.add(c);
			Entity ent = new Entity(c,members);
			list.put(c.getChainID(),ent);
		}
		
		return list;
	}
}
