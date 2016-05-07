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
package org.biojava.nbio.structure.io;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.structure.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Heuristical finding of Entities (called Compounds in legacy PDB format)
 * in a given Structure. Entities are the groups of sequence identical NCS-related polymer chains
 * in the Structure.
 *
 * This is related to {@link SeqRes2AtomAligner} but it is intended for raw PDB/mmCIF files where
 * possibly no SEQRES is given.
 *
 * @author Jose Duarte
 */
public class EntityFinder {

	private static final Logger logger = LoggerFactory.getLogger(EntityFinder.class);
	
	/**
	 * Above this ratio of mismatching residue types for same residue numbers we
	 * consider the 2 chains to have mismatching residue numbers and warn about it
	 */
	public static final double RATIO_GAPS_FOR_MISMATCH = 0.50;

	/**
	 * Identity value for 2 chains to be considered part of same entity
	 */
	public static final double IDENTITY_THRESHOLD = 0.99999;

	/**
	 * Gap coverage value (num gaps over length of sequence) for each chain of the match:
	 * 2 chains with more gap coverage than this value will not be considered part of the same entity
	 */
	public static final double GAP_COVERAGE_THRESHOLD = 0.3;

	//private Structure s;
	private List<List<Chain>> allModels;

	public EntityFinder(Structure s) {
		allModels = new ArrayList<List<Chain>>();
		for (int i=0;i<s.nrModels();i++) {
			this.allModels.add(s.getModel(i));
		}
	}
	
	public EntityFinder(List<List<Chain>> allModels) {
		this.allModels = allModels;
	}

	/**
	 * Utility method that employs some heuristics to find the {@link EntityInfo}s
	 * for this Structure in case the information is missing in PDB/mmCIF file
	 * @return
	 */
	public List<EntityInfo> findEntities() {

		TreeMap<String,EntityInfo> chainIds2entities = findCompoundsFromAlignment();

		List<EntityInfo> compounds = findUniqueCompounds(chainIds2entities);

		createPurelyNonPolyCompounds(compounds);

		return compounds;
	}

	/**
	 * Utility method to obtain a list of unique entities from the chainIds2entities map
	 * @return
	 */
	private static List<EntityInfo> findUniqueCompounds(TreeMap<String,EntityInfo> chainIds2entities) {

		List<EntityInfo> list = new ArrayList<EntityInfo>();

		for (EntityInfo cluster:chainIds2entities.values()) {
			boolean present = false;
			for (EntityInfo cl:list) {
				if (cl==cluster) {
					present = true;
					break;
				}
			}
			if (!present) list.add(cluster);
		}
		return list;
	}

	private void createPurelyNonPolyCompounds(List<EntityInfo> compounds) {

		List<Chain> pureNonPolymerChains = new ArrayList<Chain>();
		for (int i=0;i<allModels.get(0).size();i++) {
			if (StructureTools.isChainPureNonPolymer(allModels.get(0).get(i))) {
				pureNonPolymerChains.add(allModels.get(0).get(i));
			}
		}

		if (!pureNonPolymerChains.isEmpty()) {

			int maxMolId = 0; // if we have no compounds at all we just use 0, so that first assignment is 1
			if (!compounds.isEmpty()) {
				Collections.max(compounds, new Comparator<EntityInfo>() {
					@Override
					public int compare(EntityInfo o1, EntityInfo o2) {
						return new Integer(o1.getMolId()).compareTo(o2.getMolId());
					}
				}).getMolId();
			}

			int molId = maxMolId + 1;
			// for the purely non-polymeric chains we assign additional compounds
			for (Chain c: pureNonPolymerChains) {
				EntityInfo comp = new EntityInfo();
				comp.addChain(c);
				comp.setMolId(molId);
				if (StructureTools.isChainWaterOnly(c)) {
					comp.setType(EntityType.WATER);
				} else {
					comp.setType(EntityType.NONPOLYMER);
				}
				logger.warn("Chain {} is purely non-polymeric, will assign a new Compound (entity) to it (entity id {})", c.getId(), molId);
				molId++;

				compounds.add(comp);
			}

		}
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

		int countNonExisting = 0;

		for (Group g1:c1AtomGroups) {
			try {
				Group g2 = c2.getGroupByPDB(g1.getResidueNumber());
				if (!g2.getPDBName().equals(g1.getPDBName())) {
					logger.debug("Mismatch of residues between chains {},{} for residue number {}: {} {}",
							c1.getId(),c2.getId(),g1.getResidueNumber(), g1.getPDBName(), g2.getPDBName());
					return false;
				}
			} catch (StructureException e) {
				// the group doesn't exist (no density) in the chain, go on
				countNonExisting++;
				continue;
			}
		}

		if ((double)countNonExisting/(double)c1AtomGroups.size() > RATIO_GAPS_FOR_MISMATCH) {
			logger.debug("More than {} of the residues ({} out of {}) are not present in chain {} when comparing by residue numbers to chain {}.",
					RATIO_GAPS_FOR_MISMATCH, countNonExisting, c1AtomGroups.size(), c2.getId(), c1.getId());
			return false;
		}

		return true;
	}



	private TreeMap<String,EntityInfo> findCompoundsFromAlignment() {

		// first we determine which chains to consider: anything not looking
		// polymeric (protein or nucleotide chain) should be discarded
		Set<Integer> polyChainIndices = new TreeSet<Integer>();
		List<Chain> pureNonPolymerChains = new ArrayList<Chain>();
		for (int i=0;i<allModels.get(0).size();i++) {
			if (StructureTools.isChainPureNonPolymer(allModels.get(0).get(i))) {
				pureNonPolymerChains.add(allModels.get(0).get(i));
			} else {
				polyChainIndices.add(i);
			}
		}


		TreeMap<String, EntityInfo> chainIds2compounds = new TreeMap<String,EntityInfo>();

		int molId = 1;

		outer:
			for (int i:polyChainIndices) {
				for (int j:polyChainIndices) {

					if (j<=i) continue;

					Chain c1 = allModels.get(0).get(i);
					Chain c2 = allModels.get(0).get(j);

					Map<Integer,Integer> positionIndex1 = new HashMap<Integer, Integer>();
					Map<Integer,Integer> positionIndex2 = new HashMap<Integer, Integer>();
					// here we use false, which means that X will be used for unknown compounds
					String str1 = SeqRes2AtomAligner.getFullAtomSequence(c1.getAtomGroups(), positionIndex1, false);
					String str2 = SeqRes2AtomAligner.getFullAtomSequence(c2.getAtomGroups(), positionIndex2, false);

					int seq1Length = 0;
					int seq2Length = 0;

					SequencePair<?,?> pair = null;
					if (isProteinSequence(str1) && isProteinSequence(str2)) {
						ProteinSequence s1 = getProteinSequence(str1);
						ProteinSequence s2 = getProteinSequence(str2);
						seq1Length = s1.getLength();
						seq2Length = s2.getLength();

						pair = alignProtein(s1,s2);

					} else if (isDNASequence(str1) && isDNASequence(str2)) {
						DNASequence s1 = getDNASequence(str1);
						DNASequence s2 = getDNASequence(str2);
						seq1Length = s1.getLength();
						seq2Length = s2.getLength();

						pair = alignDNA(s1,s2);

					} else if (isRNASequence(str1) && isRNASequence(str2)) {
						RNASequence s1 = getRNASequence(str1);
						RNASequence s2 = getRNASequence(str2);
						seq1Length = s1.getLength();
						seq2Length = s2.getLength();

						pair = alignRNA(s1,s2);

					} else {
						logger.debug("Chains {},{} are either different kind of polymers or could not be recognized as protein or nucleotide polymers");
						continue;
					}


					int numGaps = getNumGaps(pair);
					int numGaps1 = getNumGapsQuery(pair);
					int numGaps2 = getNumGapsTarget(pair);

					int nonGaps = pair.getLength() - numGaps;

					double identity = (double)pair.getNumIdenticals()/(double)nonGaps;
					double gapCov1 = (double) numGaps1 / (double) seq1Length;
					double gapCov2 = (double) numGaps2 / (double) seq2Length;

					logger.debug("Alignment for chain pair {},{}: identity: {}, gap coverage 1: {}, gap coverage 2: {}",
							c1.getId(), c2.getId(), String.format("%4.2f",identity), String.format("%4.2f",gapCov1), String.format("%4.2f",gapCov2));
					logger.debug("\n"+pair.toString(100));

					if (identity > IDENTITY_THRESHOLD && gapCov1<GAP_COVERAGE_THRESHOLD && gapCov2<GAP_COVERAGE_THRESHOLD) {
						if (	!chainIds2compounds.containsKey(c1.getId()) &&
								!chainIds2compounds.containsKey(c2.getId())) {
							logger.debug("Creating Compound with chains {},{}",c1.getId(),c2.getId());

							EntityInfo ent = new EntityInfo();
							ent.addChain(c1);
							ent.addChain(c2);
							ent.setMolId(molId++);
							ent.setType(EntityType.POLYMER);
							chainIds2compounds.put(c1.getId(), ent);
							chainIds2compounds.put(c2.getId(), ent);


						} else {
							EntityInfo ent = chainIds2compounds.get(c1.getId());

							if (ent==null) {
								logger.debug("Adding chain {} to entity {}",c1.getId(),c2.getId());
								ent = chainIds2compounds.get(c2.getId());
								ent.addChain(c1);
								chainIds2compounds.put(c1.getId(), ent);

							} else {
								logger.debug("Adding chain {} to entity {}",c2.getId(),c1.getId());
								ent.addChain(c2);
								chainIds2compounds.put(c2.getId(), ent);

							}
						}
						if (!areResNumbersAligned(c1, c2)) {
							logger.warn("Including 100% identical chains {},{} in same Compound, although they have misaligned residue numbers",
									c1.getId(),c2.getId());
						}
					}

					if (identity>1) {
						logger.warn("Identity for chains {},{} above 1. {} identicals out of {} non-gap-aligned residues (identity {})",
								c1.getId(),c2.getId(),pair.getNumIdenticals(),nonGaps,identity);
						logger.warn("\n"+pair.toString(100));
					}

					if (chainIds2compounds.size()==polyChainIndices.size()) // we've got all chains in entities
						break outer;
				}
			}

		// anything not in a Compound will be its own Compound
		for (int i:polyChainIndices) {
			Chain c = allModels.get(0).get(i);
			if (!chainIds2compounds.containsKey(c.getId())) {
				logger.debug("Creating a 1-member Compound for chain {}",c.getId());

				EntityInfo ent = new EntityInfo();
				ent.addChain(c);
				ent.setMolId(molId++);

				chainIds2compounds.put(c.getId(),ent);
			}
		}


		return chainIds2compounds;
	}

	private SequencePair<ProteinSequence, AminoAcidCompound> alignProtein(ProteinSequence s1, ProteinSequence s2) {
		SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getIdentity();

		GapPenalty penalty = new SimpleGapPenalty(8, 1);

		PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> nw =
				Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.GLOBAL, penalty, matrix);

		return nw.getPair();
	}

	private SequencePair<DNASequence, NucleotideCompound> alignDNA(DNASequence s1, DNASequence s2) {
		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();

		GapPenalty penalty = new SimpleGapPenalty(8, 1);

		PairwiseSequenceAligner<DNASequence, NucleotideCompound> nw =
				Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.GLOBAL, penalty, matrix);

		return nw.getPair();
	}

	private SequencePair<RNASequence, NucleotideCompound> alignRNA(RNASequence s1, RNASequence s2) {
		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();

		GapPenalty penalty = new SimpleGapPenalty(8, 1);

		PairwiseSequenceAligner<RNASequence, NucleotideCompound> nw =
				Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.GLOBAL, penalty, matrix);

		return nw.getPair();
	}



	private static int getNumGaps(SequencePair<?, ?> pair) {
		int numGaps = 0;
		for (int alignmentIndex=1;alignmentIndex<=pair.getLength();alignmentIndex++) {
			if (pair.hasGap(alignmentIndex)) numGaps++;
		}
		return numGaps;
	}

	private static int getNumGapsQuery(SequencePair<?, ?> pair) {
		int numGaps = 0;
		for (int alignmentIndex=1;alignmentIndex<=pair.getLength();alignmentIndex++) {
			if (pair.getCompoundInQueryAt(alignmentIndex).getShortName().equals("-")) {
				numGaps++;
			}
		}
		return numGaps;
	}

	private static int getNumGapsTarget(SequencePair<?, ?> pair) {
		int numGaps = 0;
		for (int alignmentIndex=1;alignmentIndex<=pair.getLength();alignmentIndex++) {
			if (pair.getCompoundInTargetAt(alignmentIndex).getShortName().equals("-")) {
				numGaps++;
			}
		}
		return numGaps;
	}

	private boolean isProteinSequence(String str) {
		try {
			new ProteinSequence(str);
		} catch (CompoundNotFoundException e) {

			return false;
		}
		return true;
	}

	private boolean isDNASequence(String str) {
		try {
			new DNASequence(str);
		} catch (CompoundNotFoundException e) {

			return false;
		}
		return true;
	}

	private boolean isRNASequence(String str) {
		try {
			new RNASequence(str);
		} catch (CompoundNotFoundException e) {

			return false;
		}
		return true;

	}

	/**
	 * Returns the ProteinSequence or null if one can't be created
	 * @param str
	 * @return
	 */
	private ProteinSequence getProteinSequence(String str) {
		try {
			ProteinSequence s = new ProteinSequence(str);
			return s;
		} catch (CompoundNotFoundException e) {

			logger.error("Unexpected error when creating ProteinSequence",e);
		}
		return null;
	}

	/**
	 * Returns the DNASequence or null if one can't be created
	 * @param str
	 * @return
	 */
	private DNASequence getDNASequence(String str) {
		try {
			DNASequence s = new DNASequence(str);
			return s;
		} catch (CompoundNotFoundException e) {

			logger.error("Unexpected error when creating DNASequence ",e);
		}
		return null;
	}

	/**
	 * Returns the RNASequence or null if one can't be created
	 * @param str
	 * @return
	 */
	private RNASequence getRNASequence(String str) {
		try {
			RNASequence s = new RNASequence(str);
			return s;
		} catch (CompoundNotFoundException e) {

			logger.error("Unexpected error when creating RNASequence ",e);
		}
		return null;
	}

}
