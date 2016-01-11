/**
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
 * Created by Spencer Bliven
 *
 */
package org.biojava.nbio.structure.io;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SimpleSubstitutionMatrix;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * A utility class with methods for matching ProteinSequences with
 * Structures.
 * @author Spencer Bliven
 *
 */
public class StructureSequenceMatcher {
	
	private static final Logger logger = LoggerFactory.getLogger(StructureSequenceMatcher.class);

	/**
	 * Get a substructure of {@code wholeStructure} containing only the {@link Group Groups} that are included in
	 * {@code sequence}. The resulting structure will contain only {@code ATOM} residues; the SEQ-RES will be empty.
	 * The {@link Chain Chains} of the Structure will be new instances (cloned), but the {@link Group Groups} will not.
	 * @param sequence The input protein sequence
	 * @param wholeStructure The structure from which to take a substructure
	 * @return The resulting structure
	 * @throws StructureException
	 * @see {@link #matchSequenceToStructure(ProteinSequence, Structure)}
	 */
	public static Structure getSubstructureMatchingProteinSequence(ProteinSequence sequence, Structure wholeStructure) {
		ResidueNumber[] rns = matchSequenceToStructure(sequence, wholeStructure);
		Structure structure = wholeStructure.clone();
		structure.getChains().clear();
//		structure.getHetGroups().clear();
		Chain currentChain = null;
		for (ResidueNumber rn : rns) {
			if (rn == null) continue;
			Group group; // note that we don't clone
			try {
				group = StructureTools.getGroupByPDBResidueNumber(wholeStructure, rn);
			} catch (StructureException e) {
				throw new IllegalArgumentException("Could not find residue " + rn + " in structure", e);
			}
			Chain chain = new ChainImpl();
			chain.setChainID(group.getChainId());
			if (currentChain == null || !currentChain.getChainID().equals(chain.getChainID())) {
				structure.addChain(chain);
				chain.setCompound(group.getChain().getCompound());
				chain.setStructure(structure);
				chain.setSwissprotId(group.getChain().getSwissprotId());
				chain.setInternalChainID(group.getChain().getInternalChainID());
				chain.setId(group.getChain().getId());
				currentChain = chain;
			}
			currentChain.addGroup(group);
		}
		return structure;
	}
	
	/**
	 * Generates a ProteinSequence corresponding to the sequence of struct,
	 * and maintains a mapping from the sequence back to the original groups.
	 * 
	 * Chains are appended to one another. 'X' is used for heteroatoms.
	 * 
	 * @param struct Input structure
	 * @param groupIndexPosition An empty map, which will be populated with
	 *  (residue index in returned ProteinSequence) -> (Group within struct)
	 * @return A ProteinSequence with the full sequence of struct. Chains are
	 *  concatenated in the same order as the input structures
	 *  
	 * @see {@link SeqRes2AtomAligner#getFullAtomSequence(List, Map)}, which
	 * 	does the heavy lifting.
	 * 
	 */
	public static ProteinSequence getProteinSequenceForStructure(Structure struct, Map<Integer,Group> groupIndexPosition ) {

		if( groupIndexPosition != null) {
			groupIndexPosition.clear();
		}

		StringBuilder seqStr = new StringBuilder();

		for(Chain chain : struct.getChains()) {
			List<Group> groups = chain.getAtomGroups();
			Map<Integer,Integer> chainIndexPosition = new HashMap<Integer, Integer>();
			int prevLen = seqStr.length();

			// get the sequence for this chain
			String chainSeq = SeqRes2AtomAligner.getFullAtomSequence(groups, chainIndexPosition);
			seqStr.append(chainSeq);


			// fix up the position to include previous chains, and map the value back to a Group
			for(Integer seqIndex : chainIndexPosition.keySet()) {
				Integer groupIndex = chainIndexPosition.get(seqIndex);
				groupIndexPosition.put(prevLen + seqIndex, groups.get(groupIndex));
			}
		}

		ProteinSequence s = null;
		try {
			s = new ProteinSequence(seqStr.toString());
		} catch (CompoundNotFoundException e) {
			// I believe this can't happen, please correct this if I'm wrong - JD 2014-10-24
			// we can log an error if it does, it would mean there's a bad bug somewhere
			logger.error("Could not create protein sequence, unknown compounds in string: {}", e.getMessage());
		}
		
		return s;	
	}

	/**
	 * Given a sequence and the corresponding Structure, get the ResidueNumber
	 * for each residue in the sequence.
	 * 
	 * <p>Smith-Waterman alignment is used to match the sequences. Residues
	 * in the sequence but not the structure or mismatched between sequence
	 * and structure will have a null atom, while
	 * residues in the structure but not the sequence are ignored with a warning.
	 * @param seq The protein sequence. Should match the sequence of struct very
	 * 	closely.
	 * @param struct The corresponding protein structure
	 * @return A list of ResidueNumbers of the same length as seq, containing
	 *  either the corresponding residue or null.
	 */
	public static ResidueNumber[] matchSequenceToStructure(ProteinSequence seq, Structure struct) {

		//1. Create ProteinSequence for struct while remembering to which group each residue corresponds
		Map<Integer,Group> atomIndexPosition   = new HashMap<Integer, Group>();

		ProteinSequence structSeq = getProteinSequenceForStructure(struct,atomIndexPosition);

		// TODO This should really be semi-global alignment, though global for either sequence OR structure (we don't know which)
		//2. Run Smith-Waterman to get the alignment
		// Identity substitution matrix with +1 for match, -1 for mismatch
		// TODO
		SubstitutionMatrix<AminoAcidCompound> matrix = 
			new SimpleSubstitutionMatrix<AminoAcidCompound>(
					AminoAcidCompoundSet.getAminoAcidCompoundSet(),
					(short)1, (short)-1 );
		matrix = new SimpleSubstitutionMatrix<AminoAcidCompound>(
				AminoAcidCompoundSet.getAminoAcidCompoundSet(),
				new InputStreamReader(
						SimpleSubstitutionMatrix.class.getResourceAsStream("/matrices/blosum100.txt")),
		"blosum100");
		SequencePair<ProteinSequence, AminoAcidCompound> pair = 
			Alignments.getPairwiseAlignment(seq, structSeq,
					PairwiseSequenceAlignerType.GLOBAL, new SimpleGapPenalty(), matrix);

		//System.out.print(pair.toString());

		//3. Convert the alignment back to Atoms
		AlignedSequence<ProteinSequence,AminoAcidCompound> alignedSeq = pair.getQuery();
		AlignedSequence<ProteinSequence,AminoAcidCompound> alignedStruct = pair.getTarget();


		assert(alignedSeq.getLength() == alignedStruct.getLength());
		
//		System.out.println(pair.toString(80));
//		System.out.format("%d/min{%d,%d}\n", pair.getNumIdenticals(),
//				alignedSeq.getLength()-alignedSeq.getNumGaps(),
//				alignedStruct.getLength()-alignedStruct.getNumGaps());

		ResidueNumber[] ca = new ResidueNumber[seq.getLength()];

		for( int pos = alignedSeq.getStart().getPosition(); pos <= alignedSeq.getEnd().getPosition(); pos++ ) { // 1-indexed
			//skip missing residues from sequence. These probably represent alignment errors
			if(alignedSeq.isGap(pos)) {
				int structIndex = alignedStruct.getSequenceIndexAt(pos)-1;
				assert(structIndex > 0);//should be defined since seq gap

				Group g = atomIndexPosition.get(structIndex);

				logger.warn("Chain {} residue {} in the Structure {} has no corresponding amino acid in the sequence.",
						g.getChainId(),
						g.getResidueNumber().toString(),
						g.getChain().getStructure().getPDBCode() );
				continue;
			}

			if(! alignedStruct.isGap(pos) ) {
				int seqIndex = alignedSeq.getSequenceIndexAt(pos)-1;//1-indexed
				int structIndex = alignedStruct.getSequenceIndexAt(pos)-1;//1-indexed
				Group g = atomIndexPosition.get(structIndex);
				
				assert(0<=seqIndex && seqIndex < ca.length);
				
				ca[seqIndex] = g.getResidueNumber(); //remains null for gaps
			}
		}
		return ca;
	}
	
	
	/**
	 * Removes all gaps ('-') from a protein sequence
	 * @param gapped
	 * @return
	 */
	public static ProteinSequence removeGaps(ProteinSequence gapped) {
		final String[] gapStrings = {"-","."};
		
		StringBuilder seq = new StringBuilder();
		
		CompoundSet<AminoAcidCompound> aaSet = gapped.getCompoundSet();
		AminoAcidCompound[] gaps = new AminoAcidCompound[gapStrings.length];
		
		for(int i=0;i<gapStrings.length;i++) {
			gaps[i] = aaSet.getCompoundForString(gapStrings[i]);
		}
		
		for(int i=1; i<=gapped.getLength();i++) { //1-indexed
			AminoAcidCompound aa = gapped.getCompoundAt(i);
			boolean isGap = false;
			for(AminoAcidCompound gap : gaps) {
				if( aa.equals(gap)) {
					isGap = true;
					break;
				}
			}
			if(!isGap) {
				seq.append(aa.getShortName());
			}
		}
		
		ProteinSequence ungapped = null;
		try {
			ungapped = new ProteinSequence(seq.toString());
		} catch (CompoundNotFoundException e) {
			// this can't happen, if it does there's a bug somewhere
			logger.error("Could not create ungapped protein sequence, found unknown compounds: {}. This is most likely a bug.", e.getMessage());
		}
		
		return ungapped;
	}
	
	/**
	 * Creates a new list consisting of all columns of gapped where no row 
	 * contained a null value.
	 * 
	 * Here, "row" refers to the first index and "column" to the second, eg
	 * gapped.get(row).get(column)
	 * @param gapped A rectangular matrix containing null to mark gaps
	 * @return A new List without columns containing nulls
	 */
	public static <T> T[][] removeGaps(final T[][] gapped) {
		if(gapped == null ) return null;
		if(gapped.length < 1) return Arrays.copyOf(gapped, gapped.length);
		
		final int nProts = gapped.length;
		final int protLen = gapped[0].length; // length of gapped proteins
		
		// Verify that input is rectangular
		for(int i=0;i<nProts;i++) {
			if(gapped[i].length != protLen) {
				throw new IllegalArgumentException(String.format(
						"Expected a rectangular array, but row 0 has %d elements " +
						"while row %d has %d.", protLen,i,gapped[i].length));
				
			}
		}
		
		// determine where gaps exist in any structures
		boolean[] isGap = new boolean[protLen];
		int gaps = 0;
		for(int j=0;j<protLen;j++) {
			for(int i=0;i<nProts;i++) {
				if(gapped[i][j] == null ) {
					isGap[j] = true;
					gaps++;
					break; //go to next position
				}
			}
		}
		
		// Create ungapped array
		T[][] ungapped = Arrays.copyOf(gapped,nProts);
		final int ungappedLen = protLen-gaps;
		for(int i=0;i<nProts;i++) {
			ungapped[i] = Arrays.copyOf(gapped[i],ungappedLen);
			int k = 0;
			for(int j=0;j<protLen;j++) {
				if(!isGap[j]) { //skip gaps
					assert(gapped[i][j] != null);
					ungapped[i][k] = gapped[i][j];
					k++;
				}
			}
			assert(k == ungappedLen);
		}
		
		return ungapped;
	}
}
