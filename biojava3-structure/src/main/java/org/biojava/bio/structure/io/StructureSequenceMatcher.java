package org.biojava.bio.structure.io;

import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SimpleSubstitutionMatrix;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.template.CompoundSet;


/**
 * A utility class with methods for matching ProteinSequences with
 * Structures.
 * @author Spencer Bliven
 *
 */
public class StructureSequenceMatcher {

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


		return new ProteinSequence(seqStr.toString());	
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
						SimpleSubstitutionMatrix.class.getResourceAsStream("/blosum100.txt")),
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

				System.err.format("Warning: chain %s residue %s in the Structure %s has no corresponding amino acid in the sequence.\n",
						g.getChainId(),
						g.getResidueNumber().toString(),
						g.getChain().getParent().getPDBCode());
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
		
		ProteinSequence ungapped = new ProteinSequence(seq.toString());
		
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
