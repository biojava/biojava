/*
 *                  BioJava development code
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
 * Created on Aug 23, 2007
 *
 */

package org.biojava.nbio.structure.io;


import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.io.mmcif.chem.ResidueType;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.*;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;


/** Aligns the SEQRES residues to the ATOM residues.
 * The AminoAcids that can be matched between the two of them will be set in the SEQRES
 * chains
 *
 *
 * @author Andreas Prlic
 *
 */
public class SeqRes2AtomAligner {

	private static final Logger logger = LoggerFactory.getLogger(SeqRes2AtomAligner.class);
	


	private String alignmentString;
	
	public SeqRes2AtomAligner(){
		logger.debug("initialising SeqRes2AtomAligner");
		alignmentString = "";
	}

	public String getAlignmentString() {
		return alignmentString;
	}

	public static Chain getMatchingAtomRes(Chain seqRes, List<Chain> atomList)
	{
		Iterator<Chain> iter = atomList.iterator();
		while(iter.hasNext()){
			Chain atomChain = iter.next();
			if ( atomChain.getChainID().equals(seqRes.getChainID())){
				return atomChain;
			}
		}

		logger.info("Could not match SEQRES chainID >" + seqRes.getChainID() + "< to ATOM chains!, size of atom chain: " + atomList.size());
		return null;
	}




	public void align(Structure s, List<Chain> seqResList){

		List<Chain> atomList   = s.getModel(0);


		for (Chain seqRes: seqResList){

				Chain atomRes = getMatchingAtomRes(seqRes,atomList);
				if ( atomRes == null)
					continue;

				mapSeqresRecords(atomRes,seqRes);

 
		}

	}

	/** 
	 * Map the seqRes groups to the atomRes chain.
	 * Updates the atomRes chain object with the mapped data
	 * The seqRes chain should not be needed after this and atomRes should be further used.
	 *  
	 * @param atomRes the chain containing ATOM groups (in atomGroups slot). This chain 
	 * is modified to contain in its seqresGroups slot the mapped atom groups
	 * @param seqRes the chain containing SEQRES groups (in atomGroups slot). This chain 
	 * is not modified
	 */
	public void mapSeqresRecords(Chain atomRes, Chain seqRes) {
		List<Group> seqResGroups = seqRes.getAtomGroups();
		List<Group> atmResGroups = atomRes.getAtomGroups();


		
		logger.debug("Comparing ATOM {} ({} groups) to SEQRES {} ({} groups) ", 
				atomRes.getChainID(), atmResGroups.size(), seqRes.getChainID(), seqResGroups.size());
		

		List<Group> matchedGroups = trySimpleMatch(seqResGroups, atmResGroups);

		if ( matchedGroups != null) {
			// update the new SEQRES list			
			atomRes.setSeqResGroups(matchedGroups);
			return;
		}

		logger.debug("Could not map SEQRES to ATOM records easily, need to align...");

		int numAminosSeqres = seqRes.getAtomGroups(GroupType.AMINOACID).size();
		int numNucleotidesSeqres = seqRes.getAtomGroups(GroupType.NUCLEOTIDE).size();

		if ( numAminosSeqres < 1) {
				
			if ( numNucleotidesSeqres > 1) {
				
				logger.debug("SEQRES chain {} is a nucleotide chain ({} nucleotides), aligning nucleotides...", seqRes.getChainID(), numNucleotidesSeqres);
				
				alignNucleotideChains(seqRes,atomRes);
				return;
			} else {
				
				logger.debug("SEQRES chain {} contains {} amino acids and {} nucleotides, ignoring...", seqRes.getChainID(),numAminosSeqres,numNucleotidesSeqres);
				
				return;
			}
		}

		if ( atomRes.getAtomGroups(GroupType.AMINOACID).size() < 1) {
			logger.debug("ATOM chain {} does not contain amino acids, ignoring...", atomRes.getChainID());
			return;
		}
		
		logger.debug("Proceeding to do protein alignment for chain {}", atomRes.getChainID() );


		boolean noMatchFound = alignProteinChains(seqResGroups,atomRes.getAtomGroups());
		if ( ! noMatchFound){
			atomRes.setSeqResGroups(seqResGroups);
		}

	}

	private void alignNucleotideChains(Chain seqRes, Chain atomRes) {

		if ( atomRes.getAtomGroups(GroupType.NUCLEOTIDE).size() < 1) {
			logger.debug("ATOM chain {} does not contain nucleotides, ignoring...", atomRes.getChainID());
			
			return;
		}
		logger.debug("Alignment for chain {}", atomRes.getChainID() );

		List<Group> seqResGroups = seqRes.getAtomGroups();
		boolean noMatchFound = alignNucleotideGroups(seqResGroups,atomRes.getAtomGroups());
		if ( ! noMatchFound){
			atomRes.setSeqResGroups(seqResGroups);
		}

	}

	/** 
	 * A simple matching approach that tries to do a 1:1 mapping between SEQRES and ATOM records
	 *  
	 * @param seqRes
	 * @param atomList
	 * @return the matching or null if the matching didn't work
	 */
	private List<Group> trySimpleMatch(List<Group> seqResGroups,List<Group> atmResGroups) {
		// by default first ATOM position is 1
		//

		@SuppressWarnings("unchecked")
		List<Group> newSeqResGroups = (ArrayList<Group>)((ArrayList<Group>)seqResGroups).clone();

		boolean startAt1 = true;

		for ( int atomResPos = 0 ; atomResPos < atmResGroups.size() ; atomResPos++){

			// let's try to match this case
			Group atomResGroup = atmResGroups.get(atomResPos);

			// let's ignore waters
			if ( atomResGroup.isWater()){
				continue;
			}

			ResidueNumber atomResNum = atomResGroup.getResidueNumber();

			int seqResPos = atomResNum.getSeqNum();



			if ( seqResPos < 0) {
				logger.debug("ATOM residue number < 0 : {}", seqResPos);
				return null;
			}

			if ( seqResPos == 0){
				// make sure the first SEQRES is matching.
				Group seqResGroup = seqResGroups.get(0);
				if (  seqResGroup.getPDBName().equals(atomResGroup.getPDBName())){
					// they match! seems in this case the numbering starts with 0...
					startAt1 = false;
				} else {
					
					logger.debug("SEQRES position 1 ({}) does not match ATOM PDB res num 0 ({})",
							seqResGroup.getPDBName(), atomResGroup.getPDBName());

					
					return null;

				}
			}


			if ( startAt1 )
				seqResPos--;

			// another check that would require the alignment approach
			if ( startAt1 && seqResPos >=  seqResGroups.size()  )
			{

				// could be a HETATOM...
				if ( atomResGroup instanceof AminoAcid) {
					logger.debug(" ATOM residue nr: " + seqResPos + " > seqres! " + seqResGroups.size() + " " + atomResGroup);
					return null;
				} else if ( atomResGroup instanceof NucleotideImpl) {
					logger.debug(" NUCLEOTIDE residue nr: " + seqResPos + " > seqres! " + seqResGroups.size() + " " + atomResGroup);
					return null;
				} else {
					// we won't map HETATOM groups...
					continue;
				}
			}


			//			if ( seqResPos < 0){
			//
			//				System.err.println("What is going on??? " + atomRes.getChainID() + " " + atomResGroup);
			//			}

			if ( seqResPos >= seqResGroups.size()){
				logger.debug("seqres groups don't match atom indices " + seqResPos);
				if ( atomResGroup instanceof AminoAcid )
					return null;
				else
					continue;
			}

			Group seqResGroup = seqResGroups.get(seqResPos );

			if ( ! seqResGroup.getPDBName().trim().equals(atomResGroup.getPDBName().trim())){
				// a mismatch! something is wrong in the mapping and we need to do an alignment
				logger.debug("Mismatch of SEQRES pos " + seqResPos + " and ATOM record: " + atomResGroup + " | " + seqResGroup);
				return null;
			}

			// the two groups are identical and we can merge them
			// replace the SEQRES group with the ATOM group...

			Group replacedGroup = newSeqResGroups.set(seqResPos, atomResGroup);
			logger.debug("Merging index {}: replaced seqres group {} ({}) with atom group {} ({})", 
					seqResPos, 
					replacedGroup.getResidueNumber(), replacedGroup.getPDBName(), 
					atomResGroup.getResidueNumber(), atomResGroup.getPDBName());

		}

		// all went ok. copy over the updated list to the original one.
		// note: if something went wrong, we did not modifiy the original list.
		//seqResGroups = newSeqResGroups;

		
		//			int pos = -1;
		//			for (Group g: seqResGroups){
		//				pos++;
		//				logger.debug(pos + " " + g);
		//			}
		
		//System.out.println("I:" + seqResGroups);
		// all atom records could get matched correctly!
		return newSeqResGroups;

	}

	/** 
	 * Returns the full sequence of the Atom records of a parent
	 * with X instead of HETATMSs. The advantage of this is
	 * that it allows us to also align HETATM groups back to the SEQRES.
	 * @param groups the list of groups in a parent
	 * @param positionIndex a Map to keep track of which group is at which sequence position
	 * @return string representations
	 */
	public static String getFullAtomSequence(List<Group> groups, Map<Integer, Integer> positionIndex){

		StringBuffer sequence = new StringBuffer() ;
		int seqIndex = 0; // track sequence.length()
		for ( int i=0 ; i< groups.size(); i++){
			Group g = groups.get(i);
			
			if ( g instanceof AminoAcid ){
				AminoAcid a = (AminoAcid)g;
				char oneLetter =a.getAminoType();
				if ( oneLetter == '?')
					oneLetter = 'X';

				positionIndex.put(seqIndex,i);
				sequence.append(oneLetter );
				seqIndex++;
			} else {

				// exclude solvents
				if (  g.isWater())
					continue;

				// exclude metals
				if ( g.size() == 1 ) {

					Atom a = g.getAtom(0);
					if ( a == null)
						continue;
					if (a.getElement().isMetal())
						continue;

				}
				
				ChemComp cc = g.getChemComp();
				if ( cc == null) {
					logger.debug("No chem comp available for group {}",g.toString());
					// not sure what to do in that case!
					continue;
				}
				if (
						ResidueType.lPeptideLinking.equals(cc.getResidueType()) ||
						PolymerType.PROTEIN_ONLY.contains(cc.getPolymerType())  ||
						PolymerType.POLYNUCLEOTIDE_ONLY.contains(cc.getPolymerType())
						) {
					//System.out.println(cc.getOne_letter_code());
					String c = cc.getOne_letter_code();
					if ( c.equals("?"))
						c = "X";

					// For some unusual cases the het residue can map to 2 or more standard aas and thus give an 
					// insertion of length >1. 
					//      e.g. 1: SUI maps to DG  (in 1oew,A)
					//		e.g. 2: NRQ maps to MYG (in 3cfh,A)
					if (c.length()>1) {
						logger.info("Group '{}' maps to more than 1 standard aminoacid: {}.",
								g.toString(), c);
					}
					// because of the mapping to more than 1 aminoacid, we have 
					// to loop through it (99% of cases c will have length 1 anyway)
					for (int cIdx=0;cIdx<c.length();cIdx++) {
						positionIndex.put(seqIndex,i);
						sequence.append(c.charAt(cIdx));
						seqIndex++;
					}
				} else {
					logger.debug("Group {} is not lPeptideLinked, nor PROTEIN_ONLY, nor POLYNUCLEOTIDE_ONLY",g.toString());
					continue;
				}


				//sequence.append("X");
			}

		}

		return sequence.toString();

	}


	private boolean alignNucleotideGroups(List<Group> seqRes, List<Group> atomRes) {

		Map<Integer,Integer> seqresIndexPosition = new HashMap<Integer, Integer>();
		Map<Integer,Integer> atomIndexPosition   = new HashMap<Integer, Integer>();

		String seq1 = getFullAtomSequence(seqRes, seqresIndexPosition);
		//
		String seq2 = getFullAtomSequence(atomRes, atomIndexPosition);

		if (seq1.isEmpty() || seq2.isEmpty()) {
			logger.warn("Could not align nucleotide sequences, at least one of them is empty");
			return true;
		}

		logger.debug("align seq1 ("+ seq1.length()+") " + seq1);
		logger.debug("align seq2 ("+ seq2.length()+") " + seq2);
		


		Sequence<NucleotideCompound> s1;
		Sequence<NucleotideCompound> s2;

		try {
			s1 = new DNASequence(seq1,AmbiguityDNACompoundSet.getDNACompoundSet());
		} catch (CompoundNotFoundException e){
			///oops perhaps a RNA sequence?
			try {
				s1 = new RNASequence(seq1,AmbiguityRNACompoundSet.getRNACompoundSet());
			} catch (CompoundNotFoundException ex) {
				logger.warn("Could not determine compound set for sequence 1 " + seq1);
				return true;
			}
		}
		try {
			s2 = new DNASequence(seq2,AmbiguityDNACompoundSet.getDNACompoundSet());
		} catch (CompoundNotFoundException e){
			///oops perhaps a RNA sequence?
			try {
				s2 = new RNASequence(seq2,AmbiguityRNACompoundSet.getRNACompoundSet());
			} catch (CompoundNotFoundException ex) {
				logger.warn("Could not determine compound set for sequence 2 " + seq2);
				return true;
			}
		}

		if ( ! s1.getCompoundSet().equals(s2.getCompoundSet()) ) {
			// e.g. trying to align a DNA and an RNA sequence...
			// test PDB ID: 1A34
			// try to make both RNA sequence...
			if ( ! s1.getCompoundSet().equals(AmbiguityRNACompoundSet.getRNACompoundSet())) {
				try {
					s1 = new RNASequence(seq1,AmbiguityRNACompoundSet.getRNACompoundSet());
				} catch (CompoundNotFoundException ex) {
					logger.warn("Could not align DNA and RNA compound sets: " + seq1);
					return true;
				}
			}
			
			if( ! s2.getCompoundSet().equals(AmbiguityRNACompoundSet.getRNACompoundSet())) {
					try {
						s2 = new RNASequence(seq2,AmbiguityRNACompoundSet.getRNACompoundSet());
					} catch (CompoundNotFoundException ex) {
						logger.warn("Could not align DNA and RNA compound sets: " + seq2);
						return true;
					}
				}
		}
		
		
		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();

		GapPenalty penalty = new SimpleGapPenalty(8,1);

		PairwiseSequenceAligner<Sequence<NucleotideCompound>, NucleotideCompound> smithWaterman =
				Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.LOCAL, penalty, matrix);

		

		SequencePair<Sequence<NucleotideCompound>, NucleotideCompound> pair = smithWaterman.getPair();  



		if ( pair == null) {
			logger.warn("Could not align nucleotide sequences. ATOM and SEQRES groups will not be aligned.");
			logger.warn("Sequences: ");
			logger.warn(seq1);
			logger.warn(seq2);
			return true;

		}

		
		
		logger.debug("Alignment:\n"+pair.toString(100));
		

		boolean noMatchFound = mapDNAChains(seqRes,atomRes,pair,seqresIndexPosition, atomIndexPosition );

		return noMatchFound;

	}


	/** 
	 * Aligns two chains of groups, where the first parent is representing the
	 * list of amino acids as obtained from the SEQRES records, and the second parent
	 * represents the groups obtained from the ATOM records (and containing the actual ATOM information).
	 * This does the actual alignment and if a group can be mapped to a position in the SEQRES then the corresponding
	 * position is replaced with the group that contains the atoms.
	 *
	 * @param seqRes
	 * @param atomRes
	 * @return true if no match has been found
	 */
	private boolean alignProteinChains(List<Group> seqRes, List<Group> atomRes) {

		Map<Integer,Integer> seqresIndexPosition = new HashMap<Integer, Integer>();
		Map<Integer,Integer> atomIndexPosition   = new HashMap<Integer, Integer>();

		String seq1 = getFullAtomSequence(seqRes, seqresIndexPosition);
		//
		String seq2 = getFullAtomSequence(atomRes, atomIndexPosition);

		
		logger.debug("Protein seq1 to align (length "+ seq1.length()+"): " + seq1);
		logger.debug("Protein seq2 to align (length "+ seq2.length()+"): " + seq2);
		
		ProteinSequence s1;
		ProteinSequence s2;
		try {
			s1 = new ProteinSequence(seq1);
			s2 = new ProteinSequence(seq2);
		} catch (CompoundNotFoundException e) {
			logger.warn("Could not create protein sequences ({}) to align ATOM and SEQRES groups, they will remain unaligned.", e.getMessage());
			return true;
		}
		

		SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum65();

		GapPenalty penalty = new SimpleGapPenalty(8, 1);


		PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> smithWaterman =
				Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.LOCAL, penalty, matrix);

		SequencePair<ProteinSequence, AminoAcidCompound> pair = smithWaterman.getPair();



		if ( pair == null) {
			logger.warn("Could not align protein sequences. ATOM and SEQRES groups will not be aligned.");
			logger.warn("Sequences: ");
			logger.warn(seq1);
			logger.warn(seq2);
			return true;
		}			


		logger.debug("Alignment:\n"+pair.toString(100));
		

		boolean noMatchFound = mapChains(seqRes,atomRes,pair,seqresIndexPosition, atomIndexPosition );

		return noMatchFound;


	}


	private boolean mapChains(List<Group> seqResGroups, List<Group> atomRes,
			SequencePair<ProteinSequence, AminoAcidCompound> pair,
			Map<Integer,Integer> seqresIndexPosition,
			Map<Integer,Integer> atomIndexPosition )   {



		// at the present stage the future seqREs are still stored as Atom groups in the seqRes parent...


		int aligLength = pair.getLength();

		//		System.out.println("sequence length seqres:");
		//		System.out.println(seqresIndexPosition.keySet().size());
		//		System.out.println("alignment length: " + aligLength);
		// System.out.println(gapSymbol.getName());

		// make sure we actually find an alignment
		boolean noMatchFound = true;

		Compound gapSymbol =  AminoAcidCompoundSet.getAminoAcidCompoundSet().getCompoundForString("-");

		mainLoop:
			for (int i = 1; i <= aligLength ; i++) {

				Compound s =  pair.getCompoundAt(1, i);
				Compound a =  pair.getCompoundAt(2,i);

				// alignment is using internal index start at 1...
				int posSeq = pair.getIndexInQueryAt(i)  -1;
				int posAtom = pair.getIndexInTargetAt(i) -1;

				if (  s.equals(gapSymbol) || a.equals(gapSymbol)){
					continue;				
				}

				if ( s.equals(a)){

					// the atom record can be aligned to the SeqRes record!
					// replace the SeqRes group with the Atom group!

					Group s1 = seqResGroups.get(seqresIndexPosition.get(posSeq));
					Group a1 = atomRes.get(atomIndexPosition.get(posAtom));

					if ( s1 == null || a1 == null){
						/// can't map this position...
						logger.warn("can't map " + i + ":" + s + " " + posSeq +" " + s1 + " atom: " + posAtom + " " + a1 );
						continue mainLoop;
					}

					// need to trim the names to allow matching e.g in
					// pdb1b2m
					String pdbNameS = s1.getPDBName();
					String pdbNameA = a1.getPDBName();
					if ( pdbNameS == null || pdbNameA == null ){
						logger.warn("null value for group.getPDBName found at {} when trying to align {} and {} {}",posSeq, s1, a1, posAtom);
						logger.warn("ATOM and SEQRES sequences will not be aligned.");
						return true;
					}
					if (! pdbNameA.equals(pdbNameS)){
						if ( ! pdbNameA.trim().equals(pdbNameS.trim())) {
							logger.info(s1 + " " + posSeq + " does not align with " + a1+ " " + posAtom + " should be: " + s + " : " + a);
							if ( s1.getType().equals(HetatomImpl.type) && a1.getType().equals(HetatomImpl.type)){
								logger.info("they seem to be hetatoms, so ignoring mismatch.");
							}
							else {
								//System.out.println(lst1.seqString());
								//System.out.println(lst2.seqString());
								logger.warn("Could not match residues {} {}", s1, a1);
							}

						}
					}

					// do the actual replacing of the SEQRES group with the ATOM group
					//					if ( s1.getChain().getChainID().equals("A")) {
					//						System.out.println(" setting " + posSeq +" " + a1);
					//					}
					seqResGroups.set(seqresIndexPosition.get(posSeq),a1);
					noMatchFound = false;
				}
			}


		// now we merge the two chains into one
		// the Groups that can be aligned are now pointing to the
		// groups in the Atom records.
		if (  noMatchFound) {

			logger.debug("no alignment found!");
		}
		return noMatchFound;

	}

	private boolean mapDNAChains(List<Group> seqResGroups, List<Group> atomRes,
			SequencePair<Sequence<NucleotideCompound>, NucleotideCompound> pair,
			Map<Integer,Integer> seqresIndexPosition,
			Map<Integer,Integer> atomIndexPosition)   {


		// at the present stage the future seqREs are still stored as Atom groups in the seqRes parent...


		int aligLength = pair.getLength();

		//		System.out.println("sequence length seqres:");
		//		System.out.println(seqresIndexPosition.keySet().size());
		//		System.out.println("alignment length: " + aligLength);
		// System.out.println(gapSymbol.getName());

		// make sure we actually find an alignment
		boolean noMatchFound = true;

		Compound gapSymbol =  DNACompoundSet.getDNACompoundSet().getCompoundForString("-");

		mainLoop:
			for (int i = 1; i <= aligLength ; i++) {

				Compound s =  pair.getCompoundAt(1, i);
				Compound a =  pair.getCompoundAt(2,i);

				// alignment is using internal index start at 1...
				int posSeq = pair.getIndexInQueryAt(i)  -1;
				int posAtom = pair.getIndexInTargetAt(i) -1;

				if (  s.equals(gapSymbol) || a.equals(gapSymbol)){
					continue;				
				}

				if ( s.equals(a)){

					// the atom record can be aligned to the SeqRes record!
					// replace the SeqRes group with the Atom group!

					Group s1 = seqResGroups.get(seqresIndexPosition.get(posSeq));
					Group a1 = atomRes.get(atomIndexPosition.get(posAtom));

					if ( s1 == null || a1 == null){
						/// can't map this position...
						logger.warn("can't map " + i + ":" + s + " " + posSeq +" " + s1 + " atom: " + posAtom + " " + a1 );
						continue mainLoop;
					}

					// need to trim the names to allow matching e.g in
					// pdb1b2m
					String pdbNameS = s1.getPDBName();
					String pdbNameA = a1.getPDBName();
					if ( pdbNameS == null || pdbNameA == null ){
						logger.warn("null value found for group.getPDBName() at " + posSeq + " when trying to align " + s1 + " and " + a1 + " " + posAtom);
						logger.warn("ATOM and SEQRES sequences will not be aligned.");
					}
					if (! pdbNameA.equals(pdbNameS)){
						if ( ! pdbNameA.trim().equals(pdbNameS.trim())) {
							logger.info(s1 + " " + posSeq + " does not align with " + a1+ " " + posAtom + " should be: " + s + " : " + a);
							if ( s1.getType().equals(HetatomImpl.type) && a1.getType().equals(HetatomImpl.type)){
								logger.info("they seem to be hetatoms, so ignoring mismatch.");
							}
							else {
								//  System.exit(0);// for debug only
								//System.out.println(lst1.seqString());
								//System.out.println(lst2.seqString());
								logger.warn("Could not match residues " + s1 + " " + a1);
							}

						}
					}

					// do the actual replacing of the SEQRES group with the ATOM group
					//					if ( s1.getChain().getChainID().equals("A")) {
					//						System.out.println(" setting " + posSeq +" " + a1);
					//					}
					seqResGroups.set(seqresIndexPosition.get(posSeq),a1);
					noMatchFound = false;
				}
			}


		// now we merge the two chains into one
		// the Groups that can be aligned are now pointing to the
		// groups in the Atom records.
		if (  noMatchFound) {

			logger.debug("no alignment found!");
		}
		return noMatchFound;

	}

}
