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

package org.biojava.bio.structure.io;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;
import org.biojava.bio.structure.HetatomImpl;
import org.biojava.bio.structure.NucleotideImpl;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.io.mmcif.chem.PolymerType;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.PairwiseSequenceAligner;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.exceptions.CompoundNotFoundException;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.RNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.AmbiguityRNACompoundSet;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


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
	
	static final List<String> excludeTypes;


	String alignmentString;
	static {
		excludeTypes = new ArrayList<String>();
		excludeTypes.add("HOH"); // we don't want to align water against the SEQRES...
		excludeTypes.add("DOD"); // deuterated water


		//matrix.printMatrix();

	}

	public SeqRes2AtomAligner(){
		logger.debug("initialising SeqRes2AtomAligner");
		alignmentString = "";
	}

	public String getAlignmentString() {
		return alignmentString;
	}

	public Chain getMatchingAtomRes(Chain seqRes, List<Chain> atomList)
	{
		Iterator<Chain> iter = atomList.iterator();
		while(iter.hasNext()){
			Chain atomChain = iter.next();
			if ( atomChain.getChainID().equals(seqRes.getChainID())){
				return atomChain;
			}
		}

		logger.debug("Could not match SEQRES chainID >" + seqRes.getChainID() + "< to ATOM chains!, size of atom chain: " + atomList.size());
		return null;
	}




	public void align(Structure s, List<Chain> seqResList){

		//List<Chain> seqResList = s.getSeqRes();
		List<Chain> atomList   = s.getModel(0);


		for (Chain seqRes: seqResList){
			try {
				Chain atomRes = getMatchingAtomRes(seqRes,atomList);
				if ( atomRes == null)
					continue;
				mapSeqresRecords(atomRes,seqRes);

				//chains.add(mapped); 
			} catch (StructureException e){
				logger.error("Exception caught ", e);
			}
		}
		//s.setChains(0,chains);

	}

	/** map the seqres groups to the atomRes chain.
	 *  updates the atomRes chain object with the mapped data
	 *  seqRes chain should not be needed after this and atomRes should be continued to be used.
	 *  
	 * @param atomRes
	 * @param seqRes
	 * @throws StructureException 
	 */
	public void mapSeqresRecords(Chain atomRes, Chain seqRes) throws StructureException {
		List<Group> seqResGroups = seqRes.getAtomGroups();
		List<Group> atmResGroups = atomRes.getAtomGroups();


		
		logger.debug("Comparing ATOM {} ({}) to SEQRES {} ({}) ", 
				atomRes.getChainID(), atmResGroups.size(), seqRes.getChainID(), seqResGroups.size());
		

		List<Group> matchedGroups = trySimpleMatch(seqResGroups, atmResGroups);

		if ( matchedGroups != null) {
			// update the new SEQRES list			
			atomRes.setSeqResGroups(matchedGroups);
			return;
		}

		//System.err.println("Could not map SEQRES to ATOM records easily, need to align...");

		if ( seqRes.getAtomGroups(GroupType.AMINOACID).size() < 1) {

			if ( seqRes.getAtomGroups(GroupType.NUCLEOTIDE).size() > 1) {
				
				logger.debug("Chain {} is a nucleotide chain, aligning nucleotides...", seqRes.getChainID());
				
				align2NucleotideChains(seqRes,atomRes);
				return;
			} else {
				
				logger.debug("Chain {} does not contain amino acids, ignoring...", seqRes.getChainID());
				
				return;
			}
		}

		if ( atomRes.getAtomGroups(GroupType.AMINOACID).size() < 1) {
			logger.debug("Chain {} does not contain amino acids, ignoring...", atomRes.getChainID());
			return;
		}
		
		logger.debug("Alignment for chain {}", atomRes.getChainID() );


		boolean noMatchFound = align(seqResGroups,atomRes.getAtomGroups());
		if ( ! noMatchFound){
			atomRes.setSeqResGroups(seqResGroups);
		}

	}

	private void align2NucleotideChains(Chain seqRes, Chain atomRes) throws StructureException {

		if ( atomRes.getAtomGroups(GroupType.NUCLEOTIDE).size() < 1) {
			logger.debug("Chain {} does not contain nucleotides, ignoring...", atomRes.getChainID());
			
			return;
		}
		logger.debug("Alignment for chain {}", atomRes.getChainID() );

		List<Group> seqResGroups = seqRes.getAtomGroups();
		boolean noMatchFound = alignNucleotideGroups(seqResGroups,atomRes.getAtomGroups());
		if ( ! noMatchFound){
			atomRes.setSeqResGroups(seqResGroups);
		}

	}

	/** a simple matching approach that tries to do a 1:1 mapping between SEQRES and ATOM records
	 *  returns true if this simple matching approach worked fine
	 *  
	 * @param seqRes
	 * @param atomList
	 * @return
	 */
	private List<Group> trySimpleMatch(List<Group> seqResGroups,List<Group> atmResGroups) {
		// by default first ATOM position is 1
		//

		@SuppressWarnings("unchecked")
		List<Group> newSeqResGroups = (ArrayList<Group>)((ArrayList<Group>)seqResGroups).clone();

		boolean startAt1 = true;

		for ( int atomResPos = 0 ; atomResPos < atmResGroups.size() ; atomResPos++){

			//			if ( DEBUG)
			//				System.err.println(" trying to simple match " + atomResPos +"/" + atmResGroups.size());

			// let's try to match this case
			Group atomResGroup = atmResGroups.get(atomResPos);

			// let's ignore waters
			String threeLetterCode = atomResGroup.getPDBName();
			if ( excludeTypes.contains(threeLetterCode)){
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
					
					logger.debug("SEQRES position 1  ({}) does not match ATOM PDB res num 0 ({})",
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
			logger.debug("merging {} {}", seqResPos, atomResGroup);
			newSeqResGroups.set(seqResPos, atomResGroup);

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
			Group g = (Group) groups.get(i);

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
				if (  excludeTypes.contains(g.getPDBName()))
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

					// TODO for some unusual cases the het residue can map to 2 standard aas and thus give an 
					//      insertion of length of 2. How can we handle that? (JD - 20/08/2014)
					//      e.g. SUI maps to DG (in 1oew,A)
					positionIndex.put(seqIndex,i);
					sequence.append(c);
					seqIndex++;
				} else {
					//System.out.println(cc);
					continue;
				}


				//sequence.append("X");
			}

		}

		return sequence.toString();

	}


	private boolean alignNucleotideGroups(List<Group> seqRes, List<Group> atomRes) throws StructureException{

		Map<Integer,Integer> seqresIndexPosition = new HashMap<Integer, Integer>();
		Map<Integer,Integer> atomIndexPosition   = new HashMap<Integer, Integer>();

		String seq1 = getFullAtomSequence(seqRes, seqresIndexPosition);
		//
		String seq2 = getFullAtomSequence(atomRes, atomIndexPosition);

		

		logger.debug("align seq1 ("+ seq1.length()+") " + seq1);
		logger.debug("align seq2 ("+ seq2.length()+") " + seq2);
		


		@SuppressWarnings("rawtypes")
		Sequence s1;
		@SuppressWarnings("rawtypes")
		Sequence s2;

		try {
			s1 = new DNASequence(seq1,AmbiguityDNACompoundSet.getDNACompoundSet());
		} catch (CompoundNotFoundException e){
			///oops perhaps a RNA sequence?
			try {
				s1 = new RNASequence(seq1,AmbiguityRNACompoundSet.getRNACompoundSet());
			} catch (CompoundNotFoundException ex) {
				logger.warn("Could not determine compound set for sequence 1 " + seq1);
				return false;
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
				return false;
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
					return false;
				}
			}
			
			if( ! s2.getCompoundSet().equals(AmbiguityRNACompoundSet.getRNACompoundSet())) {
					try {
						s2 = new RNASequence(seq2,AmbiguityRNACompoundSet.getRNACompoundSet());
					} catch (CompoundNotFoundException ex) {
						logger.warn("Could not align DNA and RNA compound sets: " + seq2);
						return false;
					}
				}
		}
		
		
		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_2();

		GapPenalty penalty = new SimpleGapPenalty();

		short gop = 8;
		short extend = 1;
		penalty.setOpenPenalty(gop);
		penalty.setExtensionPenalty(extend);

		@SuppressWarnings({ "unchecked", "rawtypes" })
		PairwiseSequenceAligner smithWaterman =
				Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.LOCAL, penalty, matrix);

		@SuppressWarnings("rawtypes")
		SequencePair pair = smithWaterman.getPair();



		if ( pair == null)
			throw new StructureException("could not align objects!");


		
		logger.debug(pair.toString(60));
		

		@SuppressWarnings("unchecked")
		boolean noMatchFound = mapDNAChains(seqRes,atomRes,pair,seqresIndexPosition, atomIndexPosition );

		return noMatchFound;

	}


	/** aligns two chains of groups, where the first parent is representing the
	 * list of amino acids as obtained from the SEQRES records, and the second parent
	 * represents the groups obtained from the ATOM records (and containing the actual ATOM information).
	 * This does the actual alignment and if a group can be mapped to a position in the SEQRES then the corresponding
	 * position is replaced with the group that contains the atoms.
	 *
	 * @param seqRes
	 * @param atomRes
	 * @return true if no match has been found
	 * @throws StructureException
	 */
	public boolean align(List<Group> seqRes, List<Group> atomRes) throws StructureException{

		Map<Integer,Integer> seqresIndexPosition = new HashMap<Integer, Integer>();
		Map<Integer,Integer> atomIndexPosition   = new HashMap<Integer, Integer>();

		String seq1 = getFullAtomSequence(seqRes, seqresIndexPosition);
		//
		String seq2 = getFullAtomSequence(atomRes, atomIndexPosition);

		
		logger.debug("align seq1 ("+ seq1.length()+") " + seq1);
		logger.debug("align seq2 ("+ seq2.length()+") " + seq2);
		
		ProteinSequence s1;
		ProteinSequence s2;
		try {
			s1 = new ProteinSequence(seq1);
			s2 = new ProteinSequence(seq2);
		} catch (CompoundNotFoundException e) {
			throw new StructureException("Could not align atom to seqres sequences. Either sequence contains an unknown compound: "+e.getMessage());
		}
		

		SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum65();

		GapPenalty penalty = new SimpleGapPenalty();

		short gop = 8;
		short extend = 1;
		penalty.setOpenPenalty(gop);
		penalty.setExtensionPenalty(extend);


		PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> smithWaterman =
				Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.LOCAL, penalty, matrix);

		SequencePair<ProteinSequence, AminoAcidCompound> pair = smithWaterman.getPair();



		if ( pair == null)
			throw new StructureException("Could not align objects!");


		logger.debug(pair.toString(60));
		

		boolean noMatchFound = mapChains(seqRes,atomRes,pair,seqresIndexPosition, atomIndexPosition );

		return noMatchFound;


	}


	private boolean mapChains(List<Group> seqResGroups, List<Group> atomRes,
			SequencePair<ProteinSequence, AminoAcidCompound> pair,
			Map<Integer,Integer> seqresIndexPosition,
			Map<Integer,Integer> atomIndexPosition
			) throws StructureException{



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
						System.err.println("can't map " + i + ":" + s + " " + posSeq +" " + s1 + " atom: " + posAtom + " " + a1 );
						continue mainLoop;
					}

					// need to trim the names to allow matching e.g in
					// pdb1b2m
					String pdbNameS = s1.getPDBName();
					String pdbNameA = a1.getPDBName();
					if ( pdbNameS == null || pdbNameA == null ){
						System.err.println("null value found at " + posSeq + " when trying to align " + s1 + " and " + a1 + " " + posAtom);
						throw new StructureException("null value found at group.getPDBName()");
					}
					if (! pdbNameA.equals(pdbNameS)){
						if ( ! pdbNameA.trim().equals(pdbNameS.trim())) {
							System.err.println(s1 + " " + posSeq + " does not align with " + a1+ " " + posAtom + " should be: " + s + " : " + a);
							if ( s1.getType().equals(HetatomImpl.type) && a1.getType().equals(HetatomImpl.type)){
								System.err.println("they seem to be hetatoms, so ignoring mismatch.");
							}
							else {
								//  System.exit(0);// for debug only
								//System.out.println(lst1.seqString());
								//System.out.println(lst2.seqString());
								System.err.println("could not match residues " + s1 + " " + a1);
								//throw new StructureException("could not match residues " + s1 + " " + a1);
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
			SequencePair<DNASequence, NucleotideCompound> pair,
			Map<Integer,Integer> seqresIndexPosition,
			Map<Integer,Integer> atomIndexPosition
			) throws StructureException{



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
						System.err.println("can't map " + i + ":" + s + " " + posSeq +" " + s1 + " atom: " + posAtom + " " + a1 );
						continue mainLoop;
					}

					// need to trim the names to allow matching e.g in
					// pdb1b2m
					String pdbNameS = s1.getPDBName();
					String pdbNameA = a1.getPDBName();
					if ( pdbNameS == null || pdbNameA == null ){
						System.err.println("null value found at " + posSeq + " when trying to align " + s1 + " and " + a1 + " " + posAtom);
						throw new StructureException("null value found at group.getPDBName()");
					}
					if (! pdbNameA.equals(pdbNameS)){
						if ( ! pdbNameA.trim().equals(pdbNameS.trim())) {
							System.err.println(s1 + " " + posSeq + " does not align with " + a1+ " " + posAtom + " should be: " + s + " : " + a);
							if ( s1.getType().equals(HetatomImpl.type) && a1.getType().equals(HetatomImpl.type)){
								System.err.println("they seem to be hetatoms, so ignoring mismatch.");
							}
							else {
								//  System.exit(0);// for debug only
								//System.out.println(lst1.seqString());
								//System.out.println(lst2.seqString());
								System.err.println("could not match residues " + s1 + " " + a1);
								//throw new StructureException("could not match residues " + s1 + " " + a1);
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
