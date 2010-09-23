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
import java.util.Iterator;
import java.util.List;


import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.HetatomImpl;
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
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.template.Compound;




/** Aligns the SEQRES residues to the ATOM residues.
 * The AminoAcids that can be matched between the two of them will be set in the SEQRES
 * chains
 *
 *
 * @author Andreas Prlic
 *
 */
public class SeqRes2AtomAligner {

	boolean DEBUG = false;

	static final List<String> excludeTypes;
	

	String alignmentString;
	static {
		excludeTypes = new ArrayList<String>();
		excludeTypes.add("HOH"); // we don't want to align water against the SEQRES...
		excludeTypes.add("DOD"); // deuterated water


		//matrix.printMatrix();

	}



	public SeqRes2AtomAligner(){
		alignmentString = "";
	}

	public String getAlignmentString() {
		return alignmentString;
	}

	public boolean isDEBUG() {
		return DEBUG;
	}

	public void setDEBUG(boolean debug) {
		DEBUG = debug;
	}

	private Chain getMatchingAtomRes(Chain seqRes, List<Chain> atomList)
	throws StructureException {
		Iterator<Chain> iter = atomList.iterator();
		while(iter.hasNext()){
			Chain atom = iter.next();
			if ( atom.getName().equals(seqRes.getName())){
				return atom;
			}
		}
		throw new StructureException("could not match seqres chainID >" + seqRes.getName() + "< to ATOM chains!");
	}
	public void align(Structure s, List<Chain> seqResList){

		//List<Chain> seqResList = s.getSeqRes();
		List<Chain> atomList   = s.getModel(0);

		Iterator<Chain> iter = seqResList.iterator();
		// List<Chain> chains = new ArrayList<Chain>();
		while ( iter.hasNext()){
			Chain seqRes = iter.next();

			if ( seqRes.getAtomGroups("amino").size() < 1) {
				if (DEBUG){
					System.out.println("chain " + seqRes.getName() + " does not contain amino acids, ignoring...");
				}
				continue;
			}

			try {

				Chain atomRes = getMatchingAtomRes(seqRes,atomList);
				if ( atomRes.getAtomGroups("amino").size() < 1) {
					if (DEBUG){
						System.out.println("chain " + atomRes.getName() + " does not contain amino acids, ignoring...");
					}
					continue;
				}
				if ( DEBUG )
					System.out.println("Alignment for chain "+ atomRes.getName());

				List<Group> seqResGroups = seqRes.getAtomGroups();
				boolean noMatchFound = align(seqResGroups,atomRes.getAtomGroups());
				if ( ! noMatchFound){
					atomRes.setSeqResGroups(seqResGroups);
				}

				//chains.add(mapped);
			} catch (StructureException e){
				e.printStackTrace();
			}
		}
		//s.setChains(0,chains);



	}


	/** returns the full sequence of the Atom records of a chain
	 * with X instead of HETATMSs. The advantage of this is
	 * that it allows us to also align HETATM groups back to the SEQRES.
	 * @param groups the list of groups in a chain
	 *
	 * @return string representations
	 */
	public String getFullAtomSequence(List<Group> groups){


		StringBuffer sequence = new StringBuffer() ;
		for ( int i=0 ; i< groups.size(); i++){
			Group g = (Group) groups.get(i);
			if ( g instanceof AminoAcid ){
				AminoAcid a = (AminoAcid)g;
				char oneLetter =a.getAminoType();
				if ( oneLetter == '?')
					oneLetter = 'X';
				sequence.append(oneLetter );
			} else {
				if (  excludeTypes.contains(g.getPDBName()))
					continue;

				if ( g.size() == 1 ) {
					try {
						Atom a = g.getAtom(0);
						if (a.getElement().isMetal())
							continue;

						continue;
					} catch (StructureException e){
						e.printStackTrace();
						continue;
					}
				}

				ChemComp cc = g.getChemComp();
				if ( 
						ResidueType.lPeptideLinking.equals(cc.getResidueType()) ||
						PolymerType.PROTEIN_ONLY.contains(cc.getPolymerType())  ||
						PolymerType.POLYNUCLEOTIDE_ONLY.contains(cc.getPolymerType())
				) {
					//System.out.println(cc.getOne_letter_code());
					String c = cc.getOne_letter_code();
					if ( c.equals("?"))
						c = "X";

					sequence.append(c);
				} else {
					//System.out.println(cc);
					continue;
				}


				//sequence.append("X");
			}

		}

		return sequence.toString();

	}

	/** aligns two chains of groups, where the first chain is representing the
	 * list of amino acids as obtained from the SEQRES records, and the second chain
	 * represents the groups obtained from the ATOM records (and containing the actual ATOM information).
	 * This does the actual alignment and if a group can be mapped to a position in the SEQRES then the corresponding
	 * position is repplaced with the group that contains the atoms.
	 *
	 * @param seqRes
	 * @param atomRes
	 * @return true if no match has bee found
	 * @throws StructureException
	 */
	public boolean align(List<Group> seqRes, List<Group> atomRes) throws StructureException{
		/** int MAX_SIZE = 1000;
        if ( (seqRes.size() > MAX_SIZE)
                ||( atomRes.size() > MAX_SIZE) ) {
                    System.err.println("can not align chains, length size exceeds limits!");
                    return false;
                }
		 */
		String seq1 = getFullAtomSequence(seqRes);
		//String seq1 = seqRes.getSeqResSequence();
		String seq2 = getFullAtomSequence(atomRes);

		if ( DEBUG ) {

			System.out.println("align seq1 " + seq1);
			System.out.println("align seq2 " + seq2);
		}
		ProteinSequence s1 = new ProteinSequence(seq1);
		ProteinSequence s2 = new ProteinSequence(seq2);

		SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum65();

		GapPenalty penalty = new SimpleGapPenalty();

		short gop = 8;
		short extend = 1;
		penalty.setOpenPenalty(gop);
		penalty.setExtensionPenalty(extend);

		PairwiseSequenceAligner smithWaterman =
			Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.LOCAL, penalty, matrix);

		SequencePair<ProteinSequence, AminoAcidCompound> pair = smithWaterman.getPair();

		

		if ( pair == null)
			throw new StructureException("could not align objects!");

		
		if ( DEBUG ) {
			System.out.println(pair.toString(60));
		}

		boolean noMatchFound = mapChains(seqRes,atomRes,pair);
		return noMatchFound;

	}


	private boolean mapChains(List<Group> seqRes, List<Group> atomRes,
			 SequencePair<ProteinSequence, AminoAcidCompound> pair
	) throws StructureException{


		// at the present stage the future seqREs are still stored as Atom groups in the seqRes chain...
		List<Group> seqResGroups = seqRes;

		int aligLength = pair.getLength();

		// System.out.println(gapSymbol.getName());

		// make sure we actually find an alignment
		boolean noMatchFound = true;

		Compound gapSymbol =  AminoAcidCompoundSet.getAminoAcidCompoundSet().getCompoundForString("-");

		for (int i = 1; i <= aligLength; i++) {

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

				Group s1 = seqRes.get(posSeq);
				Group a1 = atomRes.get(posAtom);
				//System.out.println(s1.getPDBName() + " == " + a1.getPDBName());
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
				seqResGroups.set(posSeq,a1);
				noMatchFound = false;
			}
		}


		// now we merge the two chains into one
		// the Groups that can be aligned are now pointing to the
		// groups in the Atom records.
		if (  noMatchFound) {

			if ( DEBUG )
				System.out.println("no alignment found!");
		}
		return noMatchFound;

	}



}
