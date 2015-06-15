package org.biojava.nbio.structure.align.multiple;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureTools;

/**
 * Utility functions for working with {@link MultipleAlignment}. External functions related to sequence
 * alignment calculations, ... (TODO add methods to this list when implemented)
 * 
 * @author Spencer Bliven
 * @author Aleix Lafita
 *
 */
public class MultipleAlignmentTools {
	
	/**
	 * Calculate the sequence alignment Strings for all blocks in an alignment.
	 * <p>
	 * Blocks are concatenated in the order returned by {@link MultipleAlignment#getBlocks()},
	 * so sequences may not be sequential. Gaps between blocks are omitted,
	 * while gaps within blocks are represented by '-'. Separation between different Blocks is
	 * indicated by a gap in all positions, meaning that there is something unaligned inbetween.
	 * <p>
	 * It is possible to generate a mapping from the sequence alignment to the aligned positions
	 * in the structure alignment. The positions not aligned have the index -1.
	 * @param alignment input MultipleAlignment
	 * @param mapSeqToStruct provides a link from the sequence alignment position to the structure alignment 
	 * 		  position. Specially designed for the GUI. Has to be initialized previously, will be overwritten.
	 * @return a string for each row in the alignment, giving the 1-letter code 
	 *  		for each aligned residue.
	 * @throws StructureAlignmentException if the Atoms cannot be obtained
	 */
	public static List<String> getSequenceAlignment(MultipleAlignment alignment, List<Integer> mapSeqToStruct) throws StructureAlignmentException {

		//Initialize sequence variables
		List<String> alnSequences = new ArrayList<String>();
		for (int str=0; str<alignment.size(); str++) alnSequences.add("");
		mapSeqToStruct.clear();
		List<Atom[]> atoms = alignment.getEnsemble().getAtomArrays();
		int globalPos = -1;
		
		//Loop through all the alignment Blocks in the order given
		for (int b=0; b<alignment.getBlocks().size(); b++){
			if (b!=0){
				//Add a gap to all structures in order to separate visually the blocks in the alignment
				for (int str=0; str<alignment.size(); str++) alnSequences.set(str,alnSequences.get(str).concat("-"));
				mapSeqToStruct.add(-1); //means no aligned position
			}
			//Store the previous position added to the sequence alignment for this structure
			int[] previousPos = new int[alignment.size()];
			Arrays.fill(previousPos, -1);
			//Store provisional characters
			char[] provisionalChar = new char[alignment.size()];
			Arrays.fill(provisionalChar, '-');
			//Loop through all the alignment positions in the Block
			for (int pos=0; pos<alignment.getBlocks().get(b).length(); pos++){
				globalPos++;
				boolean gaps = true;  //If any structure is not consecutive with the previousPos
				//While the next position cannot be considered because there are still non consecutive residues
				while (gaps){
					gaps = false;
					//Loop through all the structures
					for (int str=0; str<alignment.size(); str++){
						//If it is the first position or before it was null
						if (previousPos[str] == -1){
							Integer residue = alignment.getBlocks().get(b).getAlignRes().get(str).get(pos);
							if (residue == null) provisionalChar[str] = '-';
							else provisionalChar[str] = StructureTools.get1LetterCode(atoms.get(str)[residue].getGroup().getPDBName());		
						}
						else{
							Integer residue = alignment.getBlocks().get(b).getAlignRes().get(str).get(pos);
							if (residue == null) provisionalChar[str] = '-';
							else if (previousPos[str]+1 == residue){
								provisionalChar[str] = StructureTools.get1LetterCode(atoms.get(str)[residue].getGroup().getPDBName());
							}
							else{
								provisionalChar[str] = ' ';  //This means there is a spacing (non-consecutive)
								gaps = true;
							}
						}
					}//End all structure analysis
					if (gaps){
						for (int str=0; str<alignment.size(); str++){
							if (provisionalChar[str] == ' ') {  //It means this residue was the non-consecutive one
								alnSequences.set(str,alnSequences.get(str).concat(""+StructureTools.get1LetterCode(atoms.get(str)[previousPos[str]+1].getGroup().getPDBName())));
								previousPos[str] += 1;
							} else { //insert a gap otherwise and do not change the index
								alnSequences.set(str,alnSequences.get(str).concat("-"));
							}
						}
						mapSeqToStruct.add(-1); //meaning that this is an unaligned position
					} 
					else {  //Append the provisional and update the indices otherwise
						for (int str=0; str<alignment.size(); str++){
							alnSequences.set(str,alnSequences.get(str).concat(""+provisionalChar[str]));
							if (provisionalChar[str] != '-') previousPos[str] = alignment.getBlocks().get(b).getAlignRes().get(str).get(pos);
						}
						mapSeqToStruct.add(globalPos);
					}
				}
			}
		}
		return alnSequences;
	}
	
	/**
	 * Calculate the sequence alignment Strings for all blocks in an alignment.
	 * <p>
	 * Blocks are concatenated in the order returned by {@link MultipleAlignment#getBlocks()},
	 * so sequences may not be sequential. Gaps between blocks are omitted,
	 * while gaps within blocks are represented by '-'. Separation between different Blocks is
	 * indicated by a gap in all positions, meaning that there is something unaligned inbetween.
	 * @param alignment input MultipleAlignment
	 * @return a string for each row in the alignment, giving the 1-letter code 
	 *  		for each aligned residue.
	 * @throws StructureAlignmentException if the Atoms cannot be obtained
	 */
	public static List<String> getSequenceAlignment(MultipleAlignment alignment) throws StructureAlignmentException {
		return getSequenceAlignment(alignment, new ArrayList<Integer>());
	}
}