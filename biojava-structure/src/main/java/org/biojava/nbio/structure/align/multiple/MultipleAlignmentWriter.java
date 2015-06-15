package org.biojava.nbio.structure.align.multiple;

import java.util.List;

/**
 * This class contains functions for the conversion of {@link MultipleAlignment} to various String outputs.
 * 
 * @author Aleix Lafita
 *
 */
public class MultipleAlignmentWriter {

	/**
	 * Converts the {@link MultipleAlignment} into a multiple sequence alignment String in FASTA format.
	 * @param alignment MultipleAlignment
	 * @return String multiple sequence alignment in FASTA format
	 * @throws StructureAlignmentException if the Atoms cannot be obtained
	 */
	public static String toFASTA(MultipleAlignment alignment) throws StructureAlignmentException {
		
		//Get the alignment sequences
		List<String> alnSequences = MultipleAlignmentTools.getSequenceAlignment(alignment);
		
		String fasta = "";
		for (int st=0; st<alignment.size(); st++){
			//Add the structure identifier as the head of the FASTA
			fasta += ">"+alignment.getEnsemble().getStructureNames().get(st)+"\n"+
					alnSequences.get(st)+"\n";
		}	
		return fasta;
	}
}