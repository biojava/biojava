package org.biojava.nbio.genome.util;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;

public class ProteinMappingTools {

    /** Converts the DNA sequence to protein sequence.
    *
    * @param dnaSequence the DNA sequence
    * 
    * @return the protein sequence
    */
	public static ProteinSequence convertDNAtoProteinSequence(String dnaSequence) throws CompoundNotFoundException {
		DNASequence dna = new DNASequence(dnaSequence);
		return convertDNAtoProteinSequence(dna);
	}

	/** Converts the DNA sequence to protein sequence.
	 *
	 * @param dnaSequence the DNA sequence
	 *
	 * @return the protein sequence
	 */
	public static ProteinSequence convertDNAtoProteinSequence(DNASequence dnaSequence) throws CompoundNotFoundException {
		RNASequence mRNA = dnaSequence.getRNASequence();
		return mRNA.getProteinSequence();
	}
}
