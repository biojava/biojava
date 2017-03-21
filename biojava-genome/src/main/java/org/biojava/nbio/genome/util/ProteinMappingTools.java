package org.biojava.nbio.genome.util;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.RNASequence;

public class ProteinMappingTools {

    /** Converts the DNA sequence to protein sequence.
    *
    * @param dnaSequence the DNA sequence
    * 
    * @return the protein sequence
    */
	public String convertDNAtoProteinSequence(String dnaSequence) throws CompoundNotFoundException {
		DNASequence dna = new DNASequence(dnaSequence);
		RNASequence mRNA = dna.getRNASequence();
		return mRNA.getProteinSequence().toString();
	}
}
