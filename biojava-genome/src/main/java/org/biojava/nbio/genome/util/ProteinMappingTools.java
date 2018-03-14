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
