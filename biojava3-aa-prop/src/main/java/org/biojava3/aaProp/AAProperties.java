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
 * Created on 2011.05.09 by kohchuanhock
 *
 */
package org.biojava3.aaProp;

import java.util.Hashtable;

import org.biojava3.core.sequence.ProteinSequence;

/**
 * An interface to generate some basic physico-chemical properties of protein
 * sequences
 * 
 * TODO: Keep the capitalization consistent (use Index instead of index)
 * 
 * TODO The javadoc for each method should contain a short description of the
 * method. The comments should make clear what the method is doing, for example
 * what is the instability index? Is there only one way to calculate it? If so
 * it should be described here.
 * 
 * TODO read on javadoc standards for example -
 * http://www.oracle.com/technetwork/java/javase/documentation/index-137868.html
 * your javadoc is not complete
 * 
 * TODO Have you thought about the input restrictions? For example do you allow
 * ambiguous characters in the input sequence? This should be documented.
 * 
 * TODO read on java naming conventions for example here
 * http://java.about.com/od/javasyntax/a/nameconventions.htm and stick to them!
 * Start from renaming the package for this interface.
 * 
 * Once you acted on todo comment, feel free to delete it.
 * 
 * @author kohchuanhock
 * @version 2011.05.09
 */
public interface AAProperties {

	/**
	 * @param sequence
	 * @return Molecular weight of sequence.
	 */
	public double getMolecularWeight(ProteinSequence sequence);

	/**
	 * @param sequence
	 * @param assumeCysReduced
	 * @return Extinction coefficient of sequence. Two possible assumption 1)
	 *         Cys are reduced 2) Cys form cystines.
	 */
	public double getExtinctionCoefficient(ProteinSequence sequence,
			boolean assumeCysReduced);

	/**
	 * @param sequence
	 * @return Instability index of sequence.
	 */
	public double getInstabilityIndex(ProteinSequence sequence);

	/**
	 * TODO name capitalization
	 * 
	 * @param sequence
	 * @return Aliphatic Index of sequence.
	 */
	public double getApliphaticindex(ProteinSequence sequence);

	/**
	 * @param sequence
	 * @return Grand Average of Hydropathy of sequence.
	 */
	public double getAvgHydropathy(ProteinSequence sequence);

	/**
	 * @param sequence
	 * @return Isoelectric point of sequence.
	 */
	public double getIsoPoint(ProteinSequence sequence);

	/**
	 * TODO Is it any different from ProteinSequence.getLenght()?
	 * 
	 * @param sequence
	 * @return Length of sequence.
	 */
	public double getLength(ProteinSequence sequence);

	/**
	 * @param sequence
	 * @return Net charge of sequence
	 */
	public double getNetCharge(ProteinSequence sequence);

	/**
	 * TODO How about using AminoAcid class instead of char? Just make sure
	 * there is a convenient way for users to define (1) single AA e.g. M, (2)
	 * subset of AA, e.g. MTA (Met-Ala-Trp). May be defining enum for letter
	 * codes AA would be good?
	 * 
	 * 
	 * 
	 * @param sequence
	 * @param aminoAcidCode
	 * @return Composition of aminoAcidCode in sequence. Total number of
	 *         aminoAcidCode / Length of sequence.
	 */
	public double getEnrichment(ProteinSequence sequence, char aminoAcidCode);

	/**
	 * TODO use interfaces as return type.
	 * 
	 * TODO FYI: do not use Hashtable - this is old badly performing class, use
	 * HashMap instead
	 * 
	 * TODO same as above RE amino acid code
	 * 
	 * @param sequence
	 * @return Hashtable that contains the composition of all amino acids.
	 */
	public Hashtable<Character, Double> getAAComposition(
			ProteinSequence sequence);
}
