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
package org.biojava3.proteinproperties;

import java.util.Hashtable;

import org.biojava.bio.structure.AminoAcid;
import org.biojava3.core.sequence.ProteinSequence;

/**
 * An interface to generate some basic physico-chemical properties of protein sequences
 * 
 * AHFU - TODO
 * What constraints are we talking about? Isnt it the same as input restrictions?
 * What input restrictions?
 * Netcharge and isoelectric points approach I found are acceptable?
 * getAAComposition - What do Peter mean by use interface as return value?
 * 
 * @author kohchuanhock
 * @version 2011.05.09
 */
public interface BasicProperties {
	/**
	 * Returns the molecular weight of sequence. 
	 * The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * <p>
	 * This method will sum the molecular weight of each amino acid in the sequence.
	 * Molecular weights are based on (http://au.expasy.org/tools/findmod/findmod_masses.html#AA).
	 *
	 * @param  sequence a protein sequence consisting of non-ambiguous characters only
	 * @return the total molecular weight of sequence
	 * @see ProteinSequence
	 */
	public double getMolecularWeight(ProteinSequence sequence);

	/**
	 * Returns the extinction coefficient of sequence. 
	 * The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * <p>
	 * The extinction coefficient indicates how much light a protein absorbs at a certain wavelength. 
	 * It is useful to have an estimation of this coefficient for following a protein which a spectrophotometer when purifying it.
	 * The computation of extinction coefficient follows the documentation in {@link http://au.expasy.org/tools/protparam-doc.html}. 
	 *
	 * @param sequence a protein sequence consisting of non-ambiguous characters only
	 * @param assumeCysReduced true if Cys are assumed to be reduced and false if Cys are assumed to form cystines
	 * @return the extinction coefficient of sequence
	 * @see ProteinSequence
	 */
	public double getExtinctionCoefficient(ProteinSequence sequence,
			boolean assumeCysReduced);
	
	/**
	 * Returns the instability index of sequence. 
	 * The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * <p>
	 * The instability index provides an estimate of the stability of your protein in a test tube.
	 * The computation of instability index follows the documentation in {@link http://au.expasy.org/tools/protparam-doc.html}. 
	 *
	 * @param sequence a protein sequence consisting of non-ambiguous characters only
	 * @return the instability index of sequence
	 * @see ProteinSequence
	 */
	public double getInstabilityIndex(ProteinSequence sequence);

	/**
	 * Returns the apliphatic index of sequence. 
	 * The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * <p>
	 * The aliphatic index of a protein is defined as the relative volume occupied by aliphatic side chains (alanine, valine, isoleucine, and leucine). 
	 * It may be regarded as a positive factor for the increase of thermostability of globular proteins.
	 * The computation of aliphatic index follows the documentation in {@link http://au.expasy.org/tools/protparam-doc.html}.
	 * 
	 * @param sequence a protein sequence consisting of non-ambiguous characters only
	 * @return the aliphatic index of sequence
	 * @see ProteinSequence
	 */
	public double getApliphaticIndex(ProteinSequence sequence);

	/**
	 * Returns the average hydropathy value of sequence. 
	 * The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * <p>
	 * The average value for a sequence is calculated as the sum of hydropathy values of all the amino acids, divided by the number of residues in the sequence. 
	 * Hydropathy values are based on 
	 * (Kyte, J. and Doolittle, R.F. (1982) A simple method for displaying the hydropathic character of a protein. J. Mol. Biol. 157, 105-132).
	 * 
	 * @param sequence a protein sequence consisting of non-ambiguous characters only
	 * @return the average hydropathy value of sequence
	 * @see ProteinSequence
	 */
	public double getAvgHydropathy(ProteinSequence sequence);

	/**
	 * Returns the isoelectric point of sequence. 
	 * The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * <p>
	 * The isoelectric point is the pH at which the protein carries no net electrical charge.
	 * The isoelectric point will be computed based on approach stated in 
	 * (http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator-notes.asp#PI)
	 * 
	 * @param sequence a protein sequence consisting of non-ambiguous characters only
	 * @return the isoelectric point of sequence
	 * @see ProteinSequence
	 */
	public double getIsoelectricPoint(ProteinSequence sequence);

	/**
	 * Returns the net charge of sequence at pH 7. 
	 * The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * <p>
	 * The net charge will be computed using the approach stated in 
	 * (http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator-notes.asp#NetCharge)
	 * 
	 * @param sequence a protein sequence consisting of non-ambiguous characters only
	 * @return the net charge of sequence at pH 7
	 * @see ProteinSequence
	 */
	public double getNetCharge(ProteinSequence sequence);

	/**
	 * Returns the composition of specified amino acid in the sequence. 
	 * The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * The aminoAcidCode must be a non-ambiguous character.
	 * <p>
	 * The composition of an amino acid is the total number of its occurrence, divided by the total length of the sequence.
	 * 
	 * @param sequence a protein sequence consisting of non-ambiguous characters only
	 * @param aminoAcidCode the code of the amino acid to compute
	 * @return the composition of specified amino acid in the sequence
	 * @see ProteinSequence, AminoAcid
	 */
	public double getEnrichment(ProteinSequence sequence, AminoAcid aminoAcidCode);

	/**
	 * Returns the composition of the 20 standard amino acid in the sequence. 
	 * The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * <p>
	 * The composition of an amino acid is the total number of its occurrence, divided by the total length of the sequence.
	 * 
	 * TODO use interfaces as return type.
	 * 
	 * TODO same as above RE amino acid code
	 * 
	 * 
	 * @param sequence a protein sequence consisting of non-ambiguous characters only
	 * @return the composition of the 20 standard amino acid in the sequence
	 * @see ProteinSequence, AminoAcid
	 */
	public Hashtable<Character, Double> getAAComposition(ProteinSequence sequence);
}
