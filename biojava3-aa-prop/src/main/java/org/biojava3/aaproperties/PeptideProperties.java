package org.biojava3.aaproperties;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;

/**
 * This is an adaptor class which enable the ease of generating protein properties.
 * At least one adaptor method is written for each available properties provided in IPeptideProperties. 
 * 
 * @author kohchuanhock
 * @version 2011.05.21
 * @see IPeptideProperties
 * @see PeptidePropertiesImpl
 */
public class PeptideProperties {
	public enum SingleLetterAACode { W, C, M, H, Y, F, Q, N, I, R, D, P, T, K, E, V, S, G, A, L}
	
	public static Set<Character> standardAASet;
	
	static{
		standardAASet = new HashSet<Character>();
		for(SingleLetterAACode c:SingleLetterAACode.values()) standardAASet.add(c.toString().charAt(0));
	}
	
	/**
	 * An adaptor method to return the molecular weight of sequence. 
	 * The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * This method will sum the molecular weight of each amino acid in the
	 * sequence. Molecular weights are based on <a href="http://au.expasy.org/tools/findmod/findmod_masses.html#AA">here</a>.
	 * 
	 * @param sequence
	 * 			a protein sequence consisting of non-ambiguous characters only
	 * @param decimals
	 * 			round off to the number of decimals
	 * @return the total molecular weight of sequence + weight of water molecule
	 */
	public static final double getMolecularWeight(String sequence, int decimals){
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return Utils.roundToDecimals(pp.getMolecularWeight(pSequence), decimals);
	}
	
	public static final double getMolecularWeight(String sequence){
		return getMolecularWeight(sequence, -1);
	}
	
	/**
	 * An adaptor method to returns the absorbance (optical density) of sequence. The sequence argument
	 * must be a protein sequence consisting of only non-ambiguous characters.
	 * The computation of absorbance (optical density) follows the
	 * documentation in <a href="http://au.expasy.org/tools/protparam-doc.html">here</a>.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param assumeCysReduced
	 *            true if Cys are assumed to be reduced and false if Cys are
	 *            assumed to form cystines
	 * @param decimals
	 * 			round off to the number of decimals
	 * @return the absorbance (optical density) of sequence
	 */
	public static final double getAbsorbance(String sequence, boolean assumeCysReduced, int decimals){
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return Utils.roundToDecimals(pp.getAbsorbance(pSequence, assumeCysReduced), decimals);
	}
	
	public static final double getAbsorbance(String sequence, boolean assumeCysReduced){
		return getAbsorbance(sequence, assumeCysReduced, -1);
	}

	/**
	 * An adaptor method to return the extinction coefficient of sequence. The sequence argument
	 * must be a protein sequence consisting of only non-ambiguous characters.
	 * The extinction coefficient indicates how much light a protein absorbs at
	 * a certain wavelength. It is useful to have an estimation of this
	 * coefficient for following a protein which a spectrophotometer when
	 * purifying it. The computation of extinction coefficient follows the
	 * documentation in <a href="http://au.expasy.org/tools/protparam-doc.html">here</a>.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param assumeCysReduced
	 *            true if Cys are assumed to be reduced and false if Cys are
	 *            assumed to form cystines
	 * @return the extinction coefficient of sequence
	 */
	public static final double getExtinctionCoefficient(String sequence, boolean assumeCysReduced) {
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getExtinctionCoefficient(pSequence, assumeCysReduced);
	}

	/**
	 * An adaptor method to return the instability index of sequence. The sequence argument must be
	 * a protein sequence consisting of only non-ambiguous characters.
	 * The instability index provides an estimate of the stability of your
	 * protein in a test tube. The computation of instability index follows the
	 * documentation in <a href="http://au.expasy.org/tools/protparam-doc.html">here</a>.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param decimals
	 * 			round off to the number of decimals
	 * @return the instability index of sequence
	 */
	public static final double getInstabilityIndex(String sequence, int decimals) {
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return Utils.roundToDecimals(pp.getInstabilityIndex(pSequence), decimals);
	}
	
	public static final double getInstabilityIndex(String sequence){
		return getInstabilityIndex(sequence, -1);
	}

	/**
	 * An adaptor method to return the apliphatic index of sequence. The sequence argument must be a
	 * protein sequence consisting of only non-ambiguous characters.
	 * The aliphatic index of a protein is defined as the relative volume
	 * occupied by aliphatic side chains (alanine, valine, isoleucine, and
	 * leucine). It may be regarded as a positive factor for the increase of
	 * thermostability of globular proteins. The computation of aliphatic index
	 * follows the documentation in <a href="http://au.expasy.org/tools/protparam-doc.html">here</a>.
	 * A protein whose instability index is smaller than 40 is predicted as stable, a value above 40 predicts that the protein may be unstable.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param decimals
	 * 			round off to the number of decimals
	 * @return the aliphatic index of sequence
	 */
	public static final double getApliphaticIndex(String sequence, int decimals) {
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return Utils.roundToDecimals(pp.getApliphaticIndex(pSequence), decimals);
	}
	
	public static final double getApliphaticIndex(String sequence){
		return getApliphaticIndex(sequence, -1);
	}
	
	/**
	 * An adaptor method to return the average hydropathy value of sequence. The sequence argument
	 * must be a protein sequence consisting of only non-ambiguous characters.
	 * The average value for a sequence is calculated as the sum of hydropathy
	 * values of all the amino acids, divided by the number of residues in the
	 * sequence. Hydropathy values are based on (Kyte, J. and Doolittle, R.F.
	 * (1982) A simple method for displaying the hydropathic character of a
	 * protein. J. Mol. Biol. 157, 105-132).
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param decimals
	 * 			round off to the number of decimals        
	 * @return the average hydropathy value of sequence
	 */
	public static final double getAvgHydropathy(String sequence, int decimals) {
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return Utils.roundToDecimals(pp.getAvgHydropathy(pSequence), decimals);
	}
	
	public static final double getAvgHydropathy(String sequence){
		return getAvgHydropathy(sequence, -1);
	}

	/**
	 * An adaptor method to return the isoelectric point of sequence. The sequence argument must be
	 * a protein sequence consisting of only non-ambiguous characters.
	 * The isoelectric point is the pH at which the protein carries no net
	 * electrical charge. The isoelectric point will be computed based on
	 * approach stated in 
	 * <a href="http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator-notes.asp#PI">here</a>
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param decimals
	 * 			round off to the number of decimals
	 * @return the isoelectric point of sequence
	 */
	public static final double getIsoelectricPoint(String sequence, int decimals) {
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return Utils.roundToDecimals(pp.getIsoelectricPoint(pSequence), decimals);
	}
	
	public static final double getIsoelectricPoint(String sequence){
		return getIsoelectricPoint(sequence, -1);
	}

	/**
	 * An adaptor method to return the net charge of sequence at pH 7. The sequence argument must be
	 * a protein sequence consisting of only non-ambiguous characters.
	 * The net charge will be computed using the approach stated in 
	 * <a href="http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator-notes.asp#NetCharge>here</a>
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param decimals
	 * 			round off to the number of decimals
	 * @return the net charge of sequence at pH 7
	 */
	public static final double getNetCharge(String sequence, int decimals) {
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return Utils.roundToDecimals(pp.getNetCharge(pSequence), decimals);
	}
	
	public static final double getNetCharge(String sequence){
		return getNetCharge(sequence, -1);
	}
	
	/**
	 * An adaptor method to return the composition of specified amino acid in the sequence. The
	 * sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters. The aminoAcidCode must be a non-ambiguous
	 * character.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param aminoAcidCode
	 *            the code of the amino acid to compute
	 * @return the composition of specified amino acid in the sequence
	 * @see SingleLetterAACode
	 */
	public static final double getEnrichment(String sequence, SingleLetterAACode aminoAcidCode) {
		return getEnrichment(sequence, aminoAcidCode.toString());
	}
	
	/**
	 * An adaptor method to return the composition of specified amino acid in the sequence. The
	 * sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters. The aminoAcidCode must be a non-ambiguous
	 * character.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param aminoAcidCode
	 *            the code of the amino acid to compute
	 * @return the composition of specified amino acid in the sequence
	 */
	public static final double getEnrichment(String sequence, char aminoAcidCode){
		return getEnrichment(sequence, aminoAcidCode + "");
	}
	
	/**
	 * An adaptor method to return the composition of specified amino acid in the sequence. The
	 * sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters. The aminoAcidCode must be a non-ambiguous
	 * character.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @param aminoAcidCode
	 *            the code of the amino acid to compute
	 * @return the composition of specified amino acid in the sequence
	 */
	public static final double getEnrichment(String sequence, String aminoAcidCode){
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		return pp.getEnrichment(pSequence, aaSet.getCompoundForString(aminoAcidCode));
	}

	/**
	 * An adaptor method to return the composition of the 20 standard amino acid in the sequence.
	 * The sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @return the composition of the 20 standard amino acid in the sequence
	 * @see AminoAcidCompound
	 */
	public static final Map<AminoAcidCompound, Double> getAAComposition(String sequence) {
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getAAComposition(pSequence);
	}
	
	/**
	 * An adaptor method to return the composition of the 20 standard amino acid in the sequence.
	 * The sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @return the composition of the 20 standard amino acid in the sequence
	 */
	public static final Map<String, Double> getAACompositionString(String sequence){
		Map<AminoAcidCompound, Double> aa2Composition = getAAComposition(sequence);
		Map<String, Double> aaString2Composition = new HashMap<String, Double>();
		for(AminoAcidCompound aaCompound:aa2Composition.keySet()){
			aaString2Composition.put(aaCompound.getShortName(), aa2Composition.get(aaCompound));
		}
		return aaString2Composition;
	}
	
	/**
	 * An adaptor method to return the composition of the 20 standard amino acid in the sequence.
	 * The sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 *            a protein sequence consisting of non-ambiguous characters only
	 * @return the composition of the 20 standard amino acid in the sequence
	 */
	public static final Map<Character, Double> getAACompositionChar(String sequence){
		Map<AminoAcidCompound, Double> aa2Composition = getAAComposition(sequence);
		Map<Character, Double> aaChar2Composition = new HashMap<Character, Double>();
		for(AminoAcidCompound aaCompound:aa2Composition.keySet()){
			aaChar2Composition.put(aaCompound.getShortName().charAt(0), aa2Composition.get(aaCompound));
		}
		return aaChar2Composition;
	}
}	
