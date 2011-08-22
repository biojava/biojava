package org.biojava3.aaproperties.profeat;

import java.util.Map;

import org.biojava3.core.sequence.ProteinSequence;

public interface IProfeatProperties {
	/**
	 * Based on Table 2 of http://nar.oxfordjournals.org/content/34/suppl_2/W32.full.pdf<br/>
	 * An interface class to generate the properties of a protein sequence based on its converted attributes.<br/>
	 * The seven different attributes are<p/>
	 * Hydrophobicity (Polar, Neutral, Hydrophobicity)<br/>
	 * Normalized van der Waals volume (Range 0 - 2.78, 2.95 - 4.0, 4.03 - 8.08)<br/>
	 * Polarity (Value 4.9 - 6.2, 8.0 - 9.2, 10.4 - 13.0)<br/>
	 * Polarizability (Value 0 - 1.08, 0.128 - 0.186, 0.219 - 0.409)<br/>
	 * Charge (Positive, Neutral, Negative)<br/>
	 * Secondary structure (Helix, Strand, Coil)<br/>
	 * Solvent accessibility (Buried, Exposed, Intermediate)<br/>
	 * 
	 * @author kohchuanhock
	 * @version 2011.06.16
	 * @since 3.0.2
	 */
	
	/**
	 * Enumeration of the seven different attributes
	 */
	public enum ATTRIBUTE {HYDROPHOBICITY, VOLUME, POLARITY, POLARIZABILITY, CHARGE, SECONDARYSTRUCTURE, SOLVENTACCESSIBILITY};
	/**
	 * Enumeration of the three different groupings for each attributes 
	 */
	public enum GROUPING {GROUP1, GROUP2, GROUP3};
	/**
	 * Enumeration of the transition between groupA and groupB 
	 */
	public enum TRANSITION {BETWEEN_11, BETWEEN_22, BETWEEN_33, BETWEEN_12, BETWEEN_13, BETWEEN_23};
	/**
	 * Enumeration of the distribution for the first, first 25%, first 50%, first 75% and 100% of the grouping
	 */
	public enum DISTRIBUTION {FIRST, FIRST25, FIRST50, FIRST75, ALL};

	/**
	 * Returns the composition of the specific grouping for the given attribute.
	 * 
	 * @param sequence
	 * 	a protein sequence consisting of non-ambiguous characters only
	 * @param attribute
	 * 	one of the seven attributes (Hydrophobicity, Volume, Polarity, Polarizability, Charge, SecondaryStructure or SolventAccessibility)
	 * @param group
	 * 	the grouping to be computed
	 * @return
	 * 	returns the composition of the specific grouping for the given attribute
	 * @throws Exception
	 * 	throws Exception if attribute or group are unknown
	 */
	public double getComposition(ProteinSequence sequence, ATTRIBUTE attribute, GROUPING group) throws Exception;
	
	public Map<GROUPING, Double> getComposition(ProteinSequence sequence, ATTRIBUTE attribute) throws Exception;
	
	public Map<ATTRIBUTE, Map<GROUPING, Double>> getComposition(ProteinSequence sequence) throws Exception;
	
	/**
	 * Returns the number of transition between the specified groups for the given attribute with respect to the length of sequence.
	 * 
	 * @param sequence
	 * 	a protein sequence consisting of non-ambiguous characters only
	 * @param attribute
	 * 	one of the seven attributes (Hydrophobicity, Volume, Polarity, Polarizability, Charge, SecondaryStructure or SolventAccessibility)
	 * @param transition
	 * 	the interested transition between the groups
	 * @return
	 *  returns the number of transition between the specified groups for the given attribute with respect to the length of sequence.
	 * @throws Exception
	 * 	throws Exception if attribute or group are unknown
	 */
	public double getTransition(ProteinSequence sequence, ATTRIBUTE attribute, TRANSITION transition) throws Exception;
	
	public Map<TRANSITION, Double> getTransition(ProteinSequence sequence, ATTRIBUTE attribute) throws Exception;
	
	public Map<ATTRIBUTE, Map<TRANSITION, Double>> getTransition(ProteinSequence sequence) throws Exception;
	
	/**
	 * Computes and return the position with respect to the sequence where the given distribution of the grouping can be found.<br/>
	 * Example: "1111122222"<br/>
	 * For the above example,<br/>
	 * position of the GROUPING.GROUP1 && DISTRIBUTION.FIRST = 0/10 (because the first occurrence of '1' is at position 0)<br/> 
	 * position of the GROUPING.GROUP1 && DISTRIBUTION.ALL = 4/10 (because all occurrences of '1' happens on and before position 4)<br/>
	 * 
	 * @param sequence
	 * 	a protein sequence consisting of non-ambiguous characters only
	 * @param attribute
	 * 	one of the seven attributes (Hydrophobicity, Volume, Polarity, Polarizability, Charge, SecondaryStructure or SolventAccessibility)
	 * @param group	
	 * 	one the three groups for the attribute
	 * @param distribution
	 * 	the distribution of the grouping
	 * 	
	 * @return
	 * 	the position with respect to the length of sequence where the given distribution of the grouping can be found.<br/> 
	 * @throws Exception
	 * 	throws Exception if attribute or group are unknown
	 */
	public double getDistributionPosition(ProteinSequence sequence, ATTRIBUTE attribute, GROUPING group, DISTRIBUTION distribution) throws Exception;
	
	public Map<DISTRIBUTION, Double> getDistributionPosition(ProteinSequence sequence, ATTRIBUTE attribute, GROUPING group) throws Exception;
	
	public Map<GROUPING, Map<DISTRIBUTION, Double>> getDistributionPosition(ProteinSequence sequence, ATTRIBUTE attribute) throws Exception;
	
	public Map<ATTRIBUTE , Map<GROUPING, Map<DISTRIBUTION, Double>>> getDistributionPosition(ProteinSequence sequence) throws Exception;
}
