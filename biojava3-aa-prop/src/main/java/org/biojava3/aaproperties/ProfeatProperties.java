package org.biojava3.aaproperties;

import java.util.Map;

import org.biojava3.aaproperties.IProfeatProperties.ATTRIBUTE;
import org.biojava3.aaproperties.IProfeatProperties.DISTRIBUTION;
import org.biojava3.aaproperties.IProfeatProperties.GROUPING;
import org.biojava3.aaproperties.IProfeatProperties.TRANSITION;
import org.biojava3.core.sequence.ProteinSequence;

/**
 * This is an adaptor class which enable the ease of generating profeat properties.
 * At least one adaptor method is written for each available properties provided in IProfeatProperties. 
 * 
 * @author kohchuanhock
 * @version 2011.06.16
 * @see IProfeatProperties
 * @see ProfeatPropertiesImpl
 */
public class ProfeatProperties {
	/**
	 * An adaptor method which returns the composition of the specific grouping for the given attribute.
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
	public static double getComposition(ProteinSequence sequence, ATTRIBUTE attribute, GROUPING group) throws Exception{
		return new ProfeatPropertiesImpl().getComposition(sequence, attribute, group);
	}
	
	public static Map<GROUPING, Double> getComposition(ProteinSequence sequence, ATTRIBUTE attribute) throws Exception{
		return new ProfeatPropertiesImpl().getComposition(sequence, attribute);
	}
	
	public static Map<ATTRIBUTE, Map<GROUPING, Double>> getComposition(ProteinSequence sequence) throws Exception{
		return new ProfeatPropertiesImpl().getComposition(sequence);
	}
	
	/**
	 * An adaptor method which returns the number of transition between the specified groups for the given attribute with respect to the length of sequence.
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
	public static double getTransition(ProteinSequence sequence, ATTRIBUTE attribute, TRANSITION transition) throws Exception{
		return new ProfeatPropertiesImpl().getTransition(sequence, attribute, transition);
	}
	
	public static Map<TRANSITION, Double> getTransition(ProteinSequence sequence, ATTRIBUTE attribute) throws Exception{
		return new ProfeatPropertiesImpl().getTransition(sequence, attribute);
	}
	
	public static Map<ATTRIBUTE, Map<TRANSITION, Double>> getTransition(ProteinSequence sequence) throws Exception{
		return new ProfeatPropertiesImpl().getTransition(sequence);
	}
	
	/**
	 * An adaptor method which computes and return the position with respect to the sequence where the given distribution of the grouping can be found.<br/>
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
	public static double getPosition(ProteinSequence sequence, ATTRIBUTE attribute, GROUPING group, DISTRIBUTION distribution) throws Exception{
		return new ProfeatPropertiesImpl().getPosition(sequence, attribute, group, distribution);
	}
	
	public static Map<DISTRIBUTION, Double> getPosition(ProteinSequence sequence, ATTRIBUTE attribute, GROUPING group) throws Exception{
		return new ProfeatPropertiesImpl().getPosition(sequence, attribute, group);
	}
	
	public static Map<GROUPING, Map<DISTRIBUTION, Double>> getPosition(ProteinSequence sequence, ATTRIBUTE attribute) throws Exception{
		return new ProfeatPropertiesImpl().getPosition(sequence, attribute);
	}
	
	public static Map<ATTRIBUTE , Map<GROUPING, Map<DISTRIBUTION, Double>>> getPosition(ProteinSequence sequence) throws Exception{
		return new ProfeatPropertiesImpl().getPosition(sequence);
	}
}
