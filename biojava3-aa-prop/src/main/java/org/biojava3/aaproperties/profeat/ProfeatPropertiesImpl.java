package org.biojava3.aaproperties.profeat;

import java.util.HashMap;
import java.util.Map;

import org.biojava3.aaproperties.profeat.convertor.Convert2Charge;
import org.biojava3.aaproperties.profeat.convertor.Convert2Hydrophobicity;
import org.biojava3.aaproperties.profeat.convertor.Convert2NormalizedVanDerWaalsVolume;
import org.biojava3.aaproperties.profeat.convertor.Convert2Polarity;
import org.biojava3.aaproperties.profeat.convertor.Convert2Polarizability;
import org.biojava3.aaproperties.profeat.convertor.Convert2SecondaryStructure;
import org.biojava3.aaproperties.profeat.convertor.Convert2SolventAccessibility;
import org.biojava3.aaproperties.profeat.convertor.Convertor;
import org.biojava3.core.sequence.ProteinSequence;

public class ProfeatPropertiesImpl implements IProfeatProperties{

	@Override
	public double getComposition(ProteinSequence sequence, ATTRIBUTE attribute, GROUPING group) throws Exception {
		Convertor convertor = getConvertor(attribute);
		String convertedSequence = convertor.convert(sequence);
		return (getTotalCount(convertedSequence, group) + 0.0) / convertedSequence.length();
	}

	private int getTotalCount(String convertedSeq, GROUPING group) throws Exception{
		char g;
		switch(group){
		case GROUP1: g = Convertor.group1; break;
		case GROUP2: g = Convertor.group2; break;
		case GROUP3: g = Convertor.group3; break;
		default: throw new Exception("Unhandled Case: " + group);
		}
		int total = 0;
		for(char c:convertedSeq.toCharArray()){
			if(c == g) total++;
		}
		return total;
	}

	@Override
	public double getTransition(ProteinSequence sequence, ATTRIBUTE attribute, TRANSITION transition) throws Exception{
		Convertor convertor = getConvertor(attribute);
		String convertedSequence = convertor.convert(sequence);
		char t1;
		char t2;
		switch(transition){
		case BETWEEN_11: t1 = '1'; t2 = '1'; break;
		case BETWEEN_22: t1 = '2'; t2 = '2'; break;
		case BETWEEN_33: t1 = '3'; t2 = '3'; break;
		case BETWEEN_12: t1 = '1'; t2 = '2'; break;
		case BETWEEN_13: t1 = '1'; t2 = '3'; break;
		case BETWEEN_23: t1 = '2'; t2 = '3'; break;
		default: throw new Exception("Unhandled Case: " + transition);
		}
		int total = 0;
		for(int i = 0; i < convertedSequence.length() - 1; i++){
			char s1 = convertedSequence.charAt(i);
			char s2 = convertedSequence.charAt(i + 1);
			if((t1 == s1 && t2 == s2) || (t2 == s1 && t1 == s2)) total++;
		}
		return total / (convertedSequence.length() - 1.0);
	}

	@Override
	public double getDistributionPosition(ProteinSequence sequence, ATTRIBUTE attribute, GROUPING group, DISTRIBUTION distribution) throws Exception{
		Convertor convertor = getConvertor(attribute);
		String convertedSequence = convertor.convert(sequence);

		int totalCount = getTotalCount(convertedSequence, group);
		int countIndex;
		switch(distribution){
		case FIRST: countIndex = 1; break;
		case FIRST25: countIndex = totalCount * 25 / 100; break;
		case FIRST50: countIndex = totalCount * 50 / 100; break;
		case FIRST75: countIndex = totalCount * 75 / 100; break;
		case ALL: countIndex = totalCount; break;
		default: throw new Exception("Unhandled Case: " + distribution);
		}

		char g;
		switch(group){
		case GROUP1: g = Convertor.group1; break;
		case GROUP2: g = Convertor.group2; break;
		case GROUP3: g = Convertor.group3; break;
		default: throw new Exception("Unhandled Case: " + group);
		}

		int currentCount = 0;
		int dist = 0;
		for(int x = 0; x < convertedSequence.length(); x++){
			if(g == convertedSequence.charAt(x)){
				currentCount++;				
				if(currentCount == countIndex){
					dist = x+1;
					break;					
				}				
			}	
		}		
		return (dist + 0.0) / convertedSequence.length();
	}

	private Convertor getConvertor(ATTRIBUTE attribute) throws Exception{
		switch(attribute){
		case HYDROPHOBICITY: return new Convert2Hydrophobicity();
		case VOLUME: return new Convert2NormalizedVanDerWaalsVolume();
		case POLARITY: return new Convert2Polarity();
		case POLARIZABILITY: return new Convert2Polarizability();
		case CHARGE: return new Convert2Charge();
		case SECONDARYSTRUCTURE: return new Convert2SecondaryStructure();
		case SOLVENTACCESSIBILITY: return new Convert2SolventAccessibility();
		default: throw new Exception("Unknown attribute: " + attribute);
		}
	}

	@Override
	public Map<GROUPING, Double> getComposition(ProteinSequence sequence, ATTRIBUTE attribute) throws Exception {
		Map<GROUPING, Double> grouping2Composition = new HashMap<GROUPING, Double>();
		for(GROUPING group:GROUPING.values()) grouping2Composition.put(group, getComposition(sequence, attribute, group));
		return grouping2Composition;
	}

	@Override
	public Map<ATTRIBUTE, Map<GROUPING, Double>> getComposition(ProteinSequence sequence) throws Exception {
		Map<ATTRIBUTE, Map<GROUPING, Double>> attribute2Grouping2Composition = new HashMap<ATTRIBUTE, Map<GROUPING, Double>>();
		for(ATTRIBUTE attribute:ATTRIBUTE.values()) attribute2Grouping2Composition.put(attribute, getComposition(sequence, attribute));
		return attribute2Grouping2Composition;
	}

	@Override
	public Map<TRANSITION, Double> getTransition(ProteinSequence sequence, ATTRIBUTE attribute) throws Exception {
		Map<TRANSITION, Double> transition2Double = new HashMap<TRANSITION, Double>();
		for(TRANSITION transition:TRANSITION.values()) transition2Double.put(transition, getTransition(sequence, attribute, transition));
		return transition2Double;
	}

	@Override
	public Map<ATTRIBUTE, Map<TRANSITION, Double>> getTransition(ProteinSequence sequence) throws Exception {
		Map<ATTRIBUTE, Map<TRANSITION, Double>> attribute2Transition2Double = new HashMap<ATTRIBUTE, Map<TRANSITION, Double>>();
		for(ATTRIBUTE attribute:ATTRIBUTE.values()) attribute2Transition2Double.put(attribute, getTransition(sequence, attribute));
		return attribute2Transition2Double;
	}

	@Override
	public Map<DISTRIBUTION, Double> getDistributionPosition(ProteinSequence sequence, ATTRIBUTE attribute, GROUPING group) throws Exception {
		Map<DISTRIBUTION, Double> distribution2Double = new HashMap<DISTRIBUTION, Double>();
		for(DISTRIBUTION distribution:DISTRIBUTION.values()) 
			distribution2Double.put(distribution, getDistributionPosition(sequence, attribute, group, distribution));
		return distribution2Double;
	}

	@Override
	public Map<GROUPING, Map<DISTRIBUTION, Double>> getDistributionPosition(ProteinSequence sequence, ATTRIBUTE attribute) throws Exception {
		Map<GROUPING, Map<DISTRIBUTION, Double>> grouping2Distribution2Double = new HashMap<GROUPING, Map<DISTRIBUTION, Double>>();
		for(GROUPING group:GROUPING.values()) 
			grouping2Distribution2Double.put(group, getDistributionPosition(sequence, attribute, group));
		return grouping2Distribution2Double;
	}

	@Override
	public Map<ATTRIBUTE, Map<GROUPING, Map<DISTRIBUTION, Double>>> getDistributionPosition(ProteinSequence sequence) throws Exception {
		Map<ATTRIBUTE, Map<GROUPING, Map<DISTRIBUTION, Double>>> attribute2Grouping2Distribution2Double =
			new HashMap<ATTRIBUTE, Map<GROUPING, Map<DISTRIBUTION, Double>>>();
		for(ATTRIBUTE attribute:ATTRIBUTE.values())
			attribute2Grouping2Distribution2Double.put(attribute, getDistributionPosition(sequence, attribute));
		return attribute2Grouping2Distribution2Double;
	}

}
