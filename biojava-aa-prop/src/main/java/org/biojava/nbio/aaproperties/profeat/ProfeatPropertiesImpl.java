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
package org.biojava.nbio.aaproperties.profeat;

import org.biojava.nbio.aaproperties.profeat.convertor.*;
import org.biojava.nbio.core.sequence.ProteinSequence;

import java.util.HashMap;
import java.util.Map;

public class ProfeatPropertiesImpl implements IProfeatProperties{

	@Override
	public double getComposition(ProteinSequence sequence, ATTRIBUTE attribute, GROUPING group) throws Exception {
		Convertor convertor = getConvertor(attribute);
		String convertedSequence = convertor.convert(sequence);
		return (getTotalCount(convertedSequence, group) + 0.0) / convertedSequence.length();
	}

	private int getTotalCount(String convertedSeq, GROUPING group) throws Exception{
		char g = getGroup(group);
		int total = 0;
		total = (int)convertedSeq.chars().filter(c ->(char) c == g) .count();
		return total;
	}

	private char getGroup(GROUPING group) {
		char g;
		switch(group){
			case GROUP1: g = Convertor.group1; break;
			case GROUP2: g = Convertor.group2; break;
			case GROUP3: g = Convertor.group3; break;
			default: throw new UnsupportedOperationException("Unhandled Case: " + group);
		}
		return g;
	}

	@Override
	public double getTransition(ProteinSequence sequence, ATTRIBUTE attribute, TRANSITION transition) throws Exception{
		Convertor convertor = getConvertor(attribute);
		String convertedSequence = convertor.convert(sequence);
		char transition1;
		char transition2;
		switch(transition){
		case BETWEEN_11: transition1 = '1'; transition2 = '1'; break;
		case BETWEEN_22: transition1 = '2'; transition2 = '2'; break;
		case BETWEEN_33: transition1 = '3'; transition2 = '3'; break;
		case BETWEEN_12: transition1 = '1'; transition2 = '2'; break;
		case BETWEEN_13: transition1 = '1'; transition2 = '3'; break;
		case BETWEEN_23: transition1 = '2'; transition2 = '3'; break;
		default: throw new Exception("Unhandled Case: " + transition);
		}
		int total = 0;
		for(int i = 0; i < convertedSequence.length() - 1; i++){
			char firstSequence = convertedSequence.charAt(i);
			char secondSequence = convertedSequence.charAt(i + 1);
			if((transition1 == firstSequence && transition2 == secondSequence) || (transition2 == firstSequence && transition1 == secondSequence)) total++;
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

		char g = getGroup(group);

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
