package org.biojava3.aaproperties;

import java.util.Map;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;

/**
 * TODO Javadoc? You already have most of it in IPeptideProperties you just need to add a reference here     
 * 
 */
public class PeptideProperties {
	public enum SingleLetterAACode { A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V}
	
	public static final double getMolecularWeight(String sequence){
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getMolecularWeight(pSequence);
	}

	public static final double getExtinctionCoefficient(String sequence, boolean assumeCysReduced) {
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getExtinctionCoefficient(pSequence, assumeCysReduced);
	}

	public static final double getInstabilityIndex(String sequence) {
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getInstabilityIndex(pSequence);
	}

	public static final double getApliphaticIndex(String sequence) {
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getApliphaticIndex(pSequence);
	}

	public static final double getAvgHydropathy(String sequence) {
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getAvgHydropathy(pSequence);
	}

	public static final double getIsoelectricPoint(String sequence) {
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getIsoelectricPoint(pSequence);
	}

	public static final double getNetCharge(String sequence) {
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getNetCharge(pSequence);
	}

	/**
	 * TODO complete
	 */
	public static final double getEnrichment(String sequence, AminoAcidCompound aminoAcid) {
		throw new UnsupportedOperationException();
	}

	public static final double getEnrichment(String sequence, SingleLetterAACode aminoAcidCode) {
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		return pp.getEnrichment(pSequence, aaSet.getCompoundForString(aminoAcidCode.toString()));
	}

	public static final Map<AminoAcidCompound, Double> getAAComposition(String sequence) {
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		return pp.getAAComposition(pSequence);
	}
	
	/**
	 * TODO complete
	 * Can you come up with a better name? 
	 */
	public static final Map<SingleLetterAACode, Double> getAAComposition2(String sequence) {
		throw new UnsupportedOperationException();
	}
	
}
