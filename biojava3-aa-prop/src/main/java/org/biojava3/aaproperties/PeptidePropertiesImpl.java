package org.biojava3.aaproperties;

import java.util.Map;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.AminoAcidImpl;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;

public class PeptidePropertiesImpl implements IPeptideProperties{

	@Override
	public double getMolecularWeight(ProteinSequence sequence) {
		double value = 0.0;
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		for(char aa:sequence.toString().toCharArray()){
			value += Constraints.aa2MolecularWeight.get(aaSet.getCompoundForString(aa + ""));
		}
		return value;
	}

	@Override
	public double getExtinctionCoefficient(ProteinSequence sequence,
			boolean assumeCysReduced) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getInstabilityIndex(ProteinSequence sequence) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getApliphaticIndex(ProteinSequence sequence) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getAvgHydropathy(ProteinSequence sequence) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getIsoelectricPoint(ProteinSequence sequence) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getNetCharge(ProteinSequence sequence) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getEnrichment(ProteinSequence sequence, AminoAcid aminoAcidCode) {
		double counter = 0.0;
		for(char aa:sequence.getSequenceAsString().toCharArray()){
			if(aminoAcidCode.getAminoType() == aa){
				counter++;
			}
		}
		return counter/sequence.getLength();
	}

	@Override
	public Map<AminoAcidCompound, Double> getAAComposition(ProteinSequence sequence) {
		// TODO Auto-generated method stub
		return null;
	}

	/*
	 * Test Method
	 */
	public static void main(String[] args){
		ProteinSequence sequence = new ProteinSequence("AAACCCAAA");
		IPeptideProperties pp = new PeptidePropertiesImpl();
		AminoAcid aa = new AminoAcidImpl();
		aa.setAminoType('A');
		System.out.println(pp.getEnrichment(sequence, aa));
		System.out.println(pp.getMolecularWeight(sequence));
	}
}
