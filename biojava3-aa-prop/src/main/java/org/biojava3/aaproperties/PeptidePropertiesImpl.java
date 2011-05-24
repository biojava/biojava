package org.biojava3.aaproperties;

import java.util.HashMap;
import java.util.Map;

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
	public double getExtinctionCoefficient(ProteinSequence sequence, boolean assumeCysReduced) {
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
//		Aliphatic index = X(Ala) + a * X(Val) + b * ( X(Ile) + X(Leu) )  
//		where X(Ala), X(Val), X(Ile), and X(Leu) are mole percent (100 X mole fraction) 
//		of alanine, valine, isoleucine, and leucine. 
//		The coefficients a and b are the relative volume of valine side chain (a = 2.9) 
//		and of Leu/Ile side chains (b = 3.9) to the side chain of alanine. 
//		Ala => A, Val => V, Ile => I, Leu => L
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		Map<AminoAcidCompound, Double> aa2Composition = getAAComposition(sequence);
		double a = 2.9;
		double b = 3.9;
		double xAla = aa2Composition.get(aaSet.getCompoundForString("A"));
		double xVal = aa2Composition.get(aaSet.getCompoundForString("V"));
		double xIle = aa2Composition.get(aaSet.getCompoundForString("I"));
		double xLeu = aa2Composition.get(aaSet.getCompoundForString("L"));
		return xAla + (a * xVal) + (b * (xIle + xLeu));
	}

	@Override
	public double getAvgHydropathy(ProteinSequence sequence) {
		double value = 0.0;
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		for(char aa:sequence.toString().toCharArray()){
			value += Constraints.aa2Hydrophathicity.get(aaSet.getCompoundForString(aa + ""));
		}
		return value / sequence.getLength();
	}

	@Override
	public double getIsoelectricPoint(ProteinSequence sequence) {
		double currentPH = 7.0;
		double changeSize = 7.0;
		Map<AminoAcidCompound, Integer> chargedAA2Count = this.getChargedAACount(sequence);
		double margin = 1.0;
		double difference = 0.00001;
		while(true){
			margin = this.getNetCharge(chargedAA2Count, currentPH) - 0.0;
			//Within allowed difference
			if(margin <= difference && margin >= -difference) break;
			changeSize /= 2.0;
			if(margin > 0){
				currentPH += changeSize;
			}else{
				currentPH -= changeSize;
			}
		}
		return currentPH;
	}

	@Override
	public double getNetCharge(ProteinSequence sequence) {
		Map<AminoAcidCompound, Integer> chargedAA2Count = this.getChargedAACount(sequence);
		return getNetCharge(chargedAA2Count, 7.0);
	}
	
	private double getNetCharge(Map<AminoAcidCompound, Integer> chargedAA2Count, double ph){
		//Lys => K, Arg => R, His => H
		//Asp => D, Glu => E, Cys => C, Tyr => Y
		//(NH2-)	9.69	(-COOH)	2.34
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		double nTerminalCharge = this.getPosCharge(9.69, ph);
		double kCharge = chargedAA2Count.get(aaSet.getCompoundForString("K")) * this.getPosCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("K")), ph);
		double rCharge = chargedAA2Count.get(aaSet.getCompoundForString("R")) * this.getPosCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("R")), ph);
		double hCharge = chargedAA2Count.get(aaSet.getCompoundForString("H")) * this.getPosCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("H")), ph);
		double dCharge = chargedAA2Count.get(aaSet.getCompoundForString("D")) * this.getNegCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("D")), ph);
		double eCharge = chargedAA2Count.get(aaSet.getCompoundForString("E")) * this.getNegCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("E")), ph);
		double cCharge = chargedAA2Count.get(aaSet.getCompoundForString("C")) * this.getNegCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("C")), ph);
		double yCharge = chargedAA2Count.get(aaSet.getCompoundForString("Y")) * this.getNegCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("Y")), ph);
		double cTerminalCharge = this.getNegCharge(2.34, ph);
		return (nTerminalCharge + kCharge + rCharge + hCharge) - (dCharge + eCharge + cCharge + yCharge + cTerminalCharge);
	}
	
	private double getPosCharge(double pka, double ph){
		return Math.pow(10, pka) / (Math.pow(10, pka) + Math.pow(10, ph));
	}
	
	private double getNegCharge(double pka, double ph){
		return Math.pow(10, ph) / (Math.pow(10, pka) + Math.pow(10, ph));
	}
	
	private Map<AminoAcidCompound, Integer> getChargedAACount(ProteinSequence sequence){
		//Lys => K, Arg => R, His => H
		//Asp => D, Glu => E, Cys => C, Tyr => Y
		int numK = 0;
		int numR = 0;
		int numH = 0;
		int numD = 0;
		int numE = 0;
		int numC = 0;
		int numY = 0;
		for(char aa:sequence.getSequenceAsString().toCharArray()){
			switch(aa){
			case 'K': numK++; break;
			case 'R': numR++; break;
			case 'H': numH++; break;
			case 'D': numD++; break;
			case 'E': numE++; break;
			case 'C': numC++; break;
			case 'Y': numY++; break;
			}
		}
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		Map<AminoAcidCompound, Integer> chargedAA2Count = new HashMap<AminoAcidCompound, Integer>();
		chargedAA2Count.put(aaSet.getCompoundForString("K"), numK);
		chargedAA2Count.put(aaSet.getCompoundForString("R"), numR);
		chargedAA2Count.put(aaSet.getCompoundForString("H"), numH);
		chargedAA2Count.put(aaSet.getCompoundForString("D"), numD);
		chargedAA2Count.put(aaSet.getCompoundForString("E"), numE);
		chargedAA2Count.put(aaSet.getCompoundForString("C"), numC);
		chargedAA2Count.put(aaSet.getCompoundForString("Y"), numY);
		return chargedAA2Count;
	}
	
	@Override
	public double getEnrichment(ProteinSequence sequence, AminoAcidCompound aminoAcidCode) {
		double counter = 0.0;
		for(char aa:sequence.getSequenceAsString().toCharArray()){
			if(aminoAcidCode.getShortName().equals(aa + "")){
				counter++;
			}
		}
		return counter/sequence.getLength();
	}

	@Override
	public Map<AminoAcidCompound, Double> getAAComposition(ProteinSequence sequence) {
		Map<AminoAcidCompound, Double> aa2Composition = new HashMap<AminoAcidCompound, Double>();
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		for(AminoAcidCompound aa:aaSet.getAllCompounds()){
			aa2Composition.put(aa, 0.0);
		}
		for(char aa:sequence.toString().toCharArray()){
			AminoAcidCompound compound = aaSet.getCompoundForString(aa + "");
			aa2Composition.put(compound, aa2Composition.get(compound) + 1.0);
		}
		for(AminoAcidCompound aa:aaSet.getAllCompounds()){
			aa2Composition.put(aa, aa2Composition.get(aa) / sequence.getLength());
		}
		return aa2Composition;
	}

	/*
	 * Test Method
	 */
	public static void main(String[] args){
		ProteinSequence sequence = new ProteinSequence("AAACCCAAA");
		IPeptideProperties pp = new PeptidePropertiesImpl();
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		System.out.println("A Composition: " + pp.getEnrichment(sequence, aaSet.getCompoundForString("A")));//CHECKED
		System.out.println("C Composition: " + pp.getEnrichment(sequence, aaSet.getCompoundForString("C")));//CHECKED
		System.out.println("Molecular Weight: " + pp.getMolecularWeight(sequence));//CHECKED
		System.out.println("All Composition: " + pp.getAAComposition(sequence));//CHECKED
		System.out.println("AvgHydro: " + pp.getAvgHydropathy(sequence));//CHECKED
		System.out.println("NetCharge: " + pp.getNetCharge(sequence));//CHECKED
		System.out.println("PI: " + pp.getIsoelectricPoint(sequence));
	}
}
