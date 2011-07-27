package org.biojava3.aaproperties;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;

import org.biojava3.aaproperties.xml.AminoAcidCompositionTable;
import org.biojava3.aaproperties.xml.ElementTable;
import org.biojava3.aaproperties.xml.MyValidationEventHandler;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;

public class PeptidePropertiesImpl implements IPeptideProperties{
	private double getWaterMoleculeWeight(){
		final double hydrogenMW = 1.0079;
		final double hydroxideMW = 17.0073;
		//H	1.0079	OH	17.0073
		return hydrogenMW + hydroxideMW;
	}
	
	private char[] getSequence(String sequence, boolean ignoreCase){
		if(ignoreCase){
			return sequence.toUpperCase().toCharArray();
		}else{
			return sequence.toCharArray();
		}
	}
	
	@Override
	public double getMolecularWeight(ProteinSequence sequence, boolean ignoreCase) {
		double value = 0.0;
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		char[] seq = getSequence(sequence.toString(), ignoreCase);
		for(char aa:seq){
			AminoAcidCompound c = aaSet.getCompoundForString(aa + "");
			if(Constraints.aa2MolecularWeight.containsKey(c)){
				value += Constraints.aa2MolecularWeight.get(c);
			}
		}
		if(value == 0)
			return value;
		else
			return value + getWaterMoleculeWeight();
	}
	
	@Override
	public double getMolecularWeight(ProteinSequence sequence, File aminoAcidCompositionFile, boolean ignoreCase) throws JAXBException, FileNotFoundException {
		File elementMassFile = new File("./src/main/resources/ElementMass.xml");
		if(elementMassFile.exists() == false){
			throw new FileNotFoundException("Cannot locate ElementMass.xml. " +
					"Please use getMolecularWeight(ProteinSequence, File, File) to specify ElementMass.xml location.");
		}
		return getMolecularWeightBasedOnXML(sequence, obtainAminoAcidCompositionTable(elementMassFile, aminoAcidCompositionFile, ignoreCase), ignoreCase);
	}
	
	@Override
	public double getMolecularWeight(ProteinSequence sequence, File elementMassFile, File aminoAcidCompositionFile, boolean ignoreCase) 
			throws JAXBException, FileNotFoundException{
		return getMolecularWeightBasedOnXML(sequence, obtainAminoAcidCompositionTable(elementMassFile, aminoAcidCompositionFile, ignoreCase), ignoreCase);
	}
	
	@Override
	public double getMolecularWeightBasedOnXML(ProteinSequence sequence, AminoAcidCompositionTable aminoAcidCompositionTable, boolean ignoreCase){
		double value = 0.0;
		char[] seq = getSequence(sequence.toString(), ignoreCase);
		for(char aa:seq){
			Double weight = aminoAcidCompositionTable.getMolecularWeight(aa);
			if(weight != null){
				value += weight;
			}
		}
		if(value == 0.0)
			return value;
		else
			return value + getWaterMoleculeWeight();
	}

	@Override
	public AminoAcidCompositionTable obtainAminoAcidCompositionTable(File aminoAcidCompositionFile, boolean ignoreCase) 
		throws JAXBException, FileNotFoundException{
		File elementMassFile = new File("./src/main/resources/ElementMass.xml");
		if(elementMassFile.exists() == false){
			throw new FileNotFoundException("Cannot locate ElementMass.xml. " +
					"Please use getMolecularWeight(ProteinSequence, File, File) to specify ElementMass.xml location.");
		}
		return obtainAminoAcidCompositionTable(elementMassFile, aminoAcidCompositionFile, ignoreCase);
	}
	
	@Override
	public AminoAcidCompositionTable obtainAminoAcidCompositionTable(File elementMassFile, File aminoAcidCompositionFile, boolean ignoreCase) 
		throws JAXBException, FileNotFoundException{
		//Parse elementMassFile
		ElementTable iTable = new ElementTable();
		// Get a JAXB Context for the object we created above
		JAXBContext jc = JAXBContext.newInstance(iTable.getClass());
		Unmarshaller u = jc.createUnmarshaller();
		u.setEventHandler(new MyValidationEventHandler()); 
		iTable = (ElementTable)u.unmarshal(new FileInputStream(elementMassFile));
		iTable.populateMaps();
		
		//Parse aminoAcidCompositionFile
		AminoAcidCompositionTable aTable = new AminoAcidCompositionTable();
		// Get a JAXB Context for the object we created above
		JAXBContext jc2 = JAXBContext.newInstance(aTable.getClass());
		Unmarshaller u2 = jc2.createUnmarshaller();
		u2.setEventHandler(new MyValidationEventHandler()); 
		aTable = (AminoAcidCompositionTable)u2.unmarshal(new FileInputStream(aminoAcidCompositionFile));
		aTable.computeMolecularWeight(iTable, ignoreCase);
		return aTable;
	}

	@Override
	public double getExtinctionCoefficient(ProteinSequence sequence, boolean assumeCysReduced, boolean ignoreCase) {
		//Tyr => Y
		//Trp => W
		//Cys => C
		//E(Prot) = Numb(Tyr)*Ext(Tyr) + Numb(Trp)*Ext(Trp) + Numb(Cystine)*Ext(Cystine)
		//where (for proteins in water measured at 280 nm): Ext(Tyr) = 1490, Ext(Trp) = 5500, Ext(Cystine) = 125;
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		Map<AminoAcidCompound, Integer> extinctAA2Count = this.getExtinctAACount(sequence, ignoreCase);
		
		double eProt;
		if(assumeCysReduced == false){
			eProt = extinctAA2Count.get(aaSet.getCompoundForString("Y")) * 
				Constraints.aa2ExtinctionCoefficient.get(aaSet.getCompoundForString("Y")) + 
				extinctAA2Count.get(aaSet.getCompoundForString("W")) * 
				Constraints.aa2ExtinctionCoefficient.get(aaSet.getCompoundForString("W")) +
				extinctAA2Count.get(aaSet.getCompoundForString("C")) * 
				Constraints.aa2ExtinctionCoefficient.get(aaSet.getCompoundForString("C"));
		}else
			eProt = extinctAA2Count.get(aaSet.getCompoundForString("Y")) * 
				Constraints.aa2ExtinctionCoefficient.get(aaSet.getCompoundForString("Y")) + 
				extinctAA2Count.get(aaSet.getCompoundForString("W")) *
				Constraints.aa2ExtinctionCoefficient.get(aaSet.getCompoundForString("W"));
		
		return eProt;
	}
	
	@Override
	public double getAbsorbance(ProteinSequence sequence, boolean assumeCysReduced, boolean ignoreCase){
		//Absorb(Prot) = E(Prot) / Molecular_weight
		double mw = this.getMolecularWeight(sequence, ignoreCase);
		double eProt = this.getExtinctionCoefficient(sequence, assumeCysReduced, ignoreCase);
		return eProt / mw;
	}
	
	private Map<AminoAcidCompound, Integer> getExtinctAACount(ProteinSequence sequence, boolean ignoreCase){
		//Cys => C, Tyr => Y, Trp => W
		int numW = 0;
		int smallW = 0;
		double numC = 0;
		double smallC = 0; 
		int numY = 0;
		int smallY = 0;
		for(char aa:sequence.getSequenceAsString().toCharArray()){
			switch(aa){
			case 'W': numW++; break;
			case 'w': smallW++; break;
			case 'C': numC += 0.5; break;
			case 'c': smallC += 0.5; break;
			case 'Y': numY++; break;
			case 'y': smallY++; break;
			}
		}
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		Map<AminoAcidCompound, Integer> extinctAA2Count = new HashMap<AminoAcidCompound, Integer>();
		if(ignoreCase){
			extinctAA2Count.put(aaSet.getCompoundForString("W"), numW + smallW);
			extinctAA2Count.put(aaSet.getCompoundForString("C"), (int) (numC + smallC));
			extinctAA2Count.put(aaSet.getCompoundForString("Y"), numY + smallY);
		}else{
			extinctAA2Count.put(aaSet.getCompoundForString("W"), numW);
			extinctAA2Count.put(aaSet.getCompoundForString("C"), (int) numC);
			extinctAA2Count.put(aaSet.getCompoundForString("Y"), numY);
		}
		return extinctAA2Count;
	}

	@Override
	public double getInstabilityIndex(ProteinSequence sequence, boolean ignoreCase) {
		double sum = 0.0;
		String s = sequence.getSequenceAsString();
		if(ignoreCase) s = s.toUpperCase();
		for(int i = 0; i < sequence.getLength() - 1; i++){
			String dipeptide = s.substring(i, i+2);
			if(Constraints.diAA2Instability.containsKey(dipeptide)){
				sum += Constraints.diAA2Instability.get(dipeptide);
			}
		}
		return sum * 10 / (s.length() - Utils.getNumberOfInvalidChar(s, null, ignoreCase));
	}

	@Override
	public double getApliphaticIndex(ProteinSequence sequence, boolean ignoreCase) {
//		Aliphatic index = X(Ala) + a * X(Val) + b * ( X(Ile) + X(Leu) )  
//		where X(Ala), X(Val), X(Ile), and X(Leu) are mole percent (100 X mole fraction) 
//		of alanine, valine, isoleucine, and leucine. 
//		The coefficients a and b are the relative volume of valine side chain (a = 2.9) 
//		and of Leu/Ile side chains (b = 3.9) to the side chain of alanine. 
//		Ala => A, Val => V, Ile => I, Leu => L
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		Map<AminoAcidCompound, Double> aa2Composition = getAAComposition(sequence, ignoreCase);
		final double a = 2.9;
		final double b = 3.9;
		double xAla = aa2Composition.get(aaSet.getCompoundForString("A"));
		double xVal = aa2Composition.get(aaSet.getCompoundForString("V"));
		double xIle = aa2Composition.get(aaSet.getCompoundForString("I"));
		double xLeu = aa2Composition.get(aaSet.getCompoundForString("L"));
		return (xAla + (a * xVal) + (b * (xIle + xLeu))) * 100;
	}

	@Override
	public double getAvgHydropathy(ProteinSequence sequence, boolean ignoreCase) {
		int validLength = 0;
		double total = 0.0;
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		char[] seq = this.getSequence(sequence.toString(), ignoreCase);
		for(char aa:seq){
			AminoAcidCompound c = aaSet.getCompoundForString(aa + "");
			if(Constraints.aa2Hydrophathicity.containsKey(c)){
				total += Constraints.aa2Hydrophathicity.get(c);
				validLength++;
			}
		}
		return total / validLength;
	}

	@Override
	public double getIsoelectricPoint(ProteinSequence sequence, boolean ignoreCase) {
		double currentPH = 7.0;
		double changeSize = 7.0;
		Map<AminoAcidCompound, Integer> chargedAA2Count = this.getChargedAACount(sequence, ignoreCase);
		double margin;
		final double difference = 0.0000001;
		while(true){
			margin = this.getNetCharge(chargedAA2Count, currentPH);
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
	public double getNetCharge(ProteinSequence sequence, boolean ignoreCase) {
		Map<AminoAcidCompound, Integer> chargedAA2Count = this.getChargedAACount(sequence, ignoreCase);
		return getNetCharge(chargedAA2Count, 7.0);
	}
	
	private double getNetCharge(Map<AminoAcidCompound, Integer> chargedAA2Count, double ph){
		//Lys => K, Arg => R, His => H
		//Asp => D, Glu => E, Cys => C, Tyr => Y
		//(NH2-)	9.69	(-COOH)	2.34
		final double pkaOfNH2 = 9.69;
		final double pkaOfCOOH = 2.34;
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		double nTerminalCharge = this.getPosCharge(pkaOfNH2, ph);
		double kCharge = chargedAA2Count.get(aaSet.getCompoundForString("K")) * this.getPosCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("K")), ph);
		double rCharge = chargedAA2Count.get(aaSet.getCompoundForString("R")) * this.getPosCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("R")), ph);
		double hCharge = chargedAA2Count.get(aaSet.getCompoundForString("H")) * this.getPosCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("H")), ph);
		double dCharge = chargedAA2Count.get(aaSet.getCompoundForString("D")) * this.getNegCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("D")), ph);
		double eCharge = chargedAA2Count.get(aaSet.getCompoundForString("E")) * this.getNegCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("E")), ph);
		double cCharge = chargedAA2Count.get(aaSet.getCompoundForString("C")) * this.getNegCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("C")), ph);
		double yCharge = chargedAA2Count.get(aaSet.getCompoundForString("Y")) * this.getNegCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("Y")), ph);
		double cTerminalCharge = this.getNegCharge(pkaOfCOOH, ph);
		if((kCharge + rCharge + hCharge) == 0.0 && (dCharge + eCharge + cCharge + yCharge) == 0.0){
			return 0.0;
		}
		return (nTerminalCharge + kCharge + rCharge + hCharge) - (dCharge + eCharge + cCharge + yCharge + cTerminalCharge);
	}
	
	private double getPosCharge(double pka, double ph){
		return Math.pow(10, pka) / (Math.pow(10, pka) + Math.pow(10, ph));
	}
	
	private double getNegCharge(double pka, double ph){
		return Math.pow(10, ph) / (Math.pow(10, pka) + Math.pow(10, ph));
	}
	
	private Map<AminoAcidCompound, Integer> getChargedAACount(ProteinSequence sequence, boolean ignoreCase){
		//Lys => K, Arg => R, His => H
		//Asp => D, Glu => E, Cys => C, Tyr => Y
		int numK = 0;
		int numR = 0;
		int numH = 0;
		int numD = 0;
		int numE = 0;
		int numC = 0;
		int numY = 0;
		char[] seq = this.getSequence(sequence.getSequenceAsString(), ignoreCase);
		for(char aa:seq){
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
	public double getEnrichment(ProteinSequence sequence, AminoAcidCompound aminoAcidCode, boolean ignoreCase) {
		double counter = 0.0;
		char[] seq = this.getSequence(sequence.getSequenceAsString(), ignoreCase);
		for(char aa:seq){
			if(aminoAcidCode.getShortName().equals(aa + "")){
				counter++;
			}
		}
		return counter/sequence.getLength();
	}

	@Override
	public Map<AminoAcidCompound, Double> getAAComposition(ProteinSequence sequence, boolean ignoreCase) {
		int validLength = 0;
		Map<AminoAcidCompound, Double> aa2Composition = new HashMap<AminoAcidCompound, Double>();
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		for(AminoAcidCompound aa:aaSet.getAllCompounds()){
			aa2Composition.put(aa, 0.0);
		}
		char[] seq = this.getSequence(sequence.toString(), ignoreCase);
		for(char aa:seq){
			if(PeptideProperties.standardAASet.contains(aa)){
				AminoAcidCompound compound = aaSet.getCompoundForString(aa + "");
				aa2Composition.put(compound, aa2Composition.get(compound) + 1.0);
				validLength++;
			}
		}
		if(validLength > 0){
			for(AminoAcidCompound aa:aaSet.getAllCompounds()){
				aa2Composition.put(aa, aa2Composition.get(aa) / validLength);
			}
		}else{
			for(AminoAcidCompound aa:aaSet.getAllCompounds()){
				aa2Composition.put(aa, 0.0);
			}
		}
		return aa2Composition;
	}
	
	/*
	 * Quick Test Method
	 */
	public static void main(String[] args){
		ProteinSequence sequence = new ProteinSequence("AAACCAAAWWTT");
		IPeptideProperties pp = new PeptidePropertiesImpl();
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		System.out.println(aaSet.getCompoundForString("1"));
		System.out.println(sequence);
		boolean ignoreCase = true;
		System.out.println("A Composition: " + pp.getEnrichment(sequence, aaSet.getCompoundForString("A"), ignoreCase));//CHECKED
		
		System.out.println("C Composition: " + pp.getEnrichment(sequence, aaSet.getCompoundForString("C"), ignoreCase));//CHECKED
		System.out.println("Molecular Weight: " + pp.getMolecularWeight(sequence, ignoreCase));//CHECKED
		System.out.println("All Composition: " + pp.getAAComposition(sequence, ignoreCase));//CHECKED
		System.out.println("AvgHydro: " + pp.getAvgHydropathy(sequence, ignoreCase));//CHECKED
		System.out.println("NetCharge: " + pp.getNetCharge(sequence, ignoreCase));//CHECKED
		System.out.println("PI: " + pp.getIsoelectricPoint(sequence, ignoreCase));//CHECKED
		System.out.println("Apliphatic: " + pp.getApliphaticIndex(sequence, ignoreCase));//CHECKED
		System.out.println("Instability: " + pp.getInstabilityIndex(sequence, ignoreCase));//CHECKED
		System.out.println("Extinct: " + pp.getExtinctionCoefficient(sequence, true, ignoreCase));//CHECKED - 
		System.out.println("Extinct: " + pp.getExtinctionCoefficient(sequence, false, ignoreCase));
	}
}
