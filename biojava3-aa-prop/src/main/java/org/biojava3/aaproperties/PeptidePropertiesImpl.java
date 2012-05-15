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

/**
 * This class contains the actual implementation of IPeptideProperties and is wrapped around by PeptideProperties for ease of use. 
 * 
 * @author kohchuanhock
 * @version 2011.08.22
 * @since 3.0.2
 * @see IPeptideProperties
 * @see PeptideProperties
 */
public class PeptidePropertiesImpl implements IPeptideProperties{
	
	/**
	 * @return the molecular weight of water
	 */
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
	public double getMolecularWeight(ProteinSequence sequence) {
		double value = 0.0;
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		char[] seq = getSequence(sequence.toString(), true);//ignore case
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
	public double getMolecularWeight(ProteinSequence sequence, File aminoAcidCompositionFile) throws JAXBException, FileNotFoundException {
		File elementMassFile = new File("./src/main/resources/ElementMass.xml");
		if(!elementMassFile.exists()){
			throw new FileNotFoundException("Cannot locate ElementMass.xml. " +
					"Please use getMolecularWeight(ProteinSequence, File, File) to specify ElementMass.xml location.");
		}
		return getMolecularWeightBasedOnXML(sequence, obtainAminoAcidCompositionTable(elementMassFile, aminoAcidCompositionFile));
	}
	
	@Override
	public double getMolecularWeight(ProteinSequence sequence, File elementMassFile, File aminoAcidCompositionFile) 
			throws JAXBException, FileNotFoundException{
		return getMolecularWeightBasedOnXML(sequence, obtainAminoAcidCompositionTable(elementMassFile, aminoAcidCompositionFile));
	}
	
	@Override
	public double getMolecularWeightBasedOnXML(ProteinSequence sequence, AminoAcidCompositionTable aminoAcidCompositionTable){
		double value = 0.0;
		char[] seq = sequence.toString().toCharArray();
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
	public AminoAcidCompositionTable obtainAminoAcidCompositionTable(File aminoAcidCompositionFile) 
		throws JAXBException, FileNotFoundException{
		File elementMassFile = new File("./src/main/resources/ElementMass.xml");
		if(!elementMassFile.exists()){
			throw new FileNotFoundException("Cannot locate ElementMass.xml. " +
					"Please use getMolecularWeight(ProteinSequence, File, File) to specify ElementMass.xml location.");
		}
		return obtainAminoAcidCompositionTable(elementMassFile, aminoAcidCompositionFile);
	}
	
	@Override
	public AminoAcidCompositionTable obtainAminoAcidCompositionTable(File elementMassFile, File aminoAcidCompositionFile) 
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
		aTable.computeMolecularWeight(iTable);
		return aTable;
	}

	@Override
	public double getExtinctionCoefficient(ProteinSequence sequence, boolean assumeCysReduced) {
		//Tyr => Y
		//Trp => W
		//Cys => C
		//E(Prot) = Numb(Tyr)*Ext(Tyr) + Numb(Trp)*Ext(Trp) + Numb(Cystine)*Ext(Cystine)
		//where (for proteins in water measured at 280 nm): Ext(Tyr) = 1490, Ext(Trp) = 5500, Ext(Cystine) = 125;
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		Map<AminoAcidCompound, Integer> extinctAA2Count = this.getExtinctAACount(sequence);
		
		double eProt;
		if(!assumeCysReduced){
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
	public double getAbsorbance(ProteinSequence sequence, boolean assumeCysReduced){
		//Absorb(Prot) = E(Prot) / Molecular_weight
		double mw = this.getMolecularWeight(sequence);
		double eProt = this.getExtinctionCoefficient(sequence, assumeCysReduced);
		return eProt / mw;
	}
	
	private Map<AminoAcidCompound, Integer> getExtinctAACount(ProteinSequence sequence){
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
		//Ignore Case is always true
		extinctAA2Count.put(aaSet.getCompoundForString("W"), numW + smallW);
		extinctAA2Count.put(aaSet.getCompoundForString("C"), (int) (numC + smallC));
		extinctAA2Count.put(aaSet.getCompoundForString("Y"), numY + smallY);
		return extinctAA2Count;
	}

	@Override
	public double getInstabilityIndex(ProteinSequence sequence) {
		double sum = 0.0;
		String s = sequence.getSequenceAsString().toUpperCase();
		for(int i = 0; i < sequence.getLength() - 1; i++){
			String dipeptide = s.substring(i, i+2);
			if(Constraints.diAA2Instability.containsKey(dipeptide)){
				sum += Constraints.diAA2Instability.get(dipeptide);
			}
		}
		return sum * 10 / (s.length() - Utils.getNumberOfInvalidChar(s, null, true));
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
		final double a = 2.9;
		final double b = 3.9;
		double xAla = aa2Composition.get(aaSet.getCompoundForString("A"));
		double xVal = aa2Composition.get(aaSet.getCompoundForString("V"));
		double xIle = aa2Composition.get(aaSet.getCompoundForString("I"));
		double xLeu = aa2Composition.get(aaSet.getCompoundForString("L"));
		return (xAla + (a * xVal) + (b * (xIle + xLeu))) * 100;
	}

	@Override
	public double getAvgHydropathy(ProteinSequence sequence) {
		int validLength = 0;
		double total = 0.0;
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		char[] seq = this.getSequence(sequence.toString(), true);
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
	public double getIsoelectricPoint(ProteinSequence sequence, boolean useExpasyValues) {
		if(useExpasyValues){
			return this.getIsoelectricPointExpasy(sequence.toString().toUpperCase());
		}else{
			return this.getIsoelectricPointInnovagen(sequence);
		}
	}
	
	private double getIsoelectricPointInnovagen(ProteinSequence sequence){
		double currentPH = 7.0;
		double changeSize = 7.0;
		String sequenceString = sequence.toString();
		char nTerminalChar = sequenceString.charAt(0);
		char cTerminalChar = sequenceString.charAt(sequenceString.length() - 1);
		
		Map<AminoAcidCompound, Integer> chargedAA2Count = this.getChargedAACount(sequence);
		double margin;
		final double difference = 0.0001;
		
		while(true){
			margin = this.getNetChargeInnovagen(chargedAA2Count, currentPH, nTerminalChar, cTerminalChar);
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
	
	/*
	 *  Pseudo code obtained from email correspondance with ExPASy Helpdesk, Gregoire Rossier
	 */
	//
	// Table of pk values :
	// Note: the current algorithm does not use the last two columns.
	// Each row corresponds to an amino acid starting with Ala. J, O and U are
	// inexistant, but here only in order to have the complete alphabet.
	//
	// Ct Nt Sm Sc Sn
	//
	private final double cPk[][] = {
			{3.55, 7.59, 0.0},  // A
			{3.55, 7.50, 0.0},  // B
			{3.55, 7.50, 9.00}, // C
//			{4.55, 7.50, 4.05}, // D
//			{4.75, 7.70, 4.45}, // E
			{3.55, 7.50, 4.05}, // D
			{3.55, 7.70, 4.45}, // E
			{3.55, 7.50, 0}, // F
			{3.55, 7.50, 0}, // G
			{3.55, 7.50, 5.98}, // H
			{3.55, 7.50, 0.0}, // I
			{0.0, 0.0, 0.0}, // J
			{3.55, 7.50, 10.00}, // K
			{3.55, 7.50, 0.0}, // L
			{3.55, 7.00, 0.0},// M
			{3.55, 7.50, 0.0},// N
			{0.00, 0.00, 0.0},// O
			{3.55, 8.36, 0.0},// P
			{3.55, 7.50, 0.0}, // Q
			{3.55, 7.50, 12.0},// R
			{3.55, 6.93, 0.0},// S
			{3.55, 6.82, 0.0}, // T
			{0.00, 0.00, 0.0}, // U
			{3.55, 7.44, 0.0},// V
			{3.55, 7.50, 0.0},// W
			{3.55, 7.50, 0.0},// X
			{3.55, 7.50, 10.00},// Y
			{3.55, 7.50, 0.0}}; // Z

	private final double PH_MIN = 0.0; /* minimum pH value */
	private final double PH_MAX = 14.0; /* maximum pH value */
	private final double MAXLOOP = 2000.0; /* maximum number of iterations */
	private final double EPSI = 0.0001; /* desired precision */
	
	private double exp10(double pka){
		return Math.pow(10, pka);
	}
	
	private double getIsoelectricPointExpasy(String sequence){
		//
		// Compute the amino-acid composition.
		//
		int comp[] = new int[26];
		for(int i = 0; i < sequence.length(); i++){
			int index = sequence.charAt(i) - 'A';
			if(index < 0 || index >= 26) continue;
			comp[index]++;
		}
		//
		// Look up N-terminal and C-terminal residue.
		//
		int nTermResidue = -1;
		int index = 0;
		while((nTermResidue < 0 || nTermResidue >= 26) && index < 25){
			nTermResidue = sequence.charAt(index++) - 'A';
		}
		
		int cTermResidue = -1;
		index = 1;
		while((cTermResidue < 0 || cTermResidue >= 26) && index < 25){
			cTermResidue = sequence.charAt(sequence.length() - index++) - 'A';
		}

		double phMin = PH_MIN;
		double phMax = PH_MAX;
		
		double phMid = 0.0;
		double charge = 1.0;
		for (int i = 0; i < MAXLOOP && (phMax - phMin) > EPSI; i++){
			phMid = phMin + (phMax - phMin) / 2.0;
	
			charge = getNetChargeExpasy(comp, nTermResidue, cTermResidue, phMid);
	
			if (charge > 0.0) phMin = phMid;
			else phMax = phMid;
		}
		return phMid;
	}
	
	@Override
	public double getIsoelectricPoint(ProteinSequence sequence){
		return getIsoelectricPoint(sequence, true);
	}

	@Override
	public double getNetCharge(ProteinSequence sequence) {
		return getNetCharge(sequence, true);
	}
	
	@Override
	public double getNetCharge(ProteinSequence sequence, boolean useExpasyValues){
		return getNetCharge(sequence, true, 7.0);
	}
	
	@Override
	public double getNetCharge(ProteinSequence sequence, boolean useExpasyValues, double pHPoint){
		if(useExpasyValues){
			return getNetChargeExpasy(sequence.toString().toUpperCase(), pHPoint);
		}else{
			return getNetChargeInnovagen(sequence, pHPoint);
		}
	}
	
	private double getNetChargeExpasy(String sequence, double pHPoint){
		//
		// Compute the amino-acid composition.
		//
		int comp[] = new int[26];
		for(int i = 0; i < sequence.length(); i++){
			int index = sequence.charAt(i) - 'A';
			if(index < 0 || index >= 26) continue;
			comp[index]++;
		}
		//
		// Look up N-terminal and C-terminal residue.
		//
		int nTermResidue = sequence.charAt(0) - 'A';
		int cTermResidue = sequence.charAt(sequence.length() - 1) - 'A';
		return getNetChargeExpasy(comp, nTermResidue, cTermResidue, pHPoint);
	}
	
	private double getNetChargeExpasy(int comp[], int nTermResidue, int cTermResidue, double ph){
		double cter = 0.0; 
		if(cTermResidue >= 0 && cTermResidue < 26) cter = exp10(-cPk[cTermResidue][0]) / (exp10(-cPk[cTermResidue][0]) + exp10(-ph));
		double nter = 0.0; 
		if(nTermResidue >= 0 && nTermResidue < 26) nter = exp10(-ph) / (exp10(-cPk[nTermResidue][1]) + exp10(-ph));

		double carg = comp['R' - 'A'] * exp10(-ph) / (exp10(-cPk['R' - 'A'][2]) + exp10(-ph)); 
		double chis = comp['H' - 'A'] * exp10(-ph) / (exp10(-cPk['H' - 'A'][2]) + exp10(-ph));
		double clys = comp['K' - 'A'] * exp10(-ph) / (exp10(-cPk['K' - 'A'][2]) + exp10(-ph));

		double casp = comp['D' - 'A'] * exp10(-cPk['D' - 'A'][2]) / (exp10(-cPk['D' - 'A'][2]) + exp10(-ph));
		double cglu = comp['E' - 'A'] * exp10(-cPk['E' - 'A'][2]) / (exp10(-cPk['E' - 'A'][2]) + exp10(-ph));

		double ccys = comp['C' - 'A'] * exp10(-cPk['C' - 'A'][2]) / (exp10(-cPk['C' - 'A'][2]) + exp10(-ph));
		double ctyr = comp['Y' - 'A'] * exp10(-cPk['Y' - 'A'][2]) / (exp10(-cPk['Y' - 'A'][2]) + exp10(-ph));

		return (carg + clys + chis + nter) - (casp + cglu + ctyr + ccys + cter);
	}
	
	private double getNetChargeInnovagen(ProteinSequence sequence, double pHPoint) {
		Map<AminoAcidCompound, Integer> chargedAA2Count = this.getChargedAACount(sequence);
		String sequenceString = sequence.getSequenceAsString();
		return getNetChargeInnovagen(chargedAA2Count, pHPoint, sequenceString.charAt(0), sequenceString.charAt(sequenceString.length() - 1));
	}
	
	private double getNetChargeInnovagen(Map<AminoAcidCompound, Integer> chargedAA2Count, double ph, char nTerminalChar, char cTerminalChar){
		//Constraints.aa2PKa is aleady reinitialized in getChargedAACount hence no need to do it again
		
		//Lys => K, Arg => R, His => H
		//Asp => D, Glu => E, Cys => C, Tyr => Y
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		
		double nTerminalCharge = 0.0; 
		AminoAcidCompound nTermCompound = aaSet.getCompoundForString(nTerminalChar + "");
		if(Constraints.aa2NTerminalPka.containsKey(nTermCompound)){
			nTerminalCharge = this.getPosCharge(Constraints.aa2NTerminalPka.get(nTermCompound), ph);
		}			
		
		double cTerminalCharge = 0.0;
		AminoAcidCompound cTermCompound = aaSet.getCompoundForString(cTerminalChar + "");
		if(Constraints.aa2CTerminalPka.containsKey(cTermCompound)){
			cTerminalCharge = this.getNegCharge(Constraints.aa2CTerminalPka.get(cTermCompound), ph);
		}
		
		double kCharge = chargedAA2Count.get(aaSet.getCompoundForString("K")) * this.getPosCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("K")), ph);
		double rCharge = chargedAA2Count.get(aaSet.getCompoundForString("R")) * this.getPosCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("R")), ph);
		double hCharge = chargedAA2Count.get(aaSet.getCompoundForString("H")) * this.getPosCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("H")), ph);
		double dCharge = chargedAA2Count.get(aaSet.getCompoundForString("D")) * this.getNegCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("D")), ph);
		double eCharge = chargedAA2Count.get(aaSet.getCompoundForString("E")) * this.getNegCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("E")), ph);
		double cCharge = chargedAA2Count.get(aaSet.getCompoundForString("C")) * this.getNegCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("C")), ph);
		double yCharge = chargedAA2Count.get(aaSet.getCompoundForString("Y")) * this.getNegCharge(Constraints.aa2PKa.get(aaSet.getCompoundForString("Y")), ph);
//		if((kCharge + rCharge + hCharge) == 0.0 && (dCharge + eCharge + cCharge + yCharge) == 0.0){
//			return 0.0;
//		}
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
		char[] seq = this.getSequence(sequence.getSequenceAsString(), true);
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
	public double getEnrichment(ProteinSequence sequence, AminoAcidCompound aminoAcidCode) {
		double counter = 0.0;
		char[] seq = this.getSequence(sequence.getSequenceAsString(), true);
		for(char aa:seq){
			if(aminoAcidCode.getShortName().equals(aa + "")){
				counter++;
			}
		}
		return counter/sequence.getLength();
	}

	@Override
	public Map<AminoAcidCompound, Double> getAAComposition(ProteinSequence sequence) {
		int validLength = 0;
		Map<AminoAcidCompound, Double> aa2Composition = new HashMap<AminoAcidCompound, Double>();
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		for(AminoAcidCompound aa:aaSet.getAllCompounds()){
			aa2Composition.put(aa, 0.0);
		}
		char[] seq = this.getSequence(sequence.toString(), true);
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
}
