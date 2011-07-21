package org.biojava3.aaproperties;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.biojava3.aaproperties.xml.AminoAcidCompositionTable;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojava3.core.sequence.io.ProteinSequenceCreator;
import org.biojava3.core.sequence.template.CompoundSet;

public class CommandPrompt {
	public static void main(String[] args) throws Exception{
		run(args);
	}
	public static void run(String[] args) throws Exception{
		/*
		 * Initialization
		 */
		Set<String> argumentsSet = new HashSet<String>();
		Character[] arguments = new Character[17];
		List<Character> propertyList = new ArrayList<Character>();
		List<Character> specificList = new ArrayList<Character>();
		String inputLocation = null;
		String outputLocation = null;
		String aminoAcidCompositionLocation = null;
		String elementMassLocation = null;
		String delimiter = ",";
		int decimalPlace = 4;
		/*
		 * a properties of 1-9
		 * 1 Molecular weight
		 * 2 Absorbance
		 * 3 Extinction coefficient
		 * 4 Instability index
		 * 5 Apliphatic index
		 * 6 Average hydropathy value
		 * 7 Isoelectric point
		 * 8 Net charge at pH 7
		 * 9 Composition of the 20 standard amino acid
		 * 0 Composition of the specific amino acid
		 */
		arguments[0] = 'i';
		arguments[1] = 'o';
		arguments[2] = 'f';
		arguments[3] = 'x';
		arguments[4] = 'y';
		arguments[5] = 'd';
		arguments[6] = 'a';
		arguments[7] = '1';
		arguments[8] = '2';
		arguments[9] = '3';
		arguments[10] = '4';
		arguments[11] = '5';
		arguments[12] = '6';
		arguments[13] = '7';
		arguments[14] = '8';
		arguments[15] = '9';
		arguments[16] = '0';
		for(char c:arguments){
			argumentsSet.add("-" + c);
		}
		
		/*
		 * Parse input arguments
		 */
		for(int i = 0; i < args.length; i++){
			if(args[i].charAt(0) != '-' || args[i].length() != 2){
				showHelp();
				throw new Error("Unknown option: " + args[i]);
			}else{
				switch(args[i].charAt(1)){
				//Required
				case 'i': inputLocation = args[++i]; break;
				//Optional
				case 'o': outputLocation = args[++i]; break;
				case 'f':
					i++;
					if(args[i].equalsIgnoreCase("csv")) delimiter = ",";
					else if(args[i].equalsIgnoreCase("tsv")) delimiter = "\t";
					else throw new Error("Invalid value for -f: " + args[i] + ". Please choose either csv or tsv only.");
					break;
				case 'x': aminoAcidCompositionLocation = args[++i]; break;
				case 'y': elementMassLocation = args[++i]; break;
				case 'd': decimalPlace = Integer.parseInt(args[++i]); break;
				//Properties
				case 'a':
					propertyList.add('1');
					propertyList.add('2');
					propertyList.add('3');
					propertyList.add('4');
					propertyList.add('5');
					propertyList.add('6');
					propertyList.add('7');
					propertyList.add('8');
					propertyList.add('9');
					break;
				case '1': propertyList.add('1'); break;
				case '2': propertyList.add('2'); break;
				case '3': propertyList.add('3'); break;
				case '4': propertyList.add('4'); break;
				case '5': propertyList.add('5'); break;
				case '6': propertyList.add('6'); break;
				case '7': propertyList.add('7'); break;
				case '8': propertyList.add('8'); break;
				case '9': propertyList.add('9'); break;
				case '0': 
					propertyList.add('0'); 
					i++;
					if(args[i].length() != 1) throw new Error("Invalid value: " + args[i] + ". Amino Acid Symbol should be of single character");
					specificList.add(args[i].toUpperCase().charAt(0)); 
					break;
				default: 
					showHelp();
					throw new Error("Unknown option: " + args[i]);
				}
			}
		}

		/*
		 * Check for validity of input arguments
		 */
		if(inputLocation == null) {
			showHelp();
			throw new Error("Please do provide location of input file.");
		}
		if(propertyList.size() == 0){
			showHelp();
			throw new Error("Please at least specify a property to compute.");
		}
		AminoAcidCompositionTable aaTable = null;
		if(aminoAcidCompositionLocation != null && elementMassLocation == null){
			aaTable = PeptideProperties.obtainAminoAcidCompositionTable(new File(aminoAcidCompositionLocation));
		}else if(aminoAcidCompositionLocation != null && elementMassLocation != null){
			aaTable = PeptideProperties.obtainAminoAcidCompositionTable(new File(aminoAcidCompositionLocation, elementMassLocation));
		}else if(aminoAcidCompositionLocation == null && elementMassLocation != null){
			throw new Error("You have define the location of Element Mass XML file. Please also define the location of Amino Acid Composition XML file");
		}
		
		/*
		 * Read input file and generate output
		 */
		PrintStream output;
		if(outputLocation != null)
			output = new PrintStream(new File(outputLocation));
		else
			output = System.out;
		printHeader(output, propertyList, specificList, delimiter);
		LinkedHashMap<String, ProteinSequence> a = readInputFile(inputLocation, aaTable);
		//Need for the last sequence
		for(Entry<String, ProteinSequence> entry:a.entrySet()){
			compute(output, entry.getValue().getOriginalHeader(), entry.getValue().getSequenceAsString().trim(), delimiter, aaTable, propertyList, specificList, decimalPlace);
		}
		output.close();
	}
	
	private static LinkedHashMap<String, ProteinSequence> readInputFile(String inputLocation, AminoAcidCompositionTable aaTable) throws Exception{
		FileInputStream inStream = new FileInputStream(inputLocation);
		FastaReader<ProteinSequence,AminoAcidCompound> fastaReader;
		if(aaTable == null)
			fastaReader = new FastaReader<ProteinSequence,AminoAcidCompound>(
					inStream, new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(), 
					new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		else
			fastaReader = new FastaReader<ProteinSequence,AminoAcidCompound>(
					inStream, new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(), 
					new ProteinSequenceCreator(aaTable.getAminoAcidCompoundSet()));
		return fastaReader.process();
	}
	
	private static void printHeader(PrintStream output, List<Character> propertyList, List<Character> specificList, String delimiter) throws IOException{
		int specificCount = 0;
		/* 
		 * 1 Molecular weight
		 * 2 Absorbance (assumed Cys reduced and assume Cys to form cystines)
		 * 3 Extinction coefficient (assumed Cys reduced and assume Cys to form cystines)
		 * 4 Instability index
		 * 5 Apliphatic index
		 * 6 Average hydropathy value
		 * 7 Isoelectric point
		 * 8 Net charge at pH 7
		 * 9 Composition of the 20 standard amino acid
		 * 0 Composition of the specific amino acid
		 */
		List<String> sList = new ArrayList<String>();
		for(Character c:propertyList){
			switch(c){
			case '1': sList.add("MolecularWeight"); break;
			case '2': 
				sList.add("Absorbance(true)");
				sList.add("Absorbance(false)"); 
				break;
			case '3': 
				sList.add("ExtinctionCoefficient(true)");
				sList.add("ExtinctionCoefficient(false)"); 
				break;
			case '4': sList.add("InstabilityIndex"); break;
			case '5': sList.add("ApliphaticIndex"); break;
			case '6': sList.add("AverageHydropathyValue"); break;
			case '7': sList.add("IsoelectricPoint"); break;
			case '8': sList.add("NetCharge@pH7"); break;
			case '9': 
				sList.add("Composition_A");
				sList.add("Composition_R");
				sList.add("Composition_N");
				sList.add("Composition_D"); 
				sList.add("Composition_C");
				sList.add("Composition_E");
				sList.add("Composition_Q");
				sList.add("Composition_G");
				sList.add("Composition_H");
				sList.add("Composition_I");
				sList.add("Composition_L");
				sList.add("Composition_K");
				sList.add("Composition_M");
				sList.add("Composition_F");
				sList.add("Composition_P");
				sList.add("Composition_S");
				sList.add("Composition_T");
				sList.add("Composition_W");
				sList.add("Composition_Y");
				sList.add("Composition_V");
				break;
			case '0': sList.add("Composition_" + specificList.get(specificCount++)); break;
			}
		}
		for(int i = 0; i < sList.size(); i++){
			if(i != 0) output.print(delimiter);
			output.print(sList.get(i));
		}
		output.println();
		output.flush();
	}

	private static void compute(PrintStream output, String header, String sequence, String delimiter, 
			AminoAcidCompositionTable aaTable, List<Character> propertyList, List<Character> specificList, int decimalPlace) throws IOException{
		/* 
		 * 1 Molecular weight
		 * 2 Absorbance (assumed Cys reduced and assume Cys to form cystines)
		 * 3 Extinction coefficient
		 * 4 Instability index
		 * 5 Apliphatic index
		 * 6 Average hydropathy value
		 * 7 Isoelectric point
		 * 8 Net charge at pH 7
		 * 9 Composition of the 20 standard amino acid
		 * 0 Composition of the specific amino acid
		 */
		ProteinSequence pSequence;
		CompoundSet<AminoAcidCompound> aaSet;
		if(aaTable != null){
			sequence = Utils.checkSequence(sequence, aaTable.getSymbolSet());
			pSequence = new ProteinSequence(sequence, aaTable.getAminoAcidCompoundSet());
			aaSet = aaTable.getAminoAcidCompoundSet(); 
		}else{
			sequence = Utils.checkSequence(sequence);
			pSequence = new ProteinSequence(sequence);
			aaSet = AminoAcidCompoundSet.getAminoAcidCompoundSet();
		}
		IPeptideProperties pp = new PeptidePropertiesImpl();
		
		int specificCount = 0;
		List<Double> dList = new ArrayList<Double>();
		for(Character c:propertyList){
			switch(c){
			case '1': 
				if(aaTable == null) 
					dList.add(pp.getMolecularWeight(pSequence));
				else 
					dList.add(pp.getMolecularWeight(pSequence));
				break;
			case '2': 
				dList.add(pp.getAbsorbance(pSequence, true)); 
				dList.add(pp.getAbsorbance(pSequence, false)); 
				break;
			case '3': 
				dList.add(pp.getExtinctionCoefficient(pSequence, true));
				dList.add(pp.getExtinctionCoefficient(pSequence, false)); 
				break;
			case '4': dList.add(pp.getInstabilityIndex(pSequence)); break;
			case '5': dList.add(pp.getApliphaticIndex(pSequence)); break;
			case '6': dList.add(pp.getAvgHydropathy(pSequence)); break;
			case '7': dList.add(pp.getIsoelectricPoint(pSequence)); break;
			case '8': dList.add(pp.getNetCharge(pSequence)); break;
			case '9': 
				Map<AminoAcidCompound, Double> aaCompound2Double = pp.getAAComposition(pSequence);
				//(A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V)
				dList.add(aaCompound2Double.get(Constraints.A));
				dList.add(aaCompound2Double.get(Constraints.R));
				dList.add(aaCompound2Double.get(Constraints.N));
				dList.add(aaCompound2Double.get(Constraints.D));
				dList.add(aaCompound2Double.get(Constraints.C));
				dList.add(aaCompound2Double.get(Constraints.E));
				dList.add(aaCompound2Double.get(Constraints.Q));
				dList.add(aaCompound2Double.get(Constraints.G));
				dList.add(aaCompound2Double.get(Constraints.H));
				dList.add(aaCompound2Double.get(Constraints.I));
				dList.add(aaCompound2Double.get(Constraints.L));
				dList.add(aaCompound2Double.get(Constraints.K));
				dList.add(aaCompound2Double.get(Constraints.M));
				dList.add(aaCompound2Double.get(Constraints.F));
				dList.add(aaCompound2Double.get(Constraints.P));
				dList.add(aaCompound2Double.get(Constraints.S));
				dList.add(aaCompound2Double.get(Constraints.T));
				dList.add(aaCompound2Double.get(Constraints.W));
				dList.add(aaCompound2Double.get(Constraints.Y));
				dList.add(aaCompound2Double.get(Constraints.V));
				break;
			case '0': dList.add(pp.getEnrichment(pSequence, aaSet.getCompoundForString("" + specificList.get(specificCount++)))); break;
			}
		}
		for(int i = 0; i < dList.size(); i++){
			if(i != 0) output.print(delimiter);
			output.print(Utils.roundToDecimals(dList.get(i), decimalPlace));
		}
		output.println();
		output.flush();
	}

	private static void showHelp(){
		System.err.println("Examples");
		System.err.println("Example 1 (Generates all possible properties): java -jar AAProperties.jar -i test.fasta -a");
		System.err.println("Example 2 (Generates only molecular weight, extinction coefficient and isoelectric point): " +
				"java -jar AAProperties.jar -i test.fasta -1 -3 -7");
		System.err.println("Example 2 (Generates composition of two specific amino acid symbol and molecular weight): " +
			"java -jar AAProperties.jar -i test.fasta -0 A -0 N -1");
		System.err.println();
		
		System.err.println("Required");
		System.err.println("-i location of input FASTA file");
		System.err.println();

		System.err.println("Optional");
		System.err.println("-o location of output file [standard output (default)]");
		System.err.println("-f output format [csv (default) or tsv]");
		System.err.println("-x location of Amino Acid Composition XML file for defining amino acid composition");
		System.err.println("-y location of Element Mass XML file for defining mass of elements");
		System.err.println("-d number of decimals (int) [4 (default)]");
		System.err.println();

		System.err.println("Provide at least one of them");
		System.err.println("-a compute properties of option 1-9");
		System.err.println("-1 compute molecular weight");
		System.err.println("-2 compute absorbance");
		System.err.println("-3 compute extinction coefficient");
		System.err.println("-4 compute instability index");
		System.err.println("-5 compute apliphatic index");
		System.err.println("-6 compute average hydropathy value");
		System.err.println("-7 compute isoelectric point");
		System.err.println("-8 compute net charge at pH 7");
		System.err.println("-9 compute composition of 20 standard amino acid (A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V)");
		System.err.println("-0 compute composition of specific amino acid symbol");
		System.err.println();
	}
}
