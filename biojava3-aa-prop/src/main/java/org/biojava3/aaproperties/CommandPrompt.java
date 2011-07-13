package org.biojava3.aaproperties;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.xml.bind.JAXBException;

import org.biojava3.aaproperties.xml.AminoAcidCompositionTable;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;

public class CommandPrompt {
	public static void main(String[] args) throws IOException, JAXBException{
		Set<String> argumentsSet = new HashSet<String>();
		Character[] arguments = new Character[16];
		List<Character> propertyList = new ArrayList<Character>();
		List<Character> specificList = new ArrayList<Character>();
		String inputLocation = null;
		String outputLocation = null;
		String aminoAcidCompositionLocation = null;
		String elementMassLocation = null;
		String delimiter = ",";
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
		arguments[5] = 'a';
		arguments[6] = '1';
		arguments[7] = '2';
		arguments[8] = '3';
		arguments[9] = '4';
		arguments[10] = '5';
		arguments[11] = '6';
		arguments[12] = '7';
		arguments[13] = '8';
		arguments[14] = '9';
		arguments[15] = '0';
		for(char c:arguments){
			argumentsSet.add("-" + c);
		}
		for(int i = 0; i < args.length; i++){
			if(args[i].charAt(0) != '-' || args[i].length() != 2){
				showOptions();
				throw new Error("Unknown option: " + args[i]);
			}else{
				switch(args[i].charAt(1)){
				//Required
				case 'i': inputLocation = args[++i]; break;
				case 'o': outputLocation = args[++i]; break;
				//Optional
				case 'f':
					i++;
					if(args[i].equalsIgnoreCase("csv")) delimiter = ",";
					else if(args[i].equalsIgnoreCase("tsv")) delimiter = "\t";
					else throw new Error("Invalid value for -f: " + args[i] + ". Please choose either csv or tsv only.");
					break;
				case 'x': aminoAcidCompositionLocation = args[++i]; break;
				case 'y': elementMassLocation = args[++i]; break;
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
					showOptions();
					throw new Error("Unknown option: " + args[i]);
				}
			}
		}

		if(inputLocation == null) {
			showOptions();
			throw new Error("Please do provide location of input file.");
		}
		if(outputLocation == null){
			showOptions();
			throw new Error("Please do provide location of output file.");
		}
		if(propertyList.size() == 0){
			showOptions();
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
			
		BufferedReader input = new BufferedReader(new FileReader(inputLocation));
		BufferedWriter output = new BufferedWriter(new FileWriter(outputLocation));
		printHeader(output, propertyList, specificList, delimiter);
		String line;
		String header = null;
		String sequence = "";
		while((line = input.readLine()) != null){
			if(line.charAt(0) == '>'){
				//header line
				if(sequence.length() > 0){
					compute(output, header, sequence.trim(), delimiter, aaTable, propertyList, specificList);
					sequence = "";
				}
				header = line;
			}else{
				//sequence line
				sequence += line;
			}
		}
		//Need for the last sequence
		compute(output, header, sequence.trim(), delimiter, aaTable, propertyList, specificList);
		output.close();
		input.close();
	}
	
	private static void printHeader(BufferedWriter output, List<Character> propertyList, List<Character> specificList, String delimiter) throws IOException{
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
		for(int i = 0; i < propertyList.size(); i++){
			Character c = propertyList.get(i);
			if(i != 0) output.write(delimiter);
			switch(c){
			case '1': output.write("MolecularWeight"); break;
			case '2': output.write("Absorbance(true)" + delimiter + "Absorbance(false)"); break;
			case '3': output.write("ExtinctionCoefficient(true)" + delimiter + "ExtinctionCoefficient(false)"); break;
			case '4': output.write("InstabilityIndex"); break;
			case '5': output.write("ApliphaticIndex"); break;
			case '6': output.write("AverageHydropathyValue"); break;
			case '7': output.write("IsoelectricPoint"); break;
			case '8': output.write("NetCharge@pH7"); break;
			case '9': output.write("Composition_A" + delimiter + "Composition_R" + delimiter + "Composition_N" + delimiter + "Composition_D" + 
					delimiter + "Composition_C" + delimiter + "Composition_E" + delimiter + "Composition_Q" + delimiter + "Composition_G" + 
					delimiter + "Composition_H" + delimiter + "Composition_I" + delimiter + "Composition_L" + delimiter + "Composition_K" + 
					delimiter + "Composition_M" + delimiter + "Composition_F" + delimiter + "Composition_P" + delimiter + "Composition_S" + 
					delimiter + "Composition_T" + delimiter + "Composition_W" + delimiter + "Composition_Y" + delimiter + "Composition_V");
				break;
			case '0': output.write("Composition_" + specificList.get(specificCount++)); break;
			}
		}
		output.newLine();
		output.flush();
	}

	private static void compute(BufferedWriter output, String header, String sequence, String delimiter, 
			AminoAcidCompositionTable aaTable, List<Character> propertyList, List<Character> specificList) throws IOException{
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
		sequence = Utils.checkSequence(sequence);
		ProteinSequence pSequence = new ProteinSequence(sequence);
		IPeptideProperties pp = new PeptidePropertiesImpl();
		AminoAcidCompoundSet aaSet = new AminoAcidCompoundSet();
		int specificCount = 0;
		for(int i = 0; i < propertyList.size(); i++){
			if(i != 0) output.write(delimiter);
			Character c = propertyList.get(i);
			switch(c){
			case '1': 
				if(aaTable == null) output.write("" + pp.getMolecularWeight(pSequence));
				else output.write("" + pp.getMolecularWeight(pSequence));
				break;
			case '2': output.write(pp.getAbsorbance(pSequence, true) + delimiter + pp.getAbsorbance(pSequence, false)); break;
			case '3': output.write(pp.getExtinctionCoefficient(pSequence, true) + delimiter + pp.getExtinctionCoefficient(pSequence, false)); break;
			case '4': output.write("" + pp.getInstabilityIndex(pSequence)); break;
			case '5': output.write("" + pp.getApliphaticIndex(pSequence)); break;
			case '6': output.write("" + pp.getAvgHydropathy(pSequence)); break;
			case '7': output.write("" + pp.getIsoelectricPoint(pSequence)); break;
			case '8': output.write("" + pp.getNetCharge(pSequence)); break;
			case '9': 
				Map<AminoAcidCompound, Double> aaCompound2Double = pp.getAAComposition(pSequence);
				//(A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V)
				output.write(aaCompound2Double.get(Constraints.A) + delimiter);
				output.write(aaCompound2Double.get(Constraints.R) + delimiter);
				output.write(aaCompound2Double.get(Constraints.N) + delimiter);
				output.write(aaCompound2Double.get(Constraints.D) + delimiter);
				output.write(aaCompound2Double.get(Constraints.C) + delimiter);
				output.write(aaCompound2Double.get(Constraints.E) + delimiter);
				output.write(aaCompound2Double.get(Constraints.Q) + delimiter);
				output.write(aaCompound2Double.get(Constraints.G) + delimiter);
				output.write(aaCompound2Double.get(Constraints.H) + delimiter);
				output.write(aaCompound2Double.get(Constraints.I) + delimiter);
				output.write(aaCompound2Double.get(Constraints.L) + delimiter);
				output.write(aaCompound2Double.get(Constraints.K) + delimiter);
				output.write(aaCompound2Double.get(Constraints.M) + delimiter);
				output.write(aaCompound2Double.get(Constraints.F) + delimiter);
				output.write(aaCompound2Double.get(Constraints.P) + delimiter);
				output.write(aaCompound2Double.get(Constraints.S) + delimiter);
				output.write(aaCompound2Double.get(Constraints.T) + delimiter);
				output.write(aaCompound2Double.get(Constraints.W) + delimiter);
				output.write(aaCompound2Double.get(Constraints.Y) + delimiter);
				output.write(aaCompound2Double.get(Constraints.V) + "");
				break;
			case '0': output.write(pp.getAAComposition(pSequence).get(aaSet.getCompoundForString(specificList.get(specificCount++) + "")) + ""); break;
			}
		}
		output.newLine();
		output.flush();
	}

	private static void showOptions(){
		System.err.println("Required");
		System.err.println("-i location of input FASTA file");
		System.err.println("-o location of output file");
		System.err.println();

		System.err.println("Optional");
		System.err.println("-f output format [csv (default) or tsv]");
		System.err.println("-x location of Amino Acid Composition XML file for defining amino acid composition");
		System.err.println("-y location of Element Mass XML file for defining mass of elements");
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
