package org.biojava3.aaproperties;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.biojava3.aaproperties.xml.AminoAcidCompositionTable;
import org.biojava3.aaproperties.xml.CaseFreeAminoAcidCompoundSet;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojava3.core.sequence.io.ProteinSequenceCreator;
import org.biojava3.core.sequence.template.CompoundSet;

public class CommandPrompt {
	
	/**
	 * The main method
	 * @param args
	 * 	See showHelp for a list of available arguments
	 * @throws Exception
	 *  To handle exception thrown by reading of XML files
	 */
	public static void main(String[] args) throws Exception{
		run(args);
	}
	
	private static AminoAcidCompositionTable checkForValidityAndObtainAATable(String inputLocation, int propertyListSize, String aminoAcidCompositionLocation,
			String elementMassLocation) throws Exception{
		if(inputLocation == null) {
			showHelp();
			throw new Error("Please do provide location of input file.");
		}
		if(propertyListSize == 0){
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
		return aaTable;
	}
	
	private static void readInputAndGenerateOutput(String outputLocation, List<Character> propertyList, List<Character> specificList,
			String delimiter, String inputLocation, AminoAcidCompositionTable aaTable, int decimalPlace) throws Exception{
		PrintStream output;
		if(outputLocation != null)
			output = new PrintStream(new File(outputLocation));
		else
			output = System.out;
		printHeader(output, propertyList, specificList, delimiter);
		LinkedHashMap<String, ProteinSequence> a = readInputFile(inputLocation, aaTable);
		//Need for the last sequence
		for(Entry<String, ProteinSequence> entry:a.entrySet()){
			compute(output, entry.getValue().getOriginalHeader(), entry.getValue().getSequenceAsString().trim(), delimiter, aaTable, propertyList, specificList,
					decimalPlace);
		}
		output.close();
	}
	
	public static void run(String[] args) throws Exception{
		/*
		 * Parse input arguments
		 */
		List<Character> propertyList = new ArrayList<Character>();
		List<Character> specificList = new ArrayList<Character>();
		String inputLocation = null;
		String outputLocation = null;
		String aminoAcidCompositionLocation = null;
		String elementMassLocation = null;
		String delimiter = ",";
		int decimalPlace = 4;
		
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
		AminoAcidCompositionTable aaTable = checkForValidityAndObtainAATable(inputLocation, propertyList.size(), aminoAcidCompositionLocation, 
				elementMassLocation);
		
		/*
		 * Read input file and generate output
		 */
		readInputAndGenerateOutput(outputLocation, propertyList, specificList, delimiter, inputLocation, aaTable, decimalPlace);
	}
	
	private static LinkedHashMap<String, ProteinSequence> readInputFile(String inputLocation, AminoAcidCompositionTable aaTable) throws Exception{
		FileInputStream inStream = new FileInputStream(inputLocation);
		FastaReader<ProteinSequence, AminoAcidCompound> fastaReader;
		if(aaTable == null){
			CompoundSet<AminoAcidCompound>	set = CaseFreeAminoAcidCompoundSet.getAminoAcidCompoundSet();
			fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
					inStream, new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(), 
					new ProteinSequenceCreator(set));
		}else{
			fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
					inStream, new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(), 
					new ProteinSequenceCreator(aaTable.getAminoAcidCompoundSet()));
		}
		return fastaReader.process();
	}
	
	public enum PropertyName{MolecularWeight, Absorbance_True, Absorbance_False, ExtinctionCoefficient_True, ExtinctionCoefficient_False, 
		InstabilityIndex, ApliphaticIndex, AverageHydropathyValue, IsoelectricPoint, NetCharge_pH_7, A, R, 
		N, D, C, E, Q, G, H, I, L,
		K, M, F, P, S, T, W, Y, V};
	
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
		sList.add("SequenceName");
		for(Character c:propertyList){
			switch(c){
			case '1': sList.add(PropertyName.MolecularWeight.toString()); break;
			case '2': sList.add(PropertyName.Absorbance_True.toString()); sList.add(PropertyName.Absorbance_False.toString()); break;
			case '3': sList.add(PropertyName.ExtinctionCoefficient_True.toString()); sList.add(PropertyName.ExtinctionCoefficient_False.toString()); break;
			case '4': sList.add(PropertyName.InstabilityIndex.toString()); break;
			case '5': sList.add(PropertyName.ApliphaticIndex.toString()); break;
			case '6': sList.add(PropertyName.AverageHydropathyValue.toString()); break;
			case '7': sList.add(PropertyName.IsoelectricPoint.toString()); break;
			case '8': sList.add(PropertyName.NetCharge_pH_7.toString()); break;
			case '9': 
				sList.add(PropertyName.A.toString()); sList.add(PropertyName.R.toString()); 
				sList.add(PropertyName.N.toString()); sList.add(PropertyName.D.toString()); 
				sList.add(PropertyName.C.toString()); sList.add(PropertyName.E.toString());
				sList.add(PropertyName.Q.toString()); sList.add(PropertyName.G.toString());
				sList.add(PropertyName.H.toString()); sList.add(PropertyName.I.toString());
				sList.add(PropertyName.L.toString()); sList.add(PropertyName.K.toString());
				sList.add(PropertyName.M.toString()); sList.add(PropertyName.F.toString());
				sList.add(PropertyName.P.toString()); sList.add(PropertyName.S.toString());
				sList.add(PropertyName.T.toString()); sList.add(PropertyName.W.toString());
				sList.add(PropertyName.Y.toString()); sList.add(PropertyName.V.toString());
				break;
			case '0': sList.add("" + specificList.get(specificCount++)); break;
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
			AminoAcidCompositionTable aaTable, List<Character> propertyList, List<Character> specificList, int decimalPlace){
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
		output.print(header.replace(delimiter, "_"));
		for(int i = 0; i < dList.size(); i++){
			output.print(delimiter + Utils.roundToDecimals(dList.get(i), decimalPlace));
		}
		output.println();
		output.flush();
	}

	private static void showHelp(){
		System.err.println("NAME");
		System.err.println("\tAn executable to generate physico-chemical properties of protein sequences.");
		System.err.println();
		
		System.err.println("EXAMPLES");
		System.err.println("\tjava -jar AAProperties.jar -i test.fasta -a");
		System.err.println("\t\tGenerates all possible properties.");
		System.err.println();
		System.err.println("\tjava -jar AAProperties.jar -i test.fasta -1 -3 -7");
		System.err.println("\t\tGenerates only molecular weight, extinction coefficient and isoelectric point.");
		System.err.println();
		System.err.println("\tjava -jar AAProperties.jar -i test.fasta -0 A -0 N -1");
		System.err.println("\t\tGenerates composition of two specific amino acid symbol and molecular weight.");
		System.err.println();
		
		System.err.println("OPTIONS");
		System.err.println("\tRequired");
		System.err.println("\t\t-i location of input FASTA file");
		System.err.println();

		System.err.println("\tOptional");
		System.err.println("\t\t-o location of output file [standard output (default)]");
		System.err.println("\t\t-f output format [csv (default) or tsv]");
		System.err.println("\t\t-x location of Amino Acid Composition XML file for defining amino acid composition");
		System.err.println("\t\t-y location of Element Mass XML file for defining mass of elements");
		System.err.println("\t\t-d number of decimals (int) [4 (default)]");
		System.err.println();
		
		System.err.println("\tProvide at least one of them");
		System.err.println("\t\t-a compute properties of option 1-9");
		System.err.println("\t\t-1 compute molecular weight");
		System.err.println("\t\t-2 compute absorbance");
		System.err.println("\t\t-3 compute extinction coefficient");
		System.err.println("\t\t-4 compute instability index");
		System.err.println("\t\t-5 compute apliphatic index");
		System.err.println("\t\t-6 compute average hydropathy value");
		System.err.println("\t\t-7 compute isoelectric point");
		System.err.println("\t\t-8 compute net charge at pH 7");
		System.err.println("\t\t-9 compute composition of 20 standard amino acid (A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V)");
		System.err.println("\t\t-0 compute composition of specific amino acid symbol");
		System.err.println();
	}
}
