package org.biojava.nbio.structure.align.multiple;

import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;

/**
 * This class contains functions for the conversion of {@link MultipleAlignment} to various String outputs.
 * <p>
 * Supported formats: FASTA, FatCat, Aligned Residues
 * 
 * @author Aleix Lafita
 *
 */
public class MultipleAlignmentWriter {

	/**
	 * Converts the {@link MultipleAlignment} into a multiple sequence alignment String in FASTA format.
	 * 
	 * @param alignment MultipleAlignment
	 * @return String multiple sequence alignment in FASTA format
	 */
	public static String toFASTA(MultipleAlignment alignment) {
		
		//Get the alignment sequences
		List<String> alnSequences = MultipleAlignmentTools.getSequenceAlignment(alignment);
		
		String fasta = "";
		for (int st=0; st<alignment.size(); st++){
			//Add the structure identifier as the head of the FASTA
			fasta += ">"+alignment.getEnsemble().getStructureNames().get(st)+"\n"+
					alnSequences.get(st)+"\n";
		}	
		return fasta;
	}
	
	/**
	 * Converts the {@link MultipleAlignment} into a FatCat String format. Includes summary information
	 * about the alignment in the top and a multiple sequence alignment at the bottom.
	 * 
	 * @param alignment MultipleAlignment
	 * @return String multiple sequence alignment in FASTA format
	 */
	public static String toFatCat(MultipleAlignment alignment) {
		
		List<Integer> mapSeqToStruc = new ArrayList<Integer>();
		String fatcat = alignment.toString();
		//Get the alignment sequences and the mapping
		List<String> alnSequences = MultipleAlignmentTools.getSequenceAlignment(alignment, mapSeqToStruc);
		
		
		
		
		
		
		return fatcat;
	}
	
	/**
	 * Prints the alignment in the simplest form: a list of groups of aligned residues.
	 * Format is one line per residue group, tab delimited:
	 * <ul><li>PDB number. Includes insertion code</li>
	 * <li>Chain.</li>
	 * <li>Amino Acid. Three letter code.</li>
	 * </ul>
	 * Example:
	 * <code>52	A	ALA	102	A	VAL	154	A	THR</code>
	 * <p>Note that this format loses information about blocks.
	 * 
	 * @param multAln MultipleAlignment object
	 * @return a String representation of the aligned residues.
	 */
	public static String toAlignedResidues(MultipleAlignment multAln) {
		StringWriter residueGroup = new StringWriter();

		//Write structure names & PDB codes
		for (int str=0; str<multAln.size(); str++){
			residueGroup.append("#Struct"+str+":\t");
			residueGroup.append(multAln.getEnsemble().getStructureNames().get(str));
			residueGroup.append("\n");
		}
		//Whrite header for columns
		for (int str=0; str<multAln.size(); str++) residueGroup.append("#Num"+str+"\tChain"+str+"\tAA"+str+"\t");
		residueGroup.append("\n");
		
		//Write optimally aligned pairs
		for(Block b:multAln.getBlocks()) {
			for(int res=0;res<b.length();res++) {
				for (int str=0; str<multAln.size(); str++) {
					Atom atom = multAln.getEnsemble().getAtomArrays().get(str)[b.getAlignRes().get(str).get(res)];
	
					residueGroup.append(atom.getGroup().getResidueNumber().toString());
					residueGroup.append('\t');
					residueGroup.append(atom.getGroup().getChain().getChainID());
					residueGroup.append('\t');
					residueGroup.append(atom.getGroup().getPDBName());
					residueGroup.append('\t');
				}
			residueGroup.append('\n');
			}
		}
		return residueGroup.toString();
	}
}