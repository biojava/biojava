package org.biojava.nbio.structure.align.multiple.util;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.core.util.PrettyXMLWriter;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsemble;
import org.biojava.nbio.structure.align.xml.MultipleAlignmentXMLConverter;

/**
 * This class contains functions for the conversion of 
 * {@link MultipleAlignment} to various String outputs.
 * <p>
 * Supported formats: FASTA, FatCat, Aligned Residues, 
 * Transformation Matrices, XML.
 * 
 * @author Aleix Lafita
 * @since 4.1.0
 *
 */
public class MultipleAlignmentWriter {

	/**
	 * Converts the {@link MultipleAlignment} into a multiple sequence 
	 * alignment String in FASTA format.
	 * 
	 * @param alignment MultipleAlignment
	 * @return String multiple sequence alignment in FASTA format
	 * @see MultipleAlignmentTools#getSequenceAlignment(MultipleAlignment)
	 */
	public static String toFASTA(MultipleAlignment alignment) {
		
		//Get the alignment sequences
		List<String> alnSequences = 
				MultipleAlignmentTools.getSequenceAlignment(alignment);
		
		String fasta = "";
		for (int st=0; st<alignment.size(); st++){
			//Add the structure identifier as the head of the FASTA
			fasta += ">"+alignment.getEnsemble().getStructureNames().get(st)+
					"\n"+ alnSequences.get(st)+"\n";
		}	
		return fasta;
	}
	
	/**
	 * Converts the {@link MultipleAlignment} into a FatCat String format. 
	 * Includes summary information about the alignment in the top and a 
	 * multiple sequence alignment at the bottom.
	 * 
	 * @param alignment MultipleAlignment
	 * @return String multiple sequence alignment in FASTA format
	 * @see MultipleAlignmentTools#getSequenceAlignment(MultipleAlignment)
	 */
	public static String toFatCat(MultipleAlignment alignment) {
		
		//Initialize the String and put the summary information
		StringWriter fatcat = new StringWriter();
		fatcat.append(alignment.toString()+"\n\n");
		
		//Get the alignment sequences and the mapping
		List<Integer> mapSeqToStruct = new ArrayList<Integer>();
		List<String> alnSequences = MultipleAlignmentTools.
				getSequenceAlignment(alignment, mapSeqToStruct);
		
		//Get the String of the Block Numbers for Position
		String blockNumbers = "";
		for (int pos=0; pos<alnSequences.get(0).length(); pos++){
			int blockNr = MultipleAlignmentTools.getBlockForSequencePosition(
					alignment, mapSeqToStruct, pos);
			if (blockNr != -1) {
				blockNumbers = blockNumbers.concat(""+(blockNr+1));
			} else blockNumbers = blockNumbers.concat(" ");
		}
		
		//Write the Sequence Alignment
		for (int str=0; str<alignment.size(); str++) {
			if (str<9) {
				fatcat.append("Chain 0"+(str+1)+
						": "+alnSequences.get(str)+"\n");
			} else {
				fatcat.append("Chain "+(str+1)+
						": "+alnSequences.get(str)+"\n");
			}
			if (str!=alignment.size()-1) {
				fatcat.append("          "+blockNumbers+"\n");			
			}
		}
		return fatcat.toString();
	}
	
	/**
	 * Converts the alignment to its simplest form: a list of groups of 
	 * aligned residues.
	 * Format is one line per residue group, tab delimited:
	 * <ul><li>PDB number (includes insertion code)
	 * <li>Chain
	 * <li>Amino Acid (three letter code)
	 * </li></ul>
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
			residueGroup.append("#Struct"+(str+1)+":\t");
			residueGroup.append(
					multAln.getEnsemble().getStructureNames().get(str));
			residueGroup.append("\n");
		}
		//Whrite header for columns
		for (int str=0; str<multAln.size(); str++) residueGroup.append(
				"#Num"+(str+1)+"\tChain"+(str+1)+"\tAA"+(str+1)+"\t");
		residueGroup.append("\n");
		
		//Write optimally aligned pairs
		for(Block b:multAln.getBlocks()) {
			for(int res=0;res<b.length();res++) {
				for (int str=0; str<multAln.size(); str++) {
					Integer residue = b.getAlignRes().get(str).get(res);
					if (residue == null){
						residueGroup.append("-");
						residueGroup.append('\t');
						residueGroup.append("-");
						residueGroup.append('\t');
						residueGroup.append("-");
						residueGroup.append('\t');
					} else {
						Atom atom = multAln.getAtomArrays().get(str)[residue];
		
						residueGroup.append(
								atom.getGroup().getResidueNumber().toString());
						residueGroup.append('\t');
						residueGroup.append(
								atom.getGroup().getChain().getChainID());
						residueGroup.append('\t');
						residueGroup.append(atom.getGroup().getPDBName());
						residueGroup.append('\t');
					}
				}
			residueGroup.append('\n');
			}
		}
		return residueGroup.toString();
	}
	
	/**
	 * Converts the transformation Matrices of the alignment into a String 
	 * output.
	 * 
	 * @param afpChain
	 * @return String transformation Matrices
	 */
	public static String toTransformMatrices(MultipleAlignment alignment) {

		StringBuffer txt = new StringBuffer();

		for (int bs = 0; bs < alignment.getBlockSets().size(); bs++){
			
			List<Matrix4d> btransforms = 
					alignment.getBlockSet(bs).getTransformations();
			if (btransforms == null || btransforms.size() < 1) continue;

			if (alignment.getBlockSets().size() > 1) {
				txt.append("Operations for block " );
				txt.append(bs+1);
				txt.append("\n");
			}

			for (int str=0; str<alignment.size(); str++){
				String origString = "ref";

				txt.append(String.format("     X"+(str+1)+ " = (%9.6f)*X"+ 
						origString +" + (%9.6f)*Y"+ 
						origString +" + (%9.6f)*Z"+ 
						origString +" + (%12.6f)",
						btransforms.get(str).getElement(0,0),
						btransforms.get(str).getElement(0,1),
						btransforms.get(str).getElement(0,2),
						btransforms.get(str).getElement(0,3)));
				txt.append( "\n");
				txt.append(String.format("     Y"+(str+1)+" = (%9.6f)*X"+ 
						origString +" + (%9.6f)*Y"+ 
						origString +" + (%9.6f)*Z"+ 
						origString +" + (%12.6f)",
						btransforms.get(str).getElement(1,0),
						btransforms.get(str).getElement(1,1), 
						btransforms.get(str).getElement(1,2), 
						btransforms.get(str).getElement(1,3)));
				txt.append( "\n");
				txt.append(String.format("     Z"+(str+1)+" = (%9.6f)*X"+ 
						origString +" + (%9.6f)*Y"+ 
						origString +" + (%9.6f)*Z"+ 
						origString +" + (%12.6f)",
						btransforms.get(str).getElement(2,0),
						btransforms.get(str).getElement(2,1), 
						btransforms.get(str).getElement(2,2), 
						btransforms.get(str).getElement(2,3)));
				txt.append("\n\n");
			}
		}
		return txt.toString();
	}
	
	/**
	 * Converts all the information of a multiple alignment ensemble into an
	 * XML String format. Cached variables, like transformation matrices and
	 * scores, are also converted.
	 * 
	 * @param ensemble the MultipleAlignmentEnsemble to convert.
	 * @return String XML representation of the ensemble
	 * @throws IOException
	 * @see MultipleAlignmentXMLConverter Helper methods for XML conversion
	 */
	public static String toXML(MultipleAlignmentEnsemble ensemble) 
			throws IOException {
		
		StringWriter result = new StringWriter();
		PrintWriter writer = new PrintWriter(result);
		PrettyXMLWriter xml = new PrettyXMLWriter(writer);
		
		MultipleAlignmentXMLConverter.printXMLensemble(xml, ensemble);
		
		writer.close();
		
		return result.toString();
	}
}