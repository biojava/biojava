package demo;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.gui.BiojavaJmol;
import org.biojava.bio.structure.io.FastaStructureParser;
import org.biojava.bio.structure.io.StructureSequenceMatcher;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.CasePreservingProteinSequenceCreator;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojava3.core.sequence.io.template.FastaHeaderParserInterface;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;

/**
 * Demo of how to use the {@link FastaStructureParser} class to read protein
 * structures from a FASTA file.
 * 
 * @author Spencer Bliven
 *
 */
public class DemoStructureFromFasta {

	public static void getAtomsFromFasta() {
		
		// Load a test sequence
		// Normally this would come from a file, eg
		// File fasta = new File("/path/to/file.fa");
		String fastaStr =
			"> 4HHB\n" +
			"vlspadktnvKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGK\n" +
			"KVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPA\n" +
			"vhasldkflasvstvltskyr\n";
		InputStream fasta;
		try {
			fasta = new ByteArrayInputStream(fastaStr.getBytes("UTF-8"));
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
			return;
		}
		
		// Create a header parser to parse the header lines into valid structure accessions.
		// The accession can be anything interpretable by AtomCache.getStructure.
		// Examples: "4HHB" (whole structure), "d4hhba_" (SCOP domain),
		//   "4HHB.A:1-15" (residue range)
		FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser;
		headerParser = new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>();
		
		// Create AtomCache to fetch structures from the PDB
		AtomCache cache = new AtomCache();
		
		// Create SequenceCreator. This converts a String to a ProteinSequence
		AminoAcidCompoundSet aaSet = AminoAcidCompoundSet.getAminoAcidCompoundSet();
		SequenceCreatorInterface<AminoAcidCompound> creator;
		creator = new CasePreservingProteinSequenceCreator(aaSet);
		
		// parse file
		FastaStructureParser parser = new FastaStructureParser(
				fasta, headerParser, creator, cache);
		try {
			parser.process();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		} catch (StructureException e) {
			e.printStackTrace();
			return;
		}

		ResidueNumber[][] residues = parser.getResidues();
		ProteinSequence[] sequences = parser.getSequences();
		Structure[] structures = parser.getStructures();
		String[] accessions = parser.getAccessions();
		
		// Optional
		// Set lowercase residues to null too
		for(int structNum = 0; structNum<sequences.length;structNum++) {
			CasePreservingProteinSequenceCreator.setLowercaseToNull(
					sequences[structNum],residues[structNum]);
		}

		// Optional, only useful for fasta files with multiple proteins
		// Remove alignment columns with a gap
		residues = StructureSequenceMatcher.removeGaps(residues);

		
		// Display the structures
		for(int i=0;i<residues.length;i++) {			
			//Display each structure
			BiojavaJmol jmol = new BiojavaJmol();
			jmol.setStructure(structures[i]);
			
			//Highlight non-null atoms
			jmol.evalString("select *; spacefill off; wireframe off; color chain; backbone 0.4;  ");
			String selectionCmd = buildJmolSelection(residues[i]);
			jmol.evalString(selectionCmd);
			jmol.evalString("backbone 1.0;");
		}
		
	}



	/**
	 * Converts an array of ResidueNumbers into a jMol selection.
	 * 
	 * <p>For example, "select 11^ :A.CA or 12^ :A.CA;" would select the
	 * CA atoms of residues 11-12 on chain A.
	 * @param residues Residues to include in the selection. Nulls are ignored.
	 * @return
	 */
	private static String buildJmolSelection(ResidueNumber[] residues) {
		StringBuilder cmd = new StringBuilder("select ");
		for(ResidueNumber res : residues) {
			if(res != null) {
				cmd.append(String.format("%d^%s:%s.CA or ", res.getSeqNum(),
						res.getInsCode()==null?" ":res.getInsCode(),
								res.getChainId()));
			}
		}
		cmd.append("none;");//easier than removing the railing 'or'
		return cmd.toString();
	}



	public static void main(String[] args) {
		getAtomsFromFasta();
	}
}
