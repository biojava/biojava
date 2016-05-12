/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package demo;

import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.gui.BiojavaJmol;
import org.biojava.nbio.structure.io.FastaStructureParser;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.sequence.io.ProteinSequenceCreator;
import org.biojava.nbio.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava.nbio.core.sequence.io.template.SequenceHeaderParserInterface;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;

/**
 * Demo of how to use the {@link FastaStructureParser} class to read protein
 * structures from a FASTA file.
 *
 * @author Spencer Bliven
 *
 */
public class DemoStructureFromFasta {

	@SuppressWarnings("unused")
	public static void getStructureFromFasta() {

		// Load a test sequence
		// Normally this would come from a file, eg
		// File fasta = new File("/path/to/file.fa");
		String fastaStr =
			"> 4HHB\n" +
			"VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGK\n" +
			"KVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPA\n" +
			"VHASLDKFLASVSTVLTSKYR\n";
		InputStream fasta;
		try {
			fasta = new ByteArrayInputStream(fastaStr.getBytes("UTF-8"));
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
			return;
		}

		// Create a header parser to parse the header lines into valid structure accessions.
		// The resulting accession can be anything interpretable by AtomCache.getStructure.
		// Possible Examples: "4HHB" (whole structure), "d4hhba_" (SCOP domain),
		//   "4HHB.A:1-15" (residue range)
		// For this example, the built-in fasta parser will extract the correct accession.
		SequenceHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser;
		headerParser = new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>();

		// Create AtomCache to fetch structures from the PDB
		AtomCache cache = new AtomCache();

		// Create SequenceCreator. This converts a String to a ProteinSequence
		AminoAcidCompoundSet aaSet = AminoAcidCompoundSet.getAminoAcidCompoundSet();
		SequenceCreatorInterface<AminoAcidCompound> creator;
		creator = new ProteinSequenceCreator(aaSet);

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

		// Get info from the parser
		ResidueNumber[][] residues = parser.getResidues();
		ProteinSequence[] sequences = parser.getSequences();
		Structure[] structures = parser.getStructures();
		String[] accessions = parser.getAccessions();

		// Use it! For example:
		// Display the structure, highlighting the sequence
		displayStructure( structures[0], residues[0]);
	}


	/**
	 * Displays the given structure and highlights the given residues.
	 *
	 * @param structure The structure to display
	 * @param residues A list of residues to highlight
	 */
	private static void displayStructure(Structure structure,
			ResidueNumber[] residues) {
		//Display each structure
		BiojavaJmol jmol = new BiojavaJmol();
		jmol.setStructure(structure);

		//Highlight non-null atoms
		jmol.evalString("select *; spacefill off; wireframe off; color chain; backbone 0.4;  ");
		String selectionCmd = buildJmolSelection(residues);
		jmol.evalString(selectionCmd);
		jmol.evalString("backbone 1.0; select none;");
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
								res.getChainName()));
			}
		}
		cmd.append("none;");//easier than removing the railing 'or'
		return cmd.toString();
	}



	public static void main(String[] args) {
		getStructureFromFasta();
	}
}
