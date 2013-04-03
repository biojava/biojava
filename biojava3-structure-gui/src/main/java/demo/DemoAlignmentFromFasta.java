package demo;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AlignmentTools;
import org.biojava.bio.structure.align.util.AtomCache;
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
public class DemoAlignmentFromFasta {

	public static void getAlignmentFromFasta() {
		
		// Load a test sequence
		// Normally this would come from a file, eg
		// File fasta = new File("/path/to/file.fa");
		String fastaStr =
			"> 1KQ1.A\n" +
			"mianeniqdkalenfkanqtevtvfflngFQ.MKGVIEEYDK.....YVVSLNsqgkQHLIYKh......\n" +
			".......................AISTYTVetegqastesee\n" +
			"> 1C4Q.D\n" +
			"............................tPDcVTGKVEYTKYndddtFTVKVG....DKELATnranlqs\n" +
			"lllsaqitgmtvtiktnachnggGFSEVIFr...........\n";
		
		
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
		
		// Set lowercase residues to null too
		for(int structNum = 0; structNum<sequences.length;structNum++) {
			CasePreservingProteinSequenceCreator.setLowercaseToNull(
					sequences[structNum],residues[structNum]);
		}

		// Remove alignment columns with a gap
		residues = StructureSequenceMatcher.removeGaps(residues);

		
		// Create AFPChain from the alignment
		Atom[] ca1 = StructureTools.getAtomCAArray(structures[0]);
		Atom[] ca2 = StructureTools.getAtomCAArray(structures[1]);
		AFPChain afp = AlignmentTools.createAFPChain(ca1, ca2, residues[0], residues[1]);
		
		
		try {
			StructureAlignmentDisplay.display(afp, ca1, ca2);
		} catch (StructureException e) {
			e.printStackTrace();
			return;
		}
	}


	public static void main(String[] args) {
		getAlignmentFromFasta();
	}
}
