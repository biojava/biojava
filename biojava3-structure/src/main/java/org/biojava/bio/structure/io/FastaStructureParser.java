package org.biojava.bio.structure.io;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.LinkedHashMap;

import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.template.FastaHeaderParserInterface;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;


/**
 * Reads a protein sequence from a fasta file and attempts to match it to a
 * 3D structure. Any gaps ('-') in the fasta file are preserved as null atoms in
 * the output, allowing structural alignments to be read from fasta files.
 * 
 * <p>Structures are loaded from an AtomCache. For this to work, the accession
 * for each protein should be parsed from the fasta header line into a form
 * understood by {@link AtomCache#getStructure(String)}.
 * 
 * <p>Lowercase letters are sometimes used to specify unaligned residues.
 * This information can be preserved by using a CasePreservingSequenceCreator,
 * which allows the case of residues to be accessed through the 
 * {@link ProteinSequence#getUserCollection()} method.
 * 
 * @author Spencer Bliven
 *
 */
public class FastaStructureParser {

	// inputs
	private FastaReader<ProteinSequence, AminoAcidCompound> reader;
	private AtomCache cache;
	
	// cache processed data
	private String[] accessions;
	private ProteinSequence[] sequences;
	private Structure[] structures;
	private ResidueNumber[][] residues;
	
	public FastaStructureParser(InputStream is,
			FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser,
			SequenceCreatorInterface<AminoAcidCompound> sequenceCreator,
			AtomCache cache)
	{
		this(new FastaReader<ProteinSequence, AminoAcidCompound>(
				is, headerParser, sequenceCreator),cache);
	}
	
	public FastaStructureParser(File file,
			FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser,
			SequenceCreatorInterface<AminoAcidCompound> sequenceCreator,
			AtomCache cache) throws FileNotFoundException
	{
		this(new FastaReader<ProteinSequence, AminoAcidCompound>(
				file, headerParser, sequenceCreator), cache);
	}
	
	public FastaStructureParser(FastaReader<ProteinSequence, AminoAcidCompound> reader,
			AtomCache cache) {
		this.reader = reader;
		this.cache = cache;
		this.accessions = null;
		this.sequences = null;
		this.structures = null;
		this.residues = null;
	}
	

	/**
	 * Parses the fasta file and loads it into memory.
	 * 
	 * Information can be subsequently accessed through
	 * {@link #getSequences()},
	 * {@link #getStructures()},
	 * {@link #getResidues()}, and
	 * {@link #getAccessions()}.
	 * 
	 * @throws IOException
	 * @throws StructureException
	 */
	public void process() throws IOException, StructureException {
		if(sequences == null) { // only process once, then return cached values
			LinkedHashMap<String, ProteinSequence> sequenceMap = reader.process();
			
			sequences = sequenceMap.values().toArray(new ProteinSequence[0]);
			accessions = new String[sequences.length];
			structures = new Structure[sequences.length];
			residues = new ResidueNumber[sequences.length][];

			// Match each sequence  to a series of PDB Residue numbers
			for(int i=0;i<sequences.length;i++) {
				accessions[i] = sequences[i].getAccession().getID();
								
				//System.out.println("Fetching "+accession);
				structures[i] = cache.getStructure(accessions[i]);

				residues[i] = StructureSequenceMatcher.matchSequenceToStructure(sequences[i], structures[i]);

				assert( residues[i].length == sequences[i].getLength());
			}
		}
	}
	

	/**
	 * Gets the protein sequences read from the Fasta file.
	 * Returns null if {@link #process()} has not been called.
	 * @return An array ProteinSequences from
	 *  parsing the fasta file, or null if process() hasn't been called.
	 */
	public ProteinSequence[] getSequences() {
		return sequences;
	}

	/**
	 * Gets the protein structures mapped from the Fasta file.
	 * Returns null if {@link #process()} has not been called.
	 * @return An array of Structures for each protein 
	 *  in the fasta file, or null if process() hasn't been called.
	 */
	public Structure[] getStructures() {
		return structures;
	}
	
	/**
	 * For each residue in the fasta file, return the ResidueNumber in the
	 * corresponding structure. If the residue cannot be found in the structure,
	 * that entry will be null. This can happen if that residue was not included
	 * in the PDB file (eg disordered residues), if the fasta sequence does not
	 * match the PDB sequence, or if errors occur during the matching process.
	 * @return A 2D array of ResidueNumbers, or null if process() hasn't been called.
	 * @see StructureSequenceMatcher#matchSequenceToStructure(ProteinSequence, Structure)
	 */
	public ResidueNumber[][] getResidues() {
		return residues;
	}
	
	/**
	 * Gets the protein accessions mapped from the Fasta file.
	 * Returns null if {@link #process()} has not been called.
	 * @return An array of Structures for each protein 
	 *  in the fasta file, or null if process() hasn't been called.
	 */
	public String[] getAccessions() {
		return accessions;
	}
}
