/**
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
 * Created on 2013-05-28
 * Created by Douglas Myers-Turnbull
 *
 * @since 3.0.6
 */
package org.biojava.bio.structure.io;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AlignmentTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.sequence.template.Sequence;

/**
 * A collection of static utilities to convert between {@link AFPChain AFPChains} and {@link FastaSequence FastaSequences}.
 * 
 * @author dmyersturnbull
 * @see StructureSequenceMatcher
 * @see FastaStructureParser
 * @see SeqRes2AtomAligner
 */
public class FastaAFPChainConverter {

	/**
	 * Reads the file {@code fastaFile}, expecting exactly two sequences which give a pairwise alignment. Uses this and two structures to create an AFPChain corresponding to the alignment.
	 * 
	 * @see #fastaToAfpChain(ProteinSequence, ProteinSequence, Structure, Structure)
	 */
	public static AFPChain fastaFileToAfpChain(File fastaFile, Structure structure1, Structure structure2)
			throws Exception {
		LinkedHashMap<String, ProteinSequence> sequences = FastaReaderHelper.readFastaProteinSequence(fastaFile);
		return fastaToAfpChain(sequences, structure1, structure2);
	}

	/**
	 * Returns an AFPChain corresponding to the alignment between {@code structure1} and {@code structure2}, which is given by the gapped protein sequences {@code sequence1} and {@code sequence2}. The
	 * sequences need not correspond to the entire structures, since local alignment is performed to match the sequences to structures.
	 */
	public static AFPChain fastaStringToAfpChain(String sequence1, String sequence2, Structure structure1,
			Structure structure2) throws Exception {
		return fastaToAfpChain(new ProteinSequence(sequence1), new ProteinSequence(sequence2), structure1, structure2);
	}

	/**
	 * Uses two sequences each with a corresponding structure to create an AFPChain corresponding to the alignment. Provided only for convenience since FastaReaders return such maps.
	 * 
	 * @param sequences
	 *            A Map containing exactly two entries from sequence names as Strings to gapped ProteinSequences; the name is ignored
	 * @see #fastaToAfpChain(ProteinSequence, ProteinSequence, Structure, Structure)
	 */
	public static AFPChain fastaToAfpChain(Map<String, ProteinSequence> sequences, Structure structure1,
			Structure structure2) throws StructureException {

		if (sequences.size() != 2) {
			throw new IllegalArgumentException("There must be exactly 2 sequences, but there were " + sequences.size());
		}

		if (structure1.getName() == null || structure2.getName() == null) {
			throw new IllegalArgumentException("A structure name is null");
		}
		if (structure1.getName().equals(structure2.getName())) {
			throw new IllegalArgumentException("");
		}

		List<ProteinSequence> seqs = new ArrayList<ProteinSequence>();
		List<String> names = new ArrayList<String>(2);
		for (Map.Entry<String, ProteinSequence> entry : sequences.entrySet()) {
			seqs.add(entry.getValue());
			names.add(entry.getKey());
		}

		return fastaToAfpChain(seqs.get(0), seqs.get(1), structure1, structure2);
	}

	/**
	 * Returns an AFPChain corresponding to the alignment between {@code structure1} and {@code structure2}, which is given by the gapped protein sequences {@code sequence1} and {@code sequence2}. The
	 * sequences need not correspond to the entire structures, since local alignment is performed to match the sequences to structures.
	 */
	public static AFPChain fastaToAfpChain(ProteinSequence sequence1, ProteinSequence sequence2, Structure structure1,
			Structure structure2) throws StructureException {

		if (sequence1 == null || sequence2 == null || structure1 == null || structure2 == null)
			return null;

		ResidueNumber[] rn1 = StructureSequenceMatcher.matchSequenceToStructure(sequence1, structure1);
		ResidueNumber[] rn2 = StructureSequenceMatcher.matchSequenceToStructure(sequence2, structure2);

		List<ResidueNumber> participating1 = new ArrayList<ResidueNumber>();
		List<ResidueNumber> participating2 = new ArrayList<ResidueNumber>();
		for (int i = 0; i < rn1.length; i++) {
			if (rn1[i] != null && rn2[i] != null) {
				participating1.add(rn1[i]);
				participating2.add(rn2[i]);
			}
		}

		ResidueNumber[] participating1Array = new ResidueNumber[participating1.size()];
		for (int i = 0; i < participating1.size(); i++) {
			participating1Array[i] = participating1.get(i);
		}
		ResidueNumber[] participating2Array = new ResidueNumber[participating2.size()];
		for (int i = 0; i < participating2.size(); i++) {
			participating2Array[i] = participating2.get(i);
		}

		Atom[] ca1 = StructureTools.getAtomCAArray(structure1);
		Atom[] ca2 = StructureTools.getAtomCAArray(structure2);

		AFPChain afpChain = AlignmentTools.createAFPChain(ca1, ca2, participating1Array, participating2Array);
		return afpChain;

	}

	/**
	 * Provided only for convenience.
	 * 
	 * @see #fastaToAfpChain(ProteinSequence, ProteinSequence, Structure, Structure)
	 */
	public static AFPChain fastaToAfpChain(SequencePair<Sequence<AminoAcidCompound>, AminoAcidCompound> alignment,
			Structure structure1, Structure structure2) throws StructureException {
		List<AlignedSequence<Sequence<AminoAcidCompound>, AminoAcidCompound>> seqs = alignment.getAlignedSequences();
		StringBuilder sb1 = new StringBuilder();
		for (AminoAcidCompound a : seqs.get(0)) {
			sb1.append(a.getBase());
		}
		ProteinSequence seq1 = new ProteinSequence(sb1.toString());
		StringBuilder sb2 = new StringBuilder();
		for (AminoAcidCompound a : seqs.get(1)) {
			sb1.append(a.getBase());
		}
		ProteinSequence seq2 = new ProteinSequence(sb2.toString());
		LinkedHashMap<String, ProteinSequence> map = new LinkedHashMap<String, ProteinSequence>();
		map.put(structure1.getName(), seq1);
		map.put(structure2.getName(), seq2);
		return fastaToAfpChain(map, structure1, structure2);
	}

	/**
	 * Prints out the XML representation of an AFPChain from a file containing exactly two FASTA sequences.
	 * 
	 * @param args
	 *            A String array of fasta-file structure-1-name structure-2-name
	 * @throws Exception
	 */
	public void main(String[] args) throws Exception {
		if (args.length != 3) {
			System.err.println("Usage: FastaAFPChainConverter fasta-file structure-1-name structure-2-name");
			return;
		}
		File fasta = new File(args[0]);
		AtomCache cache = new AtomCache();
		Structure structure1 = cache.getStructure(args[1]);
		Structure structure2 = cache.getStructure(args[2]);
		AFPChain afpChain = fastaFileToAfpChain(fasta, structure1, structure2);
		String xml = AFPChainXMLConverter.toXML(afpChain);
		System.out.println(xml);
	}

}
