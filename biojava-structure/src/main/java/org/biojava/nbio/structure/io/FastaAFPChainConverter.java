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
package org.biojava.nbio.structure.io;

import org.biojava.nbio.alignment.template.AlignedSequence;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.CasePreservingProteinSequenceCreator;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava.nbio.core.sequence.io.template.SequenceHeaderParserInterface;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.util.SequenceTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;

/**
 * A collection of static utilities to convert between {@link AFPChain AFPChains} and {@link FastaSequence FastaSequences}.
 * 
 * @author dmyersturnbull
 * @see StructureSequenceMatcher
 * @see FastaStructureParser
 * @see SeqRes2AtomAligner
 */
public class FastaAFPChainConverter {

	private final static Logger logger = LoggerFactory.getLogger(FastaAFPChainConverter.class);
	
		
	public static AFPChain cpFastaToAfpChain(String first, String second, Structure structure, int cpSite) throws StructureException, CompoundNotFoundException {
		ProteinSequence s1 = new ProteinSequence(first);
		s1.setUserCollection(getAlignedUserCollection(first));
		ProteinSequence s2 = new ProteinSequence(second);
		s2.setUserCollection(getAlignedUserCollection(second));
		return cpFastaToAfpChain(s1, s2, structure, cpSite);
	}

	/**
	 * Takes a structure and sequence corresponding to an alignment between a structure or sequence and itself (or even a structure with a sequence), where the result has a circular permutation site
	 * {@link cpSite} residues to the right.
	 * 
	 * @param fastaFile A FASTA file containing exactly 2 sequences, the first unpermuted and the second permuted
	 * @param cpSite
	 *            The number of residues from the beginning of the sequence at which the circular permutation site occurs; can be positive or negative; values greater than the length of the sequence
	 *            are acceptable
	 * @throws IOException 
	 * @throws StructureException 
	 */
	public static AFPChain cpFastaToAfpChain(File fastaFile, Structure structure, int cpSite) throws IOException, StructureException {
		InputStream inStream = new FileInputStream(fastaFile);
		SequenceCreatorInterface<AminoAcidCompound> creator = new CasePreservingProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet());
		SequenceHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser = new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>();
		FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(inStream, headerParser, creator);
		LinkedHashMap<String, ProteinSequence> sequences = fastaReader.process();
		inStream.close();
		Iterator<ProteinSequence> iter = sequences.values().iterator();
		ProteinSequence first = iter.next();
		ProteinSequence second = iter.next();
		return cpFastaToAfpChain(first, second, structure, cpSite);
	}

	/**
	 * Takes a structure and sequence corresponding to an alignment between a structure or sequence and itself (or even a structure with a sequence), where the result has a circular permutation site
	 * {@link cpSite} residues to the right.
	 * 
	 * @param first The unpermuted sequence
	 * @param second The sequence permuted by cpSite
	 * @param cpSite
	 *            The number of residues from the beginning of the sequence at which the circular permutation site occurs; can be positive or negative; values greater than the length of the sequence
	 *            are acceptable
	 * @throws StructureException
	 */
	public static AFPChain cpFastaToAfpChain(ProteinSequence first, ProteinSequence second, Structure structure, int cpSite)
			throws StructureException {

		if (structure == null) {
			throw new IllegalArgumentException("The structure is null");
		}

		if (first == null) {
			throw new IllegalArgumentException("The sequence is null");
		}

		// we need to find the ungapped CP site
		int gappedCpShift = 0;
		int ungappedCpShift = 0;
		while (ungappedCpShift < Math.abs(cpSite)) {
			char c;
			try {
				if (cpSite <= 0) {
					c = second.getSequenceAsString().charAt(gappedCpShift);
				} else {
					c = second.getSequenceAsString().charAt(first.getLength()-1 - gappedCpShift);
				}
			} catch (StringIndexOutOfBoundsException e) {
				throw new IllegalArgumentException("CP site of " + cpSite + " is wrong");
			}
			if (c != '-') {
				ungappedCpShift++;
			}
			gappedCpShift++;
		}

		Atom[] ca1 = StructureTools.getRepresentativeAtomArray(structure);
		Atom[] ca2 =  StructureTools.getRepresentativeAtomArray(structure); // can't use cloneCAArray because it doesn't set parent group.chain.structure
		
		ProteinSequence antipermuted = null;
		try {
			antipermuted = new ProteinSequence(SequenceTools.permuteCyclic(second.getSequenceAsString(), gappedCpShift));
		} catch (CompoundNotFoundException e) {
			// this can't happen, the original sequence comes from a ProteinSequence
			logger.error("Unexpected error while creating protein sequence: {}. This is most likely a bug.",e.getMessage() );
		}

		ResidueNumber[] residues = StructureSequenceMatcher.matchSequenceToStructure(first, structure);
		ResidueNumber[] antipermutedResidues = StructureSequenceMatcher.matchSequenceToStructure(antipermuted, structure);

		ResidueNumber[] nonpermutedResidues = new ResidueNumber[antipermutedResidues.length];
		SequenceTools.permuteCyclic(antipermutedResidues, nonpermutedResidues, -gappedCpShift);

		// nullify ResidueNumbers that have a lowercase sequence character
		if (first.getUserCollection() != null) {
			CasePreservingProteinSequenceCreator.setLowercaseToNull(first, residues);
		}
		if (second.getUserCollection() != null) {
			CasePreservingProteinSequenceCreator.setLowercaseToNull(second, nonpermutedResidues);
		}

//		for (int i = 0; i < residues.length; i++) {
//			if (residues[i] == null) {
//				System.out.print("=");
//			} else {
//				System.out.print(sequence.getSequenceAsString().charAt(i));
//			}
//		}
//		System.out.println();
//		for (int i = 0; i < residues.length; i++) {
//			if (nonpermutedResidues[i] == null) {
//				System.out.print("=");
//			} else {
//				System.out.print(second.getSequenceAsString().charAt(i));
//			}
//		}
//		System.out.println();

		return buildAlignment(ca1, ca2, residues, nonpermutedResidues);

	}

	/**
	 * Reads the file {@code fastaFile}, expecting exactly two sequences which give a pairwise alignment. Uses this and two structures to create an AFPChain corresponding to the alignment. Uses a
	 * {@link CasePreservingProteinSequenceCreator} and assumes that a residue is aligned if and only if it is given by an uppercase letter.
	 * 
	 * @see #fastaToAfpChain(ProteinSequence, ProteinSequence, Structure, Structure)
	 * @throws IOException
	 * @throws StructureException
	 */
	public static AFPChain fastaFileToAfpChain(File fastaFile, Structure structure1, Structure structure2)
			throws IOException, StructureException {
		InputStream inStream = new FileInputStream(fastaFile);
		SequenceCreatorInterface<AminoAcidCompound> creator = new CasePreservingProteinSequenceCreator(
				AminoAcidCompoundSet.getAminoAcidCompoundSet());
		SequenceHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser = new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>();
		FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
				inStream, headerParser, creator);
		LinkedHashMap<String, ProteinSequence> sequences = fastaReader.process();
		inStream.close();
		return fastaToAfpChain(sequences, structure1, structure2);
	}

	/**
	 * Returns an AFPChain corresponding to the alignment between {@code structure1} and {@code structure2}, which is given by the gapped protein sequences {@code sequence1} and {@code sequence2}. The
	 * sequences need not correspond to the entire structures, since local alignment is performed to match the sequences to structures.
	 * @throws StructureException
	 * @throws CompoundNotFoundException
	 */
	public static AFPChain fastaStringToAfpChain(String sequence1, String sequence2, Structure structure1,
			Structure structure2) throws StructureException, CompoundNotFoundException {
		ProteinSequence seq1 = new ProteinSequence(sequence1);
		ProteinSequence seq2 = new ProteinSequence(sequence2);
		return fastaToAfpChain(seq1, seq2, structure1, structure2);
	}

	/**
	 * Uses two sequences each with a corresponding structure to create an AFPChain corresponding to the alignment. Provided only for convenience since FastaReaders return such maps.
	 * 
	 * @param sequences
	 *            A Map containing exactly two entries from sequence names as Strings to gapped ProteinSequences; the name is ignored
	 * @see #fastaToAfpChain(ProteinSequence, ProteinSequence, Structure, Structure)
	 * @throws StructureException
	 */
	public static AFPChain fastaToAfpChain(Map<String, ProteinSequence> sequences, Structure structure1,
			Structure structure2) throws StructureException {

		if (sequences.size() != 2) {
			throw new IllegalArgumentException("There must be exactly 2 sequences, but there were " + sequences.size());
		}

		if (structure1 == null || structure2 == null) {
			throw new IllegalArgumentException("A structure is null");
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
	 * TODO Write comment
	 * @param sequence1
	 * @param sequence2
	 * @param structure1
	 * @param structure2
	 * @return
	 * @throws StructureException
	 * @throws CompoundNotFoundException 
	 */
	public static AFPChain fastaToAfpChain(String sequence1, String sequence2, Structure structure1,
			Structure structure2) throws StructureException, CompoundNotFoundException {
		ProteinSequence s1 = new ProteinSequence(sequence1);
		s1.setUserCollection(getAlignedUserCollection(sequence1));
		ProteinSequence s2 = new ProteinSequence(sequence2);
		s2.setUserCollection(getAlignedUserCollection(sequence2));
		return fastaToAfpChain(s1, s2, structure1, structure2);
	}

	/**
	 * Returns an AFPChain corresponding to the alignment between {@code structure1} and {@code structure2}, which is given by the gapped protein sequences {@code sequence1} and {@code sequence2}. The
	 * sequences need not correspond to the entire structures, since local alignment is performed to match the sequences to structures. Assumes that a residue is aligned if and only if it is given by
	 * an uppercase letter.
	 * @param sequence1 <em>Must</em> have {@link ProteinSequence#getUserCollection()} set to document upper- and lower-case as aligned and unaligned; see {@link #getAlignedUserCollection(String)}
	 * @throws StructureException
	 */
	public static AFPChain fastaToAfpChain(ProteinSequence sequence1, ProteinSequence sequence2, Structure structure1,
			Structure structure2) throws StructureException {

		if (structure1 == null || structure2 == null) {
			throw new IllegalArgumentException("A structure is null");
		}

		if (sequence1 == null || sequence2 == null) {
			throw new IllegalArgumentException("A sequence is null");
		}

		Atom[] ca1 = StructureTools.getRepresentativeAtomArray(structure1);
		Atom[] ca2 = StructureTools.getRepresentativeAtomArray(structure2);

		ResidueNumber[] residues1 = StructureSequenceMatcher.matchSequenceToStructure(sequence1, structure1);
		ResidueNumber[] residues2 = StructureSequenceMatcher.matchSequenceToStructure(sequence2, structure2);

		// nullify ResidueNumbers that have a lowercase sequence character
		if (sequence1.getUserCollection() != null) {
			CasePreservingProteinSequenceCreator.setLowercaseToNull(sequence1, residues1);
		}
		if (sequence2.getUserCollection() != null) {
			CasePreservingProteinSequenceCreator.setLowercaseToNull(sequence2, residues2);
		}

		return buildAlignment(ca1, ca2, residues1, residues2);

	}

	/**
	 * Provided only for convenience.
	 * 
	 * @see #fastaToAfpChain(ProteinSequence, ProteinSequence, Structure, Structure)
	 * @throws StructureException
	 */
	public static AFPChain fastaToAfpChain(SequencePair<Sequence<AminoAcidCompound>, AminoAcidCompound> alignment,
			Structure structure1, Structure structure2) throws StructureException {
		List<AlignedSequence<Sequence<AminoAcidCompound>, AminoAcidCompound>> seqs = alignment.getAlignedSequences();
		StringBuilder sb1 = new StringBuilder();
		for (AminoAcidCompound a : seqs.get(0)) {
			sb1.append(a.getBase());
		}
		try {
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
		} catch (CompoundNotFoundException e) {
			logger.error("Unexpected error while creating protein sequences: {}. This is most likely a bug.",e.getMessage());
			return null;
		}
	}

	/**
	 * Builds an {@link AFPChain} from already-matched arrays of atoms and residues.
	 * 
	 * @param ca1
	 *            An array of atoms in the first structure
	 * @param ca2
	 *            An array of atoms in the second structure
	 * @param residues1
	 *            An array of {@link ResidueNumber ResidueNumbers} in the first structure that are aligned. Only null ResidueNumbers are considered to be unaligned
	 * @param residues2
	 *            An array of {@link ResidueNumber ResidueNumbers} in the second structure that are aligned. Only null ResidueNumbers are considered to be unaligned
	 * @throws StructureException
	 */
	private static AFPChain buildAlignment(Atom[] ca1, Atom[] ca2, ResidueNumber[] residues1, ResidueNumber[] residues2)
			throws StructureException {

		// remove any gap
		// this includes the ones introduced by the nullifying above
		List<ResidueNumber> alignedResiduesList1 = new ArrayList<ResidueNumber>();
		List<ResidueNumber> alignedResiduesList2 = new ArrayList<ResidueNumber>();
		for (int i = 0; i < residues1.length; i++) {
			if (residues1[i] != null && residues2[i] != null) {
				alignedResiduesList1.add(residues1[i]);
				alignedResiduesList2.add(residues2[i]);
			}
		}

		ResidueNumber[] alignedResidues1 = alignedResiduesList1.toArray(new ResidueNumber[alignedResiduesList1.size()]);
		ResidueNumber[] alignedResidues2 = alignedResiduesList2.toArray(new ResidueNumber[alignedResiduesList2.size()]);

		AFPChain afpChain = AlignmentTools.createAFPChain(ca1, ca2, alignedResidues1, alignedResidues2);
		afpChain.setAlgorithmName("unknown");

		AlignmentTools.updateSuperposition(afpChain, ca1, ca2);

		afpChain.setBlockSize(new int[] {afpChain.getNrEQR()});
		afpChain.setBlockRmsd(new double[] {afpChain.getTotalRmsdOpt()});
		afpChain.setBlockGap(new int[] {afpChain.getGapLen()});

		return afpChain;

	}

	/**
	 * Takes a protein sequence string with capital and lowercase letters and sets its {@link ProteinSequence#getUserCollection() user collection} to record which letters are uppercase (aligned) and which are lowercase (unaligned).
	 * @param sequence Make sure <em>not</em> to use {@link ProteinSequence#getSequenceAsString()} for this, as it won't preserve upper- and lower-case
	 */
	public static List<Object> getAlignedUserCollection(String sequence) {
		List<Object> aligned = new ArrayList<Object>(sequence.length());
		for (char c : sequence.toCharArray()) {
			aligned.add(Character.isUpperCase(c));
		}
		return aligned;
	}

	/**
	 * Prints out the XML representation of an AFPChain from a file containing exactly two FASTA sequences.
	 * 
	 * @param args
	 *            A String array of fasta-file structure-1-name structure-2-name
	 * @throws StructureException
	 * @throws IOException
	 */
	public static void main(String[] args) throws StructureException, IOException {
		if (args.length != 3) {
			System.err.println("Usage: FastaAFPChainConverter fasta-file structure-1-name structure-2-name");
			return;
		}
		File fasta = new File(args[0]);
		Structure structure1 = StructureTools.getStructure(args[1]);
		Structure structure2 = StructureTools.getStructure(args[2]);
		if (structure1 == null) throw new IllegalArgumentException("No structure for " + args[1] + " was found");
		if (structure2 == null) throw new IllegalArgumentException("No structure for " + args[2] + " was found");
		AFPChain afpChain = fastaFileToAfpChain(fasta, structure1, structure2);
		String xml = AFPChainXMLConverter.toXML(afpChain);
		System.out.println(xml);
	}
	
}
