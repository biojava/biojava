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
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFP;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AlignmentTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.CasePreservingProteinSequenceCreator;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojava3.core.sequence.io.template.FastaHeaderParserInterface;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
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
	 * Uses a {@link CasePreservingProteinSequenceCreator} and assumes that a residue is aligned if and only if it is given by an uppercase letter.
	 * @see #fastaToAfpChain(ProteinSequence, ProteinSequence, Structure, Structure)
	 */
	public static AFPChain fastaFileToAfpChain(File fastaFile, Structure structure1, Structure structure2)
			throws Exception {
		InputStream inStream = new FileInputStream(fastaFile);
		SequenceCreatorInterface<AminoAcidCompound> creator = new CasePreservingProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet());
		FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser = new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>();
		FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(inStream, headerParser, creator);
		LinkedHashMap<String, ProteinSequence> sequences = fastaReader.process();
		inStream.close();
		return fastaToAfpChain(sequences, structure1, structure2);
	}

	private static void setUserCollection(ProteinSequence sequence) {
		if (sequence.getUserCollection() != null) return;
		List<Object> aligned = new ArrayList<Object>(sequence.getLength());
		for (char c : sequence.getSequenceAsString().toCharArray()) {
			if (Character.isUpperCase(c)) {
				aligned.add(true);
			} else {
				aligned.add(false);
			}
		}
		sequence.setUserCollection(aligned);
	}
	
	/**
	 * Returns an AFPChain corresponding to the alignment between {@code structure1} and {@code structure2}, which is given by the gapped protein sequences {@code sequence1} and {@code sequence2}. The
	 * sequences need not correspond to the entire structures, since local alignment is performed to match the sequences to structures.
	 */
	public static AFPChain fastaStringToAfpChain(String sequence1, String sequence2, Structure structure1,
			Structure structure2) throws Exception {
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

	private static Atom[] getAligned(Atom[] ca, List<ResidueNumber> residues) {
		List<Atom> alignedList = new ArrayList<Atom>();
		for (Atom atom : ca) {
			if (residues.contains(atom.getGroup().getResidueNumber())) {
				alignedList.add(atom);
			}
		}
		return alignedList.toArray(new Atom[alignedList.size()]);
	}

	/**
	 * Returns an AFPChain corresponding to the alignment between {@code structure1} and {@code structure2}, which is given by the gapped protein sequences {@code sequence1} and {@code sequence2}. The
	 * sequences need not correspond to the entire structures, since local alignment is performed to match the sequences to structures.
	 * Assumes that a residue is aligned if and only if it is given by an uppercase letter. As a side effect, sets the {@link ProteinSequence#getUserCollection()} to a list of booleans.
	 */
	public static AFPChain fastaToAfpChain(ProteinSequence sequence1, ProteinSequence sequence2, Structure structure1,
			Structure structure2) throws StructureException {

		setUserCollection(sequence1);
		setUserCollection(sequence2);
		
		if (structure1 == null || structure2 == null) {
			throw new IllegalArgumentException("A structure is null");
		}

		if (sequence1 == null || sequence2 == null) {
			throw new IllegalArgumentException("A sequence is null");
		}

		Atom[] ca1 = StructureTools.getAtomCAArray(structure1);
		Atom[] ca2 = StructureTools.getAtomCAArray(structure2);

		ResidueNumber[] residues1 = StructureSequenceMatcher.matchSequenceToStructure(sequence1, structure1);
		ResidueNumber[] residues2 = StructureSequenceMatcher.matchSequenceToStructure(sequence2, structure2);

		// nullify ResidueNumbers that have a lowercase sequence character
		if (sequence1.getUserCollection() != null) {
			CasePreservingProteinSequenceCreator.setLowercaseToNull(sequence1, residues1);
		}
		if (sequence2.getUserCollection() != null) {
			CasePreservingProteinSequenceCreator.setLowercaseToNull(sequence2, residues2);
		}

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

		if (alignedResidues1.length > 0 && alignedResidues2.length > 0) {

			// Perform singular value decomposition to find translation and rotation
			// This allows us to superimpose the aligned structures
			// This is only important for visualization of an alignment
			Atom[] alignedAtoms1 = getAligned(ca1, alignedResiduesList1);
			Atom[] alignedAtoms2 = getAligned(ca2, alignedResiduesList2);
			SVDSuperimposer svd = new SVDSuperimposer(alignedAtoms1, alignedAtoms2);
			Matrix matrix = svd.getRotation();
			Atom shift = svd.getTranslation();

			// apply the rotation and translation to a clone of the second structure
			// do NOT actually apply rotation to ca2
			Atom[] aligned2Clone = StructureTools.cloneCAArray(alignedAtoms2); // do this so we don't modify the method argument as a side effect
			for (Atom atom : aligned2Clone) {
				Calc.rotate(atom, matrix);
				Calc.shift(atom, shift);
			}

			// calculate RMSD and TM-score
			double rmsd = SVDSuperimposer.getRMS(alignedAtoms1, aligned2Clone);
			double tmScore = SVDSuperimposer.getTMScore(alignedAtoms1, aligned2Clone, ca1.length, ca2.length);

			for (AFP afp : afpChain.getAfpSet()) {
				afp.setFragLen(1);
				afp.setM(matrix);
				afp.setRmsd(rmsd);
				afp.setScore(tmScore);
			}

			afpChain.setBlockShiftVector(new Atom[] {shift});
			afpChain.setBlockRotationMatrix(new Matrix[] {matrix});
			afpChain.setBlockRmsd(new double[] {rmsd});
			afpChain.setBlockGap(new int[] {afpChain.getGapLen()});
			afpChain.setBlockSize(new int[] {alignedAtoms1.length});
			afpChain.setTotalRmsdOpt(rmsd);
			afpChain.setTMScore(tmScore);

		} else {
			afpChain.setTMScore(0);
		}

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
