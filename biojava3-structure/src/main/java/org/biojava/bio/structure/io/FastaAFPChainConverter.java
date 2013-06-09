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

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFP;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.jama.Matrix;
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
 * 
 */
public class FastaAFPChainConverter {

	private boolean isFastaIncomplete;

	public FastaAFPChainConverter(boolean isFastaIncomplete) {
		super();
		this.isFastaIncomplete = isFastaIncomplete;
	}

	/**
	 * Reads the file {@code fastaFile}, expecting exactly two sequences which give a pairwise alignment. Uses this and two arrays of C-alpha {@link Atom Atoms} to create an AFPChain.
	 */
	public AFPChain fastaFileToAfpChain(File fastaFile, Structure structure1, Structure structure2) throws Exception {
		LinkedHashMap<String, ProteinSequence> sequences = FastaReaderHelper.readFastaProteinSequence(fastaFile);
		return fastaToAfpChain(sequences, structure1, structure2);
	}

	public AFPChain fastaStringToAfpChain(String fastaString, Structure structure1, Structure structure2) throws Exception {
		InputStream is = null;
		LinkedHashMap<String, ProteinSequence> sequences;
		try {
			is = new ByteArrayInputStream(fastaString.getBytes()); // default encoding
			sequences = FastaReaderHelper.readFastaProteinSequence(is);
		} finally {
			if (is != null) {
				is.close();
			}
		}
		return fastaToAfpChain(sequences, structure1, structure2);
	}

	public AFPChain fastaToAfpChain(ProteinSequence sequence1, ProteinSequence sequence2, Structure structure1, Structure structure2)
			throws StructureException {

		if (sequence1 == null || sequence2 == null || structure1 == null || structure2 == null) return null;

		String seqString1 = sequence1.getSequenceAsString();
		String seqString2 = sequence2.getSequenceAsString();

		if (isFastaIncomplete) {
			structure1 = StructureSequenceMatcher.getSubstructureMatchingProteinSequence(sequence1, structure1);
			structure2 = StructureSequenceMatcher.getSubstructureMatchingProteinSequence(sequence2, structure2);
		}

		Atom[] ca1 = StructureTools.getAtomCAArray(structure1);
		Atom[] ca2 = StructureTools.getAtomCAArray(structure2);

		String ungapped1 = StructureSequenceMatcher.removeGaps(sequence1).getSequenceAsString();
		String ungapped2 = StructureSequenceMatcher.removeGaps(sequence2).getSequenceAsString();

		if (seqString1.length() != seqString2.length()) {
			throw new IllegalArgumentException("Sequence lengths differ: " + seqString1.length() + " is not equal to " + seqString2.length());
		}
		if (ungapped1.length() != ca1.length) {
			throw new IllegalArgumentException(
					"ca1 has " + ca1.length + " atoms but the sequence has " + ungapped1.length() + " residues");
		}
		if (ungapped2.length() != ca2.length) {
			throw new IllegalArgumentException(
					"ca2 has " + ca2.length + " atoms but the sequence has " + ungapped2.length() + " residues");
		}

		AFPChain afpChain = new AFPChain();
		afpChain.setCa1Length(ca1.length);
		afpChain.setCa2Length(ca2.length);

		List<Atom> participating1 = new ArrayList<Atom>();
		List<Atom> participating2 = new ArrayList<Atom>();
		
		/*
		 * Example:
		 * ACC---CGGC--TTCGCAA
		 * ACC---CGGCCCTTC--AA
		 * Gap length = 3 + 2 + 2.
		 * Notice that the first gap must be counted only once.
		 * That is, if pos1 moves, pos2 moves, or both move, we increment gapLength by 1.
		 * Whenever we move neither pos1 nor pos2, we increment alignedLength by 1.
		 */
		int pos1 = 0; // the number of gaps in the first sequence
		int pos2 = 0; // the number of gaps in the second sequence
		int alignedLength = 0;
		int gapLength = 0;
		for (int i = 0; i < seqString1.length(); i++) {
			char aa1 = seqString1.charAt(i);
			char aa2 = seqString2.charAt(i);
			boolean movedForward = false;
			if (aa1 == '-') {
				pos1++;
				movedForward = true;
			}
			if (aa2 == '-') {
				pos2++;
				movedForward = true;
			}
			if (movedForward) {
				gapLength++;
			} else {
				participating1.add(ca1[pos1]);
				participating2.add(ca2[pos2]);
				alignedLength++;
				AFP afp = new AFP();
				afp.setP1(pos1);
				afp.setP2(pos2);
				afp.setFragLen(1); // every fragment contains just one residue (a little weird, but ok)
				afp.setId(i);
				afpChain.getAfpSet().add(afp);
			}
		}

		Atom[] participating1Array = new Atom[participating1.size()];
		for (int i = 0; i < participating1.size(); i++) {
			participating1Array[i] = participating1.get(i);
		}
		Atom[] participating2Array = new Atom[participating2.size()];
		for (int i = 0; i < participating2.size(); i++) {
			participating2Array[i] = participating2.get(i);
		}

		// the best possible alignment we can get matches the two sequences together without gaps
		afpChain.setOptLength(Math.min(ca1.length, ca2.length));
		int[] optLen = new int[] { afpChain.getOptLength() };
		afpChain.setOptLen(optLen);

		afpChain.setAlnLength(alignedLength);
		afpChain.setGapLen(gapLength);

		if (participating1Array.length != participating2Array.length) {
			System.err.println(
					"ca1 has " + participating1Array.length + " atoms but ca2 has " + participating2Array.length + " atoms");
		} else {

			// Perform singular value decomposition to find translation and rotation
			// This allows us to superimpose the aligned structures
			// This is only important for visualization of an alignment
			SVDSuperimposer svd = new SVDSuperimposer(participating1Array, participating2Array);
			Matrix matrix = svd.getRotation();
			Atom shift = svd.getTranslation();

			// apply the rotation and translation to the second structure ONLY
			// From now on, ca2 is rotated
//			Atom[] ca2Clone = StructureTools.cloneCAArray(ca2); // do this so we don't modify the method argument as a side effect
			for (Atom atom : participating2Array) {
				Calc.rotate(atom, matrix);
				Calc.shift(atom, shift);
			}

			afpChain.setBlockNum(1);
			// TODO set identity, similarity, etc.
			
			// calculate RMSD and TM-score
			double rmsd = SVDSuperimposer.getRMS(participating1Array, participating2Array);
			double tmScore = SVDSuperimposer.getTMScore(participating1Array, participating2Array, ca1.length, ca2.length);

			for (AFP afp : afpChain.getAfpSet()) {
				afp.setFragLen(1);
				afp.setM(matrix);
				afp.setRmsd(rmsd);
				afp.setScore(tmScore);
			}

		}

		return afpChain;
	}

	public AFPChain fastaToAfpChain(Map<String, ProteinSequence> sequences, Structure structure1,
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

	public AFPChain fastaToAfpChain(SequencePair<Sequence<AminoAcidCompound>, AminoAcidCompound> alignment,
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
	 * Prints out an AFPChain from a file of two FASTA sequences
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
