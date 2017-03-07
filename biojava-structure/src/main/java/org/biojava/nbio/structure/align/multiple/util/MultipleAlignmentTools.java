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
package org.biojava.nbio.structure.align.multiple.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.MultipleSequenceAlignment;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.phylo.DistanceMatrixCalculator;
import org.biojava.nbio.phylo.TreeConstructor;
import org.biojava.nbio.phylo.TreeConstructorType;
import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.jama.Matrix;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.phylogeny.Phylogeny;

/**
 * Utility functions for working with {@link MultipleAlignment}.
 * <p>
 * Supported functions:
 * <ul>
 * <li>Obtain the alignment as sequence strings
 * <li>Map from sequence alignment position to structure Atom
 * <li>Map from sequence alignment position to Block number
 * <li>Transform the aligned Atoms of a MultipleAlignment
 * <li>Get all the core alignment positions of the alignment
 * <li>Calculate the average residue distance of all aligned positions
 * <li>Sort Blocks in a MultipleAlignment by a specified row
 * <li>Convert a MultipleAlignment to a MultipleSequenceAlignment
 * </ul>
 *
 * @author Spencer Bliven
 * @author Aleix Lafita
 * @since 4.1.0
 *
 */
public class MultipleAlignmentTools {

	/**
	 * Calculate the sequence alignment Strings for the whole alignment. This
	 * method creates a sequence alignment where aligned residues are in
	 * uppercase and unaligned residues are in lowercase, thus providing a more
	 * compact way to represent the alignment.
	 * <p>
	 * Blocks are concatenated in the order returned by
	 * {@link MultipleAlignment#getBlocks()}, so sequences may not be
	 * sequential. Gaps are represented by '-'. Separation between different
	 * Blocks is indicated by a gap in all positions, meaning that there is a
	 * possible discontinuity.
	 *
	 * @param alignment
	 *            input MultipleAlignment
	 * @param mapSeqToStruct
	 *            provides a link from the sequence alignment position to the
	 *            structure alignment position. Specially designed for the GUI.
	 *            Has to be initialized previously and will be overwritten.
	 * @return a string for each row in the alignment, giving the 1-letter code
	 *         for each aligned residue.
	 */
	public static List<String> getSequenceAlignment(
			MultipleAlignment alignment, final List<Integer> mapSeqToStruct) {

		// Initialize sequence variables
		List<String> alnSequences = new ArrayList<String>();
		for (int str = 0; str < alignment.size(); str++)
			alnSequences.add("");
		mapSeqToStruct.clear();
		List<Atom[]> atoms = alignment.getAtomArrays();
		int globalPos = -1;

		// Initialize helper variables in constucting the sequence alignment
		List<SortedSet<Integer>> freePool = new ArrayList<SortedSet<Integer>>();
		List<SortedSet<Integer>> blockStarts = new ArrayList<SortedSet<Integer>>();
		List<List<Integer>> aligned = new ArrayList<List<Integer>>();

		// Generate freePool residues from the ones not aligned
		for (int i = 0; i < alignment.size(); i++) {
			List<Integer> residues = new ArrayList<Integer>();
			freePool.add(new TreeSet<Integer>());
			blockStarts.add(new TreeSet<Integer>());
			for (BlockSet bs : alignment.getBlockSets()) {
				for (Block b : bs.getBlocks()) {
					boolean first = true;
					for (int l = 0; l < b.length(); l++) {
						Integer residue = b.getAlignRes().get(i).get(l);
						if (residue != null) {
							if (first)
								blockStarts.get(i).add(residue);
							residues.add(residue);
							first = false;
						}
					}
				}
			}
			aligned.add(residues);
		}
		// Add any residue not aligned to the free pool for every structure
		for (int i = 0; i < alignment.size(); i++) {
			for (int k = 0; k < atoms.get(i).length; k++) {
				if (!aligned.get(i).contains(k))
					freePool.get(i).add(k);
			}
		}

		for (int b = 0; b < alignment.getBlocks().size(); b++) {
			if (b != 0) {
				// Add a gap to all structures to separate visually the Blocks
				for (int str = 0; str < alignment.size(); str++)
					alnSequences.set(str, alnSequences.get(str).concat("-"));
				mapSeqToStruct.add(-1); // means unaligned position
			}
			// Store the previous position added to the sequence alignment
			int[] previousPos = new int[alignment.size()];
			Arrays.fill(previousPos, -1);
			// Store provisional characters
			char[] provisionalChar = new char[alignment.size()];
			Arrays.fill(provisionalChar, '-');

			for (int pos = 0; pos < alignment.getBlocks().get(b).length(); pos++) {
				globalPos++;
				boolean gaps = true; // true if consecutive with the previous
				while (gaps) {
					gaps = false;
					// Loop through all the structures
					for (int str = 0; str < alignment.size(); str++) {
						// If it is the first position or before it was null
						if (previousPos[str] == -1) {
							Integer residue = alignment.getBlocks().get(b)
									.getAlignRes().get(str).get(pos);
							if (residue == null)
								provisionalChar[str] = '-';
							else {
								Atom a = atoms.get(str)[residue];
								String group = a.getGroup().getPDBName();
								provisionalChar[str] = StructureTools
										.get1LetterCode(group);
							}
						} else {
							Integer residue = alignment.getBlocks().get(b)
									.getAlignRes().get(str).get(pos);
							int nextPos = previousPos[str] + 1;
							if (residue == null) {
								if (freePool.get(str).contains(nextPos)) {
									Atom a = atoms.get(str)[nextPos];
									String g = a.getGroup().getPDBName();
									char aa = StructureTools.get1LetterCode(g);
									provisionalChar[str] = Character
											.toLowerCase(aa);
								} else
									provisionalChar[str] = '-';
							} else if (nextPos == residue) {
								Atom a = atoms.get(str)[nextPos];
								String group = a.getGroup().getPDBName();
								provisionalChar[str] = StructureTools
										.get1LetterCode(group);
							} else {
								// This means non-consecutive
								provisionalChar[str] = ' ';
								gaps = true;
							}
						}
					}// End all structure analysis

					if (gaps) {
						for (int str = 0; str < alignment.size(); str++) {
							if (provisionalChar[str] == ' ') {
								// It means this residue was non-consecutive
								Atom a = atoms.get(str)[previousPos[str] + 1];
								String group = a.getGroup().getPDBName();
								char aa = StructureTools.get1LetterCode(group);
								alnSequences
										.set(str,
												alnSequences
														.get(str)
														.concat(String.valueOf(Character
																		.toLowerCase(aa))) );
								previousPos[str]++;
							} else {
								// Insert a gap otherwise
								alnSequences.set(str, alnSequences.get(str)
										.concat("-"));
							}
						}
						mapSeqToStruct.add(-1); // unaligned position
					} else {
						// Add provisional and update the indices
						for (int str = 0; str < alignment.size(); str++) {
							alnSequences.set(
									str,
									alnSequences.get(str).concat(
											String.valueOf(provisionalChar[str])));

							if (provisionalChar[str] != '-') {
								if (alignment.getBlocks().get(b).getAlignRes()
										.get(str).get(pos) == null) {
									previousPos[str]++;
								} else {
									previousPos[str] = alignment.getBlocks()
											.get(b).getAlignRes().get(str)
											.get(pos);
								}
							}
						}
						mapSeqToStruct.add(globalPos); // alignment index
					}
				}
			} // All positions in the Block considered so far

			// Calculate the index of the next Block for every structure
			int[] blockEnds = new int[alignment.size()];
			for (int str = 0; str < alignment.size(); str++) {
				for (int res : blockStarts.get(str)) {
					if (previousPos[str] > res)
						blockEnds[str] = res;
					else {
						blockEnds[str] = res;
						break;
					}
				}
			}

			// Add the unaligned residues in between Blocks (lowercase)
			boolean allGaps = false; // true means no more residues to add
			while (!allGaps) {
				allGaps = true;
				for (int str = 0; str < alignment.size(); str++) {
					if (previousPos[str] + 1 < blockEnds[str]) {
						Atom a = atoms.get(str)[previousPos[str] + 1];
						String group = a.getGroup().getPDBName();
						char letter = StructureTools.get1LetterCode(group);

						provisionalChar[str] = Character.toLowerCase(letter);
						previousPos[str]++;
						allGaps = false;
					} else
						provisionalChar[str] = '-';
				}
				if (!allGaps) {
					for (int str = 0; str < alignment.size(); str++) {
						alnSequences.set(
								str,
								alnSequences.get(str).concat(
										String.valueOf(provisionalChar[str])) );
					}
					mapSeqToStruct.add(-1); // unaligned position
				}
			}
		}
		return alnSequences;
	}

	/**
	 * Calculate the sequence alignment Strings for the whole alignment. This
	 * method creates a sequence alignment where aligned residues are in
	 * uppercase and unaligned residues are in lowercase, thus providing a more
	 * compact way to represent the alignment.
	 * <p>
	 * Blocks are concatenated in the order returned by
	 * {@link MultipleAlignment#getBlocks()}, so sequences may not be
	 * sequential. Gaps are represented by '-'. Separation between different
	 * Blocks is indicated by a gap in all positions, meaning that there is a
	 * possible discontinuity.
	 *
	 * @param alignment
	 *            input MultipleAlignment
	 * @return String for each row in the alignment, giving the 1-letter code
	 *         for each aligned residue.
	 */
	public static List<String> getSequenceAlignment(MultipleAlignment msa) {
		return getSequenceAlignment(msa, new ArrayList<Integer>());
	}

	/**
	 * Calculate the sequence alignment Strings for the alignment Blocks in an
	 * alignment. This method creates a sequence alignment where all residues
	 * are in uppercase and a residue alone with gaps in all the other
	 * structures represents unaligned residues. Because of this latter
	 * constraint only the residues within the Blocks are represented, for a
	 * more compact alignment. For a sequence alignment of the full protein use
	 * {@link #getSequenceAlignment(MultipleAlignment)}.
	 * <p>
	 * Blocks are concatenated in the order returned by
	 * {@link MultipleAlignment#getBlocks()}, so sequences may not be
	 * sequential. Gaps between blocks are omitted, while gaps within blocks are
	 * represented by '-'. Separation between different Blocks is indicated by a
	 * gap in all positions, meaning that there is something unaligned
	 * inbetween.
	 *
	 * @param alignment
	 *            input MultipleAlignment
	 * @param mapSeqToStruct
	 *            provides a link from the sequence alignment position to the
	 *            structure alignment position. Specially designed for the GUI.
	 *            Has to be initialized previously and will be overwritten.
	 * @return a string for each row in the alignment, giving the 1-letter code
	 *         for each aligned residue.
	 */
	public static List<String> getBlockSequenceAlignment(
			MultipleAlignment alignment, List<Integer> mapSeqToStruct) {

		// Initialize sequence variables
		List<String> alnSequences = new ArrayList<String>();
		for (int str = 0; str < alignment.size(); str++)
			alnSequences.add("");
		mapSeqToStruct.clear();
		List<Atom[]> atoms = alignment.getAtomArrays();
		int globalPos = -1;

		// Loop through all the alignment Blocks in the order given
		for (int b = 0; b < alignment.getBlocks().size(); b++) {
			if (b != 0) {
				// Add a gap to all structures to separate Blocks
				for (int str = 0; str < alignment.size(); str++)
					alnSequences.set(str, alnSequences.get(str).concat("-"));
				mapSeqToStruct.add(-1); // means unaligned position
			}

			// Store the previous position added to the sequence alignment
			int[] previousPos = new int[alignment.size()];
			Arrays.fill(previousPos, -1);
			// Store provisional characters
			char[] provisionalChar = new char[alignment.size()];
			Arrays.fill(provisionalChar, '-');

			for (int pos = 0; pos < alignment.getBlocks().get(b).length(); pos++) {
				globalPos++;
				boolean gaps = true;
				while (gaps) {
					gaps = false;
					// Loop through all the structures
					for (int str = 0; str < alignment.size(); str++) {
						// If it is the first position or before it was null
						if (previousPos[str] == -1) {
							Integer residue = alignment.getBlocks().get(b)
									.getAlignRes().get(str).get(pos);
							if (residue == null)
								provisionalChar[str] = '-';
							else {
								Atom a = atoms.get(str)[residue];
								String g = a.getGroup().getPDBName();
								char aa = StructureTools.get1LetterCode(g);
								provisionalChar[str] = aa;
							}
						} else {
							Integer residue = alignment.getBlocks().get(b)
									.getAlignRes().get(str).get(pos);
							if (residue == null)
								provisionalChar[str] = '-';
							else if (previousPos[str] + 1 == residue) {
								Atom a = atoms.get(str)[residue];
								String g = a.getGroup().getPDBName();
								char aa = StructureTools.get1LetterCode(g);
								provisionalChar[str] = aa;
							} else {
								provisionalChar[str] = ' ';
								gaps = true;
							}
						}
					}// End all structures analysis

					if (gaps) {
						for (int str = 0; str < alignment.size(); str++) {
							if (provisionalChar[str] == ' ') {
								// It means this residue was non-consecutive
								for (int s2 = 0; s2 < alignment.size(); s2++) {
									if (str == s2) {
										int next = previousPos[str] + 1;
										Atom a = atoms.get(s2)[next];
										String g = a.getGroup().getPDBName();
										char aa = StructureTools
												.get1LetterCode(g);
										alnSequences.set(
												s2,
												alnSequences.get(s2).concat(
														String.valueOf(aa)) );
									} else {
										alnSequences.set(s2,
												alnSequences.get(s2)
														.concat("-"));
									}
								}
								mapSeqToStruct.add(-1); // unaligned
								previousPos[str] += 1;
							}
						}
					} else { // Append the provisional and update the indices
						for (int str = 0; str < alignment.size(); str++) {
							alnSequences.set(
									str,
									alnSequences.get(str).concat(
											String.valueOf(provisionalChar[str])) );
							if (provisionalChar[str] != '-') {
								previousPos[str] = alignment.getBlocks().get(b)
										.getAlignRes().get(str).get(pos);
							}
						}
						mapSeqToStruct.add(globalPos);
					}
				}
			}
		}
		return alnSequences;
	}

	/**
	 * Calculate the sequence alignment Strings for the alignment Blocks in an
	 * alignment. This method creates a sequence alignment where all residues
	 * are in uppercase and a residue alone with gaps in all the other
	 * structures represents unaligned residues. Because of this latter
	 * constraint only the residues within the Blocks are represented, for a
	 * more compact alignment. For a sequence alignment of the full protein use
	 * {@link #getSequenceAlignment(MultipleAlignment)}.
	 * <p>
	 * Blocks are concatenated in the order returned by
	 * {@link MultipleAlignment#getBlocks()}, so sequences may not be
	 * sequential. Gaps between blocks are omitted, while gaps within blocks are
	 * represented by '-'. Separation between different Blocks is indicated by a
	 * gap in all positions, meaning that there is something unaligned
	 * inbetween.
	 *
	 * @param alignment
	 *            input MultipleAlignment
	 * @return String for each row in the alignment, giving the 1-letter code
	 *         for each aligned residue.
	 */
	public static List<String> getBlockSequenceAlignment(MultipleAlignment ma) {
		return getBlockSequenceAlignment(ma, new ArrayList<Integer>());
	}

	/**
	 * Returns the Atom of the specified structure that is aligned in the
	 * sequence alignment position specified.
	 *
	 * @param multAln
	 *            the MultipleAlignment object from where the sequence alignment
	 *            has been generated
	 * @param mapSeqToStruct
	 *            the mapping between sequence and structure generated with the
	 *            sequence alignment
	 * @param str
	 *            the structure index of the alignment (row)
	 * @param sequencePos
	 *            the sequence alignment position (column)
	 * @return Atom the atom in that position or null if there is a gap
	 */
	public static Atom getAtomForSequencePosition(MultipleAlignment msa,
			List<Integer> mapSeqToStruct, int str, int sequencePos) {

		int seqPos = mapSeqToStruct.get(sequencePos);

		// Check if the position selected is an aligned position
		if (seqPos == -1)
			return null;
		else {
			Atom a = null;
			// Calculate the corresponding structure position
			int sum = 0;
			for (Block b : msa.getBlocks()) {
				if (sum + b.length() <= seqPos) {
					sum += b.length();
					continue;
				} else {
					for (Integer p : b.getAlignRes().get(str)) {
						if (sum == seqPos) {
							if (p != null) {
								a = msa.getAtomArrays().get(str)[p];
							}
							break;
						}
						sum++;
					}
					break;
				}
			}
			return a;
		}
	}

	/**
	 * Returns the block number of a specified position in the sequence
	 * alignment, given the mapping from structure to function.
	 *
	 * @param multAln
	 *            the MultipleAlignment object from where the sequence alignment
	 *            has been generated.
	 * @param mapSeqToStruct
	 *            the mapping between sequence and structure generated with the
	 *            sequence alignment
	 * @param sequencePos
	 *            the position in the sequence alignment (column)
	 * @return int the block index, or -1 if the position is not aligned
	 */
	public static int getBlockForSequencePosition(MultipleAlignment multAln,
			List<Integer> mapSeqToStruct, int sequencePos) {

		int seqPos = mapSeqToStruct.get(sequencePos);
		// Check if the position selected is an aligned position
		if (seqPos == -1)
			return -1;
		else {
			// Calculate the corresponding block (by iterating all Blocks)
			int sum = 0;
			int block = 0;
			for (Block b : multAln.getBlocks()) {
				if (sum + b.length() <= seqPos) {
					sum += b.length();
					block++;
					continue;
				} else
					break;
			}
			return block;
		}
	}

	/**
	 * The average residue distance Matrix contains the average distance from
	 * each residue to all other residues aligned with it.
	 * <p>
	 * Complexity: T(n,l) = O(l*n^2), if n=number of structures and l=alignment
	 * length.
	 *
	 * @param alignment
	 *            MultipleAlignment
	 * @return Matrix containing all average residue distances
	 */
	public static Matrix getAverageResidueDistances(MultipleAlignment msa) {
		List<Atom[]> transformed = transformAtoms(msa);
		return getAverageResidueDistances(transformed);
	}

	/**
	 * The average residue distance Matrix contains the average distance from
	 * each residue to all other residues aligned with it.
	 * <p>
	 * Complexity: T(n,l) = O(l*n^2), if n=number of structures and l=alignment
	 * length.
	 *
	 * @param transformed
	 *            List of Atom arrays containing only the aligned atoms of each
	 *            structure, or null if there is a gap.
	 * @return Matrix containing all average residue distances. Entry -1 means
	 *         there is a gap in the position.
	 */
	public static Matrix getAverageResidueDistances(List<Atom[]> transformed) {

		int size = transformed.size();
		int length = transformed.get(0).length;
		Matrix resDist = new Matrix(size, length, -1);

		// Calculate the average residue distances
		for (int r1 = 0; r1 < size; r1++) {
			for (int c = 0; c < transformed.get(r1).length; c++) {
				Atom refAtom = transformed.get(r1)[c];
				if (refAtom == null)
					continue;

				for (int r2 = r1 + 1; r2 < size; r2++) {
					Atom atom = transformed.get(r2)[c];
					if (atom != null) {
						double distance = Calc.getDistance(refAtom, atom);

						if (resDist.get(r1, c) == -1) {
							resDist.set(r1, c, 1 + distance);
						} else {
							resDist.set(r1, c, resDist.get(r1, c) + distance);
						}

						if (resDist.get(r2, c) == -1) {
							resDist.set(r2, c, 1 + distance);
						} else {
							resDist.set(r2, c, resDist.get(r2, c) + distance);
						}
					}
				}
			}
		}

		for (int c = 0; c < length; c++) {
			int nonNullRes = 0;
			for (int r = 0; r < size; r++) {
				if (resDist.get(r, c) != -1)
					nonNullRes++;
			}
			for (int r = 0; r < size; r++) {
				if (resDist.get(r, c) != -1) {
					resDist.set(r, c, resDist.get(r, c) / nonNullRes);
				}
			}
		}
		return resDist;
	}

	/**
	 * Transforms atoms according to the superposition stored in the alignment.
	 * <p>
	 * For each structure in the alignment, returns an atom for each
	 * representative atom in the aligned columns, omitting unaligned residues
	 * (i.e. an array of length <tt>alignment.length()</tt> ).
	 * <p>
	 * All blocks are concatenated together, so Atoms may not appear in the same
	 * order as in their parent structure. If the alignment blocks contain null
	 * residues (gaps), then the returned array will also contain null Atoms in
	 * the same positions.
	 *
	 * @param alignment
	 *            MultipleAlignment
	 * @return List of Atom arrays of only the aligned atoms of every structure
	 *         (null Atom if a gap position)
	 */
	public static List<Atom[]> transformAtoms(MultipleAlignment alignment) {
		if (alignment.getEnsemble() == null) {
			throw new NullPointerException("No ensemble set for this alignment");
		}

		List<Atom[]> atomArrays = alignment.getAtomArrays();
		List<Atom[]> transformed = new ArrayList<Atom[]>(atomArrays.size());

		// Loop through structures
		for (int i = 0; i < atomArrays.size(); i++) {

			Matrix4d transform = null;
			Atom[] curr = atomArrays.get(i); // all CA atoms from structure

			// Concatenated list of all blocks for this structure
			Atom[] transformedAtoms = new Atom[alignment.length()];
			int transformedAtomsLength = 0;

			// Each blockset gets transformed independently
			for (BlockSet bs : alignment.getBlockSets()) {

				Atom[] blocksetAtoms = new Atom[bs.length()];
				int blockPos = 0;

				for (Block blk : bs.getBlocks()) {
					if (blk.size() != atomArrays.size()) {
						throw new IllegalStateException(String.format(
								"Mismatched block size. Expected %d "
										+ "structures, found %d.",
								atomArrays.size(), blk.size()));
					}
					// Extract aligned atoms
					for (int j = 0; j < blk.length(); j++) {
						Integer alignedPos = blk.getAlignRes().get(i).get(j);
						if (alignedPos != null) {
							blocksetAtoms[blockPos] = (Atom) curr[alignedPos]
									.clone();
						}
						blockPos++;
					}
				}

				// Transform according to the blockset or alignment matrix
				Matrix4d blockTrans = null;
				if (bs.getTransformations() != null)
					blockTrans = bs.getTransformations().get(i);
				if (blockTrans == null) {
					blockTrans = transform;
				}

				for (Atom a : blocksetAtoms) {
					if (a != null)
						Calc.transform(a, blockTrans);
					transformedAtoms[transformedAtomsLength] = a;
					transformedAtomsLength++;
				}
			}
			assert (transformedAtomsLength == alignment.length());

			transformed.add(transformedAtoms);
		}
		return transformed;
	}

	/**
	 * Calculate a List of alignment indicies that correspond to the core of a
	 * Block, which means that all structures have a residue in that positon.
	 *
	 * @param block
	 *            alignment Block
	 * @return List of positions in the core of the alignment
	 */
	public static List<Integer> getCorePositions(Block block) {

		List<Integer> corePositions = new ArrayList<Integer>();

		for (int col = 0; col < block.length(); col++) {
			boolean core = true;
			for (int str = 0; str < block.size(); str++) {
				if (block.getAlignRes().get(str).get(col) == null) {
					core = false;
					break;
				}
			}
			if (core)
				corePositions.add(col);
		}
		return corePositions;
	}

	/**
	 * Sort blocks so that the specified row is in sequential order. The sort
	 * happens in place.
	 *
	 * @param blocks
	 *            List of blocks
	 * @param sortedIndex
	 *            Index of the row to be sorted
	 */
	public static void sortBlocks(List<Block> blocks, final int sortedIndex) {
		Collections.sort(blocks, new Comparator<Block>() {
			@Override
			public int compare(Block o1, Block o2) {
				// Compare the first non-null residue of each block
				List<Integer> alignres1 = o1.getAlignRes().get(sortedIndex);
				List<Integer> alignres2 = o2.getAlignRes().get(sortedIndex);
				Integer res1 = null;
				Integer res2 = null;
				for (Integer r : alignres1) {
					if (r != null) {
						res1 = r;
						break;
					}
				}
				for (Integer r : alignres2) {
					if (r != null) {
						res2 = r;
						break;
					}
				}
				return res1.compareTo(res2);
			}
		});
	}

	/**
	 * Convert a MultipleAlignment into a MultipleSequenceAlignment of AminoAcid
	 * residues. This method is only valid for protein structure alignments.
	 *
	 * @param msta
	 *            Multiple Structure Alignment
	 * @return MultipleSequenceAlignment of protein sequences
	 * @throws CompoundNotFoundException
	 */
	public static MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound> toProteinMSA(
			MultipleAlignment msta) throws CompoundNotFoundException {

		// Check that the alignment is of protein structures
		Group g = msta.getAtomArrays().get(0)[0].getGroup();
		if (!(g instanceof AminoAcid)) {
			throw new IllegalArgumentException(
					"Cannot convert to multiple sequence alignment: "
							+ "the structures aligned are not proteins");
		}

		MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound> msa = new MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound>();

		Map<String, Integer> uniqueID = new HashMap<String, Integer>();
		List<String> seqs = getSequenceAlignment(msta);
		for (int i = 0; i < msta.size(); i++) {
			// Make sure the identifiers are unique (required by AccessionID)
			String id = msta.getStructureIdentifier(i).toString();
			if (uniqueID.containsKey(id)) {
				uniqueID.put(id, uniqueID.get(id) + 1);
				id += "_" + uniqueID.get(id);
			} else
				uniqueID.put(id, 1);

			AccessionID acc = new AccessionID(id);
			ProteinSequence pseq = new ProteinSequence(seqs.get(i));
			pseq.setAccession(acc);
			msa.addAlignedSequence(pseq);
		}
		return msa;
	}
	
	public static Structure toMultimodelStructure(MultipleAlignment multAln, List<Atom[]> transformedAtoms) throws StructureException {
		PDBHeader header = new PDBHeader();
		String title = multAln.getEnsemble().getAlgorithmName() + " V."
				+ multAln.getEnsemble().getVersion() + " : ";

		for (StructureIdentifier name : multAln.getEnsemble()
				.getStructureIdentifiers()) {
			title += name.getIdentifier() + " ";
		}
		Structure artificial = getAlignedStructure(transformedAtoms);

		artificial.setPDBHeader(header);
		header.setTitle(title);
		return artificial;
	}

	/**
	 * Get an artificial Structure containing a different model for every
	 * input structure, so that the alignment result can be viewed in Jmol.
	 * The Atoms have to be rotated beforehand.
	 *
	 * @param atomArrays an array of Atoms for every aligned structure
	 * @return a structure object containing a set of models,
	 * 			one for each input array of Atoms.
	 * @throws StructureException
	 */
	public static final Structure getAlignedStructure(List<Atom[]> atomArrays)
			throws StructureException {

		Structure s = new StructureImpl();
		for (int i=0; i<atomArrays.size(); i++){
			List<Chain> model = AlignmentTools.getAlignedModel(atomArrays.get(i));
			s.addModel(model);
		}
		return s;
	}
	
	/**
	 * Calculate the RMSD matrix of a MultipleAlignment, that is, entry (i,j) of
	 * the matrix contains the RMSD between structures i and j.
	 *
	 * @param msa
	 *            Multiple Structure Alignment
	 * @return Matrix of RMSD with size the number of structures squared
	 */
	public static Matrix getRMSDMatrix(MultipleAlignment msa) {

		Matrix rmsdMat = new Matrix(msa.size(), msa.size());
		List<Atom[]> superposed = transformAtoms(msa);

		for (int i = 0; i < msa.size(); i++) {
			for (int j = i; j < msa.size(); j++) {
				if (i == j)
					rmsdMat.set(i, j, 0.0);
				List<Atom[]> compared = new ArrayList<Atom[]>();
				compared.add(superposed.get(i));
				compared.add(superposed.get(j));
				double rmsd = MultipleAlignmentScorer.getRMSD(compared);
				rmsdMat.set(i, j, rmsd);
				rmsdMat.set(j, i, rmsd);
			}
		}
		return rmsdMat;
	}

	/**
	 * Calculate a phylogenetic tree of the MultipleAlignment using Kimura
	 * distances and the Neighbor Joining algorithm from forester.
	 *
	 * @param msta
	 *            MultipleAlignment of protein structures
	 * @return Phylogeny phylogenetic tree
	 * @throws CompoundNotFoundException
	 * @throws IOException
	 */
	public static Phylogeny getKimuraTree(MultipleAlignment msta)
			throws CompoundNotFoundException, IOException {
		MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound> msa = MultipleAlignmentTools
				.toProteinMSA(msta);
		BasicSymmetricalDistanceMatrix distmat = (BasicSymmetricalDistanceMatrix) DistanceMatrixCalculator
				.kimuraDistance(msa);
		Phylogeny tree = TreeConstructor.distanceTree(distmat,
				TreeConstructorType.NJ);
		tree.setName("Kimura Tree");
		return tree;
	}

	/**
	 * Calculate a phylogenetic tree of the MultipleAlignment using
	 * dissimilarity scores (DS), based in SDM Substitution Matrix (ideal for
	 * distantly related proteins, structure-derived) and the Neighbor Joining
	 * algorithm from forester.
	 *
	 * @param msta
	 *            MultipleAlignment of protein structures
	 * @return Phylogeny phylogenetic tree
	 * @throws CompoundNotFoundException
	 * @throws IOException
	 */
	public static Phylogeny getHSDMTree(MultipleAlignment msta)
			throws CompoundNotFoundException, IOException {
		MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound> msa = MultipleAlignmentTools
				.toProteinMSA(msta);
		BasicSymmetricalDistanceMatrix distmat = (BasicSymmetricalDistanceMatrix) DistanceMatrixCalculator
				.dissimilarityScore(msa, SubstitutionMatrixHelper.getAminoAcidSubstitutionMatrix("PRLA000102"));
		Phylogeny tree = TreeConstructor.distanceTree(distmat,
				TreeConstructorType.NJ);
		tree.setName("HSDM Tree");
		return tree;
	}

	/**
	 * Calculate a phylogenetic tree of the MultipleAlignment using RMSD
	 * distances and the Neighbor Joining algorithm from forester.
	 *
	 * @param msta
	 *            MultipleAlignment of protein structures
	 * @return Phylogeny phylogenetic tree
	 * @throws CompoundNotFoundException
	 */
	public static Phylogeny getStructuralTree(MultipleAlignment msta) {
		double[][] rmsdMat = MultipleAlignmentTools.getRMSDMatrix(msta)
				.getArray();
		BasicSymmetricalDistanceMatrix rmsdDist = (BasicSymmetricalDistanceMatrix) DistanceMatrixCalculator
				.structuralDistance(rmsdMat, 1, 5, 0.4);
		// Set the identifiers of the matrix
		Map<String, Integer> alreadySeen = new HashMap<String, Integer>();
		for (int i = 0; i < msta.size(); i++) {
			// Make sure the identifiers are unique
			String id = msta.getStructureIdentifier(i).toString();
			if (alreadySeen.containsKey(id)) {
				alreadySeen.put(id, alreadySeen.get(id) + 1);
				id += "_" + alreadySeen.get(id);
			} else
				alreadySeen.put(id, 1);
			rmsdDist.setIdentifier(i, id);
		}
		Phylogeny tree = TreeConstructor.distanceTree(rmsdDist,
				TreeConstructorType.NJ);
		tree.setName("Structural Tree");
		return tree;
	}

}
