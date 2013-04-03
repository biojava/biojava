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

package org.biojava.bio.alignment;

import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.impl.SimpleGappedSequence;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;

/**
 * Smith and Waterman developed an efficient dynamic programming algorithm to
 * perform local sequence alignments, which returns the most conserved region of
 * two sequences (longest common substring with modifications). This algorithm
 * is performed by the method <code>pairwiseAlignment</code> of this class. It
 * uses affine gap penalties if and only if the expenses of a delete or insert
 * operation are unequal to the expenses of gap extension. This uses
 * significantly more memory (four times as much) and increases the runtime if
 * swapping is performed.
 * 
 * @author Andreas Dr&auml;ger
 * @author Gero Greiner
 * @author Mark Schreiber
 * @since 1.5
 */
public class SmithWaterman extends AlignmentAlgorithm {

	private short match, replace, insert, delete, gapExt;
	private SubstitutionMatrix subMatrix;

	/**
	 * Default constructor.
	 */
	public SmithWaterman() {
		this.subMatrix = null;
	}

	/**
	 * Constructs the new SmithWaterman alignment object. Alignments are only
	 * performed, if the alphabet of the given <code>SubstitutionMatrix</code>
	 * equals the alphabet of both the query and the target
	 * <code>Sequence</code>. The alignment parameters here are expenses and not
	 * scores as they are in the <code>NeedlemanWunsch</code> object. scores are
	 * just given by multiplying the expenses with <code>(-1)</code>. For
	 * example you could use parameters like "-2, 5, 3, 3, 0". If the expenses
	 * for gap extension are equal to the cost of starting a gap (delete or
	 * insert), no affine gap penalties are used, which saves memory.
	 * 
	 * @param match
	 *            expenses for a match
	 * @param replace
	 *            expenses for a replace operation
	 * @param insert
	 *            expenses for a gap opening in the query sequence
	 * @param delete
	 *            expenses for a gap opening in the target sequence
	 * @param gapExtend
	 *            expenses for the extension of a gap which was started earlier.
	 * @param matrix
	 *            the <code>SubstitutionMatrix</code> object to use.
	 */
	public SmithWaterman(short match, short replace, short insert,
			short delete, short gapExtend, SubstitutionMatrix matrix) {
		this.match = (short) -match;
		this.replace = (short) -replace;
		this.insert = (short) -insert;
		this.delete = (short) -delete;
		this.gapExt = (short) -gapExtend;
		this.subMatrix = matrix;
	}

	/**
	 * @return the delete
	 */
	public short getDelete() {
		return delete;
	}

	/**
	 * @return the gapExt
	 */
	public short getGapExt() {
		return gapExt;
	}

	/**
	 * @return the insert
	 */
	public short getInsert() {
		return insert;
	}

	/**
	 * @return the match
	 */
	public short getMatch() {
		return match;
	}

	/**
	 * @return the replace
	 */
	public short getReplace() {
		return replace;
	}

	/**
	 * @return the subMatrix
	 */
	public SubstitutionMatrix getSubMatrix() {
		return subMatrix;
	}

	/**
	 * This method computes the scores for the substitution of the i-th symbol
	 * of query by the j-th symbol of subject.
	 * 
	 * @param query
	 *            The query sequence
	 * @param subject
	 *            The target sequence
	 * @param i
	 *            The position of the symbol under consideration within the
	 *            query sequence (starting from one)
	 * @param j
	 *            The position of the symbol under consideration within the
	 *            target sequence
	 * @return The score for the given substitution.
	 */
	private short matchReplace(Sequence query, Sequence subject, int i, int j) {
		try {
			if (subMatrix != null)
				return subMatrix.getValueAt(query.symbolAt(i), subject
						.symbolAt(j));
		} catch (Exception exc) {
		}
		if (query.symbolAt(i).equals(subject.symbolAt(j)))
			return match;
		return replace;
	}

	/**
	 * This just computes the maximum of four integers.
	 * 
	 * @param w
	 * @param x
	 * @param y
	 * @param z
	 * @return the maximum of four <code>int</code>s.
	 */
	private int max(int w, int x, int y, int z) {
		if ((w > x) && (w > y) && (w > z))
			return w;
		if ((x > y) && (x > z))
			return x;
		if ((y > z))
			return y;
		return z;
	}

	/**
	 * Overrides the method inherited from the NeedlemanWunsch and performs only
	 * a local alignment. It finds only the longest common subsequence. This is
	 * good for the beginning, but it might be better to have a system to find
	 * more than only one hit within the score matrix. Therefore, one should
	 * only define the k-th best hit, where k is somehow related to the number
	 * of hits.
	 * 
	 * @see AlignmentAlgorithm#pairwiseAlignment(org.biojava.bio.symbol.SymbolList,
	 *      org.biojava.bio.symbol.SymbolList)
	 */
	public AlignmentPair pairwiseAlignment(SymbolList query, SymbolList subject)
			throws BioRuntimeException {
		int[][] scoreMatrix;
		Sequence squery = null;
		Sequence ssubject = null;

		if (query instanceof Sequence) {
			squery = (Sequence) query;
		} else {
			// make it a sequence
			squery = new SimpleSequence(query, "", "query",
					new SimpleAnnotation());
		}
		if (subject instanceof Sequence) {
			ssubject = (Sequence) subject;
		} else {
			// make it a sequence
			ssubject = new SimpleSequence(subject, "", "subject",
					new SimpleAnnotation());
		}

		if (squery.getAlphabet().equals(ssubject.getAlphabet())
				&& (subMatrix == null || squery.getAlphabet().equals(
						subMatrix.getAlphabet()))) {
			StringBuffer[] align = { new StringBuffer(), new StringBuffer() };
			long time = System.currentTimeMillis();
			int i, j, maxI = 0, maxJ = 0, queryStart = 0, subjectStart = 0;
			scoreMatrix = new int[squery.length() + 1][ssubject.length() + 1];

			/*
			 * Variables needed for traceback
			 */

			SymbolTokenization st;
			try {
				st = squery.getAlphabet().getTokenization("default");
			} catch (BioException exc) {
				throw new BioRuntimeException(exc);
			}

			/*
			 * Use affine gap panalties.
			 */
			if ((gapExt != delete) || (gapExt != insert)) {

				int[][] E = new int[squery.length() + 1][ssubject.length() + 1]; // Inserts
				int[][] F = new int[squery.length() + 1][ssubject.length() + 1]; // Deletes

				scoreMatrix[0][0] = 0;
				E[0][0] = F[0][0] = Integer.MIN_VALUE; // Double.NEGATIVE_INFINITY;
				for (i = 1; i <= squery.length(); i++) {
					scoreMatrix[i][0] = F[i][0] = 0;
					E[i][0] = Integer.MIN_VALUE; // Double.NEGATIVE_INFINITY;
				}
				for (j = 1; j <= ssubject.length(); j++) {
					scoreMatrix[0][j] = E[0][j] = 0;
					F[0][j] = Integer.MIN_VALUE; // Double.NEGATIVE_INFINITY;
				}
				for (i = 1; i <= squery.length(); i++)
					for (j = 1; j <= ssubject.length(); j++) {
						E[i][j] = Math.max(E[i][j - 1], scoreMatrix[i][j - 1]
								+ insert)
								+ gapExt;
						F[i][j] = Math.max(F[i - 1][j], scoreMatrix[i - 1][j]
								+ delete)
								+ gapExt;
						scoreMatrix[i][j] = max(0, E[i][j], F[i][j],
								scoreMatrix[i - 1][j - 1]
										+ matchReplace(squery, ssubject, i, j));

						if (scoreMatrix[i][j] > scoreMatrix[maxI][maxJ]) {
							maxI = i;
							maxJ = j;
						}
					}
				// System.out.println(printCostMatrix(G,
				// query.seqString().toCharArray(),
				// subject.seqString().toCharArray()));

				/*
				 * Here starts the traceback for affine gap penalities
				 */
				try {
					boolean[] gap_extend = { false, false };
					j = maxJ;
					for (i = maxI; i > 0;) {
						do {
							// only Deletes or Inserts or Replaces possible.
							// That's not what we want to have.
							if (scoreMatrix[i][j] == 0) {
								queryStart = i;
								subjectStart = j;
								i = j = 0;

								// Match/Replace
							} else if ((scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1]
									+ matchReplace(squery, ssubject, i, j))
									&& !(gap_extend[0] || gap_extend[1])) {

								align[0].insert(0, st.tokenizeSymbol(squery
										.symbolAt(i--)));
								align[1].insert(0, st.tokenizeSymbol(ssubject
										.symbolAt(j--)));

								// Insert || finish gap if extended gap is
								// opened
							} else if (scoreMatrix[i][j] == E[i][j]
									|| gap_extend[0]) {
								// check if gap has been extended or freshly
								// opened
								gap_extend[0] = (E[i][j] != scoreMatrix[i][j - 1]
										+ insert + gapExt);

								align[0].insert(0, '-');
								align[1].insert(0, st.tokenizeSymbol(ssubject
										.symbolAt(j--)));

								// Delete || finish gap if extended gap is
								// opened
							} else {
								// check if gap has been extended or freshly
								// opened
								gap_extend[1] = (F[i][j] != scoreMatrix[i - 1][j]
										+ delete + gapExt);

								align[0].insert(0, st.tokenizeSymbol(squery
										.symbolAt(i--)));
								align[1].insert(0, '-');
							}
						} while (j > 0);
					}
				} catch (BioException exc) {
					throw new BioRuntimeException(exc);
				}

				/*
				 * No affine gap penalties to save memory.
				 */
			} else {

				for (i = 0; i <= squery.length(); i++)
					scoreMatrix[i][0] = 0;
				for (j = 0; j <= ssubject.length(); j++)
					scoreMatrix[0][j] = 0;
				for (i = 1; i <= squery.length(); i++)
					for (j = 1; j <= ssubject.length(); j++) {

						scoreMatrix[i][j] = max(0, scoreMatrix[i - 1][j]
								+ delete, scoreMatrix[i][j - 1] + insert,
								scoreMatrix[i - 1][j - 1]
										+ matchReplace(squery, ssubject, i, j));

						if (scoreMatrix[i][j] > scoreMatrix[maxI][maxJ]) {
							maxI = i;
							maxJ = j;
						}
					}

				/*
				 * Here starts the traceback for non-affine gap penalities
				 */
				try {
					j = maxJ;
					for (i = maxI; i > 0;) {
						do {
							// only Deletes or Inserts or Replaces possible.
							// That's not what
							// we want to have.
							if (scoreMatrix[i][j] == 0) {
								queryStart = i;
								subjectStart = j;
								i = j = 0;

								// Match/Replace
							} else if (scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1]
									+ matchReplace(squery, ssubject, i, j)) {

								align[0].insert(0, st.tokenizeSymbol(squery
										.symbolAt(i--)));
								align[1].insert(0, st.tokenizeSymbol(ssubject
										.symbolAt(j--)));

								// Insert
							} else if (scoreMatrix[i][j] == scoreMatrix[i][j - 1]
									+ insert) {
								align[0].insert(0, '-');
								align[1].insert(0, st.tokenizeSymbol(ssubject
										.symbolAt(j--)));

								// Delete
							} else {
								align[0].insert(0, st.tokenizeSymbol(squery
										.symbolAt(i--)));
								align[1].insert(0, '-');
							}
						} while (j > 0);
					}
				} catch (BioException exc) {
					throw new BioRuntimeException(exc);
				}
			}

			/*
			 * From here both cases are equal again.
			 */

			try {

				/*
				 * Make equal length alignments!
				 */
				for (i = 0; i < queryStart; i++)
					align[0].insert(0, '-');
				for (i = align[0].length(); i < Math.max(maxI, maxJ); i++)
					align[0].append('-');
				for (i = 0; i < subjectStart; i++)
					align[1].insert(0, '-');
				for (i = align[1].length(); i < Math.max(maxJ, align[0]
						.length()); i++)
					align[1].append('-');
				for (i = align[0].length(); i < align[1].length(); i++)
					align[0].append('-');

				squery = new SimpleGappedSequence(
						new SimpleSequence(new SimpleSymbolList(squery
								.getAlphabet().getTokenization("token"),
								align[0].toString()), squery.getURN(), squery
								.getName(), squery.getAnnotation()));
				ssubject = new SimpleGappedSequence(
						new SimpleSequence(new SimpleSymbolList(ssubject
								.getAlphabet().getTokenization("token"),
								align[1].toString()), ssubject.getURN(),
								ssubject.getName(), ssubject.getAnnotation()));

				AlignmentPair pairalign = new AlignmentPair(squery, ssubject,
						queryStart + 1, maxI, subjectStart + 1, maxJ, subMatrix);
				pairalign.setComputationTime(System.currentTimeMillis() - time);
				pairalign.setScore(scoreMatrix[maxI][maxJ]);
				// Don't waste any memory.
				scoreMatrix = null;
				return pairalign;

			} catch (BioException exc) {
				throw new BioRuntimeException(exc);
			}
		} else
			throw new BioRuntimeException(
					"The alphabets of the sequences and the substitution matrix have to be equal.");
	}

	/**
	 * Overrides the method inherited from the NeedlemanWunsch and sets the
	 * penalty for a delete operation to the specified value. Reason: internally
	 * scores are used instead of penalties so that the value is muliplyed with
	 * -1.
	 * 
	 * @param del
	 *            costs for a single deletion operation
	 */
	public void setDelete(short del) {
		this.delete = (short) -del;
	}

	/**
	 * Overrides the method inherited from the NeedlemanWunsch and sets the
	 * penalty for an extension of any gap (insert or delete) to the specified
	 * value. Reason: internally scores are used instead of penalties so that
	 * the value is muliplyed with -1.
	 * 
	 * @param ge
	 *            costs for any gap extension
	 */
	public void setGapExt(short ge) {
		this.gapExt = (short) -ge;
	}

	/**
	 * Overrides the method inherited from the NeedlemanWunsch and sets the
	 * penalty for an insert operation to the specified value. Reason:
	 * internally scores are used instead of penalties so that the value is
	 * muliplyed with -1.
	 * 
	 * @param ins
	 *            costs for a single insert operation
	 */
	public void setInsert(short ins) {
		this.insert = (short) -ins;
	}

	/**
	 * Overrides the method inherited from the NeedlemanWunsch and sets the
	 * penalty for a match operation to the specified value. Reason: internally
	 * scores are used instead of penalties so that the value is muliplyed with
	 * -1.
	 * 
	 * @param ma
	 *            costs for a single match operation
	 */
	public void setMatch(short ma) {
		this.match = (short) -ma;
	}

	/**
	 * Overrides the method inherited from the NeedlemanWunsch and sets the
	 * penalty for a replace operation to the specified value. Reason:
	 * internally scores are used instead of penalties so that the value is
	 * muliplyed with -1.
	 * 
	 * @param rep
	 *            costs for a single replace operation
	 */
	public void setReplace(short rep) {
		this.replace = (short) -rep;
	}

	/**
	 * @param subMatrix
	 *            the subMatrix to set
	 */
	public void setSubMatrix(SubstitutionMatrix subMatrix) {
		this.subMatrix = subMatrix;
	}

}
