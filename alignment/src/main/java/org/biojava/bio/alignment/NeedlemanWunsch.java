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

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;

import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleGappedSequence;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;

/**
 * Needleman and Wunsch defined the problem of global sequence alignments, from
 * the first till the last symbol of a sequence. This class is able to perform
 * such global sequence comparisons efficiently by dynamic programming. If
 * inserts and deletes are equally expensive and as expensive as the extension
 * of a gap, the alignment method of this class does not use affine gap
 * penalties. Otherwise it does. Those costs need four times as much memory,
 * which has significant effects on the run time, if the computer needs to swap.
 * 
 * @author Andreas Dr&auml;ger
 * @author Gero Greiner
 * @author Mark Schreiber
 * @since 1.5
 */

public class NeedlemanWunsch extends SequenceAlignment {

	/**
	 * A matrix with the size length(sequence1) times length(sequence2)
	 */
	protected int[][] CostMatrix;

	/**
	 * A matrix with the size length(alphabet) times length(alphabet)
	 */
	protected SubstitutionMatrix subMatrix;

	/**
	 * The result of a successful alignment
	 */
	protected SimpleAlignment pairalign;

	/**
	 * The result of a successful alignment as a simple String.
	 */
	protected String alignment;

	/**
	 * Expenses for inserts.
	 */
	private short insert;

	/**
	 * Expenses for deletes.
	 */
	private short delete;

	/**
	 * Expenses for the extension of a gap.
	 */
	private short gapExt;

	/**
	 * Expenses for matches.
	 */
	private short match;

	/**
	 * Expenses for replaces.
	 */
	private short replace;

	/**
	 * Constructs a new Object with the given parameters based on the
	 * Needleman-Wunsch algorithm The alphabet of sequences to be aligned will
	 * be taken from the given substitution matrix.
	 * 
	 * @param match
	 *            This gives the costs for a match operation. It is only used,
	 *            if there is no entry for a certain match of two symbols in the
	 *            substitution matrix (default value).
	 * @param replace
	 *            This is like the match parameter just the default, if there is
	 *            no entry in the substitution matrix object.
	 * @param insert
	 *            The costs of a single insert operation.
	 * @param delete
	 *            The expenses of a single delete operation.
	 * @param gapExtend
	 *            The expenses of an extension of a existing gap (that is a
	 *            previous insert or delete. If the costs for insert and delete
	 *            are equal and also equal to gapExtend, no affine gap penalties
	 *            will be used, which saves a significant amount of memory.
	 * @param subMat
	 *            The substitution matrix object which gives the costs for
	 *            matches and replaces.
	 */
	public NeedlemanWunsch(short match, short replace, short insert,
			short delete, short gapExtend, SubstitutionMatrix subMat) {
		this.subMatrix = subMat;
		this.insert = insert;
		this.delete = delete;
		this.gapExt = gapExtend;
		this.match = match;
		this.replace = replace;
		this.alignment = "";
	}

	/**
	 * Sets the substitution matrix to be used to the specified one. Afterwards
	 * it is only possible to align sequences of the alphabet of this
	 * substitution matrix.
	 * 
	 * @param matrix
	 *            an instance of a substitution matrix.
	 */
	public void setSubstitutionMatrix(SubstitutionMatrix matrix) {
		this.subMatrix = matrix;
	}

	/**
	 * Sets the penalty for an insert operation to the specified value.
	 * 
	 * @param ins
	 *            costs for a single insert operation
	 */
	public void setInsert(short ins) {
		this.insert = ins;
	}

	/**
	 * Sets the penalty for a delete operation to the specified value.
	 * 
	 * @param del
	 *            costs for a single deletion operation
	 */
	public void setDelete(short del) {
		this.delete = del;
	}

	/**
	 * Sets the penalty for an extension of any gap (insert or delete) to the
	 * specified value.
	 * 
	 * @param ge
	 *            costs for any gap extension
	 */
	public void setGapExt(short ge) {
		this.gapExt = ge;
	}

	/**
	 * Sets the penalty for a match operation to the specified value.
	 * 
	 * @param ma
	 *            costs for a single match operation
	 */
	public void setMatch(short ma) {
		this.match = ma;
	}

	/**
	 * Sets the penalty for a replace operation to the specified value.
	 * 
	 * @param rep
	 *            costs for a single replace operation
	 */
	public void setReplace(short rep) {
		this.replace = rep;
	}

	/**
	 * Returns the current expenses of a single insert operation.
	 * 
	 * @return insert
	 */
	public short getInsert() {
		return insert;
	}

	/**
	 * Returns the current expenses of a single delete operation.
	 * 
	 * @return delete
	 */
	public short getDelete() {
		return delete;
	}

	/**
	 * Returns the current expenses of any extension of a gap operation.
	 * 
	 * @return gapExt
	 */
	public short getGapExt() {
		return gapExt;
	}

	/**
	 * Returns the current expenses of a single match operation.
	 * 
	 * @return match
	 */
	public short getMatch() {
		return match;
	}

	/**
	 * Returns the current expenses of a single replace operation.
	 * 
	 * @return replace
	 */
	public short getReplace() {
		return replace;
	}

	/**
	 * Prints a String representation of the CostMatrix for the given Alignment
	 * on the screen. This can be used to get a better understanding of the
	 * algorithm. There is no other purpose. This method also works for all
	 * extensions of this class with all kinds of matrices.
	 * 
	 * @param CostMatrix
	 *            The matrix that contains all expenses for swapping symbols.
	 * @param queryChar
	 *            a character representation of the query sequence (
	 *            <code>mySequence.seqString().toCharArray()</code>).
	 * @param targetChar
	 *            a character representation of the target sequence.
	 * @return a String representation of the matrix.
	 */
	public static String printCostMatrix(int[][] CostMatrix, char[] queryChar,
			char[] targetChar) {
		int line, col;
		StringBuilder output = new StringBuilder('\t');
		String ls = System.getProperty("line.separator");

		for (col = 0; col <= targetChar.length; col++)
			if (col == 0) {
				output.append('['); 
				output.append(col);
				output.append("]\t");
			} else {
				output.append('[');
				output.append(targetChar[col - 1]);
				output.append("]\t");
			}
		for (line = 0; line <= queryChar.length; line++) {
			if (line == 0) {
				output.append(ls);
				output.append('[');
				output.append(line);
				output.append("]\t");
			} else {
				output.append(ls);
				output.append('[');
				output.append(queryChar[line - 1]);
				output.append("]\t");
			}
			for (col = 0; col <= targetChar.length; col++) {
				output.append(CostMatrix[line][col]);
				output.append('\t');
			}
		}
		output.append(ls);
		output.append("delta[Edit] = ");
		output.append(CostMatrix[line - 1][col - 1]);
		output.append(ls);
		return output.toString();
	}

	/**
	 * prints the alignment String on the screen (standard output).
	 * 
	 * @param align
	 *            The parameter is typically given by the
	 *            {@link #getAlignmentString() getAlignmentString()} method.
	 */
	public static void printAlignment(String align) {
		System.out.print(align);
	}

	/**
	 * This method is good if one wants to reuse the alignment calculated by
	 * this class in another BioJava class. It just performs
	 * {@link #pairwiseAlignment(SymbolList, SymbolList) pairwiseAlignment} and
	 * returns an <code>Alignment</code> instance containing the two aligned
	 * sequences.
	 * 
	 * @return Alignment object containing the two gapped sequences constructed
	 *         from query and target.
	 * @throws Exception
	 */
	public Alignment getAlignment(SymbolList query, SymbolList target)
			throws Exception {
		pairwiseAlignment(query, target);
		return pairalign;
	}

	/**
	 * This gives the edit distance according to the given parameters of this
	 * certain object. It returns just the last element of the internal cost
	 * matrix (left side down). So if you extend this class, you can just do the
	 * following:
	 * <code>int myDistanceValue = foo; this.CostMatrix = new int[1][1]; this.CostMatrix[0][0] = myDistanceValue;</code>
	 * 
	 * @return returns the edit_distance computed with the given parameters.
	 */
	public int getEditDistance() {
		return CostMatrix[CostMatrix.length - 1][CostMatrix[CostMatrix.length - 1].length - 1];
	}

	/**
	 * This just computes the minimum of three integer values.
	 * 
	 * @param x
	 * @param y
	 * @param z
	 * @return Gives the minimum of three integers
	 */
	protected static int min(int x, int y, int z) {
		if ((x < y) && (x < z))
			return x;
		if (y < z)
			return y;
		return z;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see toolbox.align.SequenceAlignment#getAlignment()
	 */
	public String getAlignmentString() throws BioException {
		return alignment;
	}

	public Alignment getAlignment() {
		return pairalign;
	}
	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * toolbox.align.SequenceAlignment#alignAll(org.biojava.bio.seq.SequenceIterator
	 * , org.biojava.bio.seq.db.SequenceDB)
	 */
	public List<Alignment> alignAll(SequenceIterator source,
			SequenceDB subjectDB) throws NoSuchElementException, BioException {
		List<Alignment> l = new LinkedList<Alignment>();
		while (source.hasNext()) {
			Sequence query = source.nextSequence();
			// compare all the sequences of both sets.
			SequenceIterator target = subjectDB.sequenceIterator();
			while (target.hasNext())
				try {
					l.add(getAlignment(query, target.nextSequence()));
					// pairwiseAlignment(query, target.nextSequence());
				} catch (Exception exc) {
					exc.printStackTrace();
				}
		}
		return l;
	}

	/**
	 * Global pairwise sequence alignment of two BioJava-Sequence objects
	 * according to the Needleman-Wunsch-algorithm.
	 * 
	 * @see org.biojava.bio.alignment.SequenceAlignment#pairwiseAlignment(org.biojava.bio.symbol.SymbolList,
	 *      org.biojava.bio.symbol.SymbolList)
	 */
	public int pairwiseAlignment(SymbolList query, SymbolList subject)
			throws BioException {
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

		SymbolTokenization st = null;
		st = subMatrix.getAlphabet().getTokenization("default");

		if (squery.getAlphabet().equals(ssubject.getAlphabet())
				&& squery.getAlphabet().equals(subMatrix.getAlphabet())) {

			StringBuffer[] align = { new StringBuffer(), new StringBuffer() };

			long time = System.currentTimeMillis();
			int i, j;
			this.CostMatrix = new int[squery.length() + 1][ssubject.length() + 1];

			/*
			 * Variables for the traceback
			 */

			StringBuffer path = new StringBuffer();

			// construct the matrix:
			CostMatrix[0][0] = 0;

			/*
			 * If we want to have affine gap penalties, we have to initialise
			 * additional matrices: If this is not necessary, we won't do that
			 * (because it's expensive).
			 */
			if ((gapExt != delete) || (gapExt != insert)) {

				int[][] E = new int[squery.length() + 1][ssubject.length() + 1]; // Inserts
				int[][] F = new int[squery.length() + 1][ssubject.length() + 1]; // Deletes

				E[0][0] = F[0][0] = Integer.MAX_VALUE; // Double.MAX_VALUE;
				for (i = 1; i <= squery.length(); i++) {
					// CostMatrix[i][0] = CostMatrix[i-1][0] + delete;
					E[i][0] = Integer.MAX_VALUE; // Double.POSITIVE_INFINITY;
					CostMatrix[i][0] = F[i][0] = delete + i * gapExt;
				}
				for (j = 1; j <= ssubject.length(); j++) {
					// CostMatrix[0][j] = CostMatrix[0][j - 1] + insert;
					F[0][j] = Integer.MAX_VALUE; // Double.POSITIVE_INFINITY;
					CostMatrix[0][j] = E[0][j] = insert + j * gapExt;
				}
				for (i = 1; i <= squery.length(); i++)
					for (j = 1; j <= ssubject.length(); j++) {
						E[i][j] = Math.min(E[i][j - 1], CostMatrix[i][j - 1]
								+ insert)
								+ gapExt;
						F[i][j] = Math.min(F[i - 1][j], CostMatrix[i - 1][j]
								+ delete)
								+ gapExt;
						CostMatrix[i][j] = min(E[i][j], F[i][j],
								CostMatrix[i - 1][j - 1]
										- matchReplace(squery, ssubject, i, j));
					}

				/*
				 * Traceback for affine gap penalties.
				 */

				boolean[] gap_extend = { false, false };
				j = this.CostMatrix[CostMatrix.length - 1].length - 1;

				for (i = this.CostMatrix.length - 1; i > 0;) {
					do {
						// only Insert.
						if (i == 0) {
							align[0].insert(0, '~');
							align[1].insert(0, st.tokenizeSymbol(ssubject
									.symbolAt(j--)));
							path.insert(0, ' ');

							// only Delete.
						} else if (j == 0) {
							align[0].insert(0, st.tokenizeSymbol(squery
									.symbolAt(i--)));
							align[1].insert(0, '~');
							path.insert(0, ' ');

							// Match/Replace
						} else if ((CostMatrix[i][j] == CostMatrix[i - 1][j - 1]
								- matchReplace(squery, ssubject, i, j))
								&& !(gap_extend[0] || gap_extend[1])) {
							if (squery.symbolAt(i) == ssubject.symbolAt(j))
								path.insert(0, '|');
							else
								path.insert(0, ' ');
							align[0].insert(0, st.tokenizeSymbol(squery
									.symbolAt(i--)));
							align[1].insert(0, st.tokenizeSymbol(ssubject
									.symbolAt(j--)));

							// Insert || finish gap if extended gap is
							// opened
						} else if (CostMatrix[i][j] == E[i][j] || gap_extend[0]) {
							// check if gap has been extended or freshly
							// opened
							gap_extend[0] = (E[i][j] != CostMatrix[i][j - 1]
									+ insert + gapExt);

							align[0].insert(0, '-');
							align[1].insert(0, st.tokenizeSymbol(ssubject
									.symbolAt(j--)));
							path.insert(0, ' ');

							// Delete || finish gap if extended gap is
							// opened
						} else {
							// check if gap has been extended or freshly
							// opened
							gap_extend[1] = (F[i][j] != CostMatrix[i - 1][j]
									+ delete + gapExt);

							align[0].insert(0, st.tokenizeSymbol(squery
									.symbolAt(i--)));
							align[1].insert(0, '-');
							path.insert(0, ' ');
						}
					} while (j > 0);
				}

				/*
				 * No affine gap penalties, constant gap penalties, which is
				 * much faster and needs less memory.
				 */

			} else {

				for (i = 1; i <= squery.length(); i++)
					CostMatrix[i][0] = CostMatrix[i - 1][0] + delete;
				for (j = 1; j <= ssubject.length(); j++)
					CostMatrix[0][j] = CostMatrix[0][j - 1] + insert;
				for (i = 1; i <= squery.length(); i++)
					for (j = 1; j <= ssubject.length(); j++) {
						CostMatrix[i][j] = min(CostMatrix[i - 1][j] + delete,
								CostMatrix[i][j - 1] + insert,
								CostMatrix[i - 1][j - 1]
										- matchReplace(squery, ssubject, i, j));
					}

				/*
				 * Traceback for constant gap penalties.
				 */
				j = this.CostMatrix[CostMatrix.length - 1].length - 1;

				// System.out.println(printCostMatrix(CostMatrix,
				// query.seqString().toCharArray(),
				// subject.seqString().toCharArray()));

				for (i = this.CostMatrix.length - 1; i > 0;) {
					do {
						// only Insert.
						if (i == 0) {
							align[0].insert(0, '~');
							align[1].insert(0, st.tokenizeSymbol(ssubject
									.symbolAt(j--)));
							path.insert(0, ' ');

							// only Delete.
						} else if (j == 0) {
							align[0].insert(0, st.tokenizeSymbol(squery
									.symbolAt(i--)));
							align[1].insert(0, '~');
							path.insert(0, ' ');

							// Match/Replace
						} else if (CostMatrix[i][j] == CostMatrix[i - 1][j - 1]
								- matchReplace(squery, ssubject, i, j)) {

							if (squery.symbolAt(i) == ssubject.symbolAt(j))
								path.insert(0, '|');
							else
								path.insert(0, ' ');
							align[0].insert(0, st.tokenizeSymbol(squery
									.symbolAt(i--)));
							align[1].insert(0, st.tokenizeSymbol(ssubject
									.symbolAt(j--)));

							// Insert
						} else if (CostMatrix[i][j] == CostMatrix[i][j - 1]
								+ insert) {
							align[0].insert(0, '-');
							align[1].insert(0, st.tokenizeSymbol(ssubject
									.symbolAt(j--)));
							path.insert(0, ' ');

							// Delete
						} else {
							align[0].insert(0, st.tokenizeSymbol(squery
									.symbolAt(i--)));
							align[1].insert(0, '-');
							path.insert(0, ' ');
						}
					} while (j > 0);
				}

			}

			/*
			 * From here both cases are equal again.
			 */
			squery = new SimpleGappedSequence(new SimpleSequence(
					new SimpleSymbolList(squery.getAlphabet().getTokenization(
							"token"), align[0].toString()), squery.getURN(),
					squery.getName(), squery.getAnnotation()));
			ssubject = new SimpleGappedSequence(new SimpleSequence(
					new SimpleSymbolList(ssubject.getAlphabet()
							.getTokenization("token"), align[1].toString()),
					ssubject.getURN(), ssubject.getName(), ssubject
							.getAnnotation()));
			Map<String, SymbolList> m = new HashMap<String, SymbolList>();
			m.put(squery.getName(), squery);
			m.put(ssubject.getName(), ssubject);
			pairalign = new SimpleAlignment(m);

			// this.printCostMatrix(queryChar, targetChar); // only for
			// tests
			// important
			this.alignment = formatOutput(squery.getName(), // name of the
					// query
					// sequence
					ssubject.getName(), // name of the target sequence
					align, // the String representation of the alignment
					path, // String match/missmatch representation
					0, // Start position of the alignment in the query
					// sequence
					CostMatrix.length - 1, // End position of the alignment
					// in the query sequence
					CostMatrix.length - 1, // length of the query sequence
					0, // Start position of the alignment in the target
					// sequence
					CostMatrix[0].length - 1, // End position of the
					// alignment in the target sequence
					CostMatrix[0].length - 1, // length of the target
					// sequence
					getEditDistance(), // the edit distance
					System.currentTimeMillis() - time, subMatrix, st)
					+ System.getProperty("line.separator"); // time
			// consumption

			// System.out.println(printCostMatrix(CostMatrix,
			// query.seqString().toCharArray(),
			// subject.seqString().toCharArray()));
			int score = getEditDistance();
			pairalign.setScore(score);
			return score;
		} else
			throw new BioException(
					"Alphabet missmatch occured: sequences with different alphabet cannot be aligned.");
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
	private int matchReplace(Sequence query, Sequence subject, int i, int j) {
		try {
			return subMatrix.getValueAt(query.symbolAt(i), subject.symbolAt(j));
		} catch (Exception exc) {
			if (query.symbolAt(i).getMatches().contains(subject.symbolAt(j))
					|| subject.symbolAt(j).getMatches().contains(
							query.symbolAt(i)))
				return -match;
			return -replace;
		}
	}

}
