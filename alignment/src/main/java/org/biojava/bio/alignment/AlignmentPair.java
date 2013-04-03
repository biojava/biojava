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

import java.util.Formatter;
import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * This class stores the result of an alignment procedure that creates a
 * pairwise alignment of two sequences. Currently, these sequences must have the
 * identical length (this may be changed in the future). A format routine
 * produces a BLAST-like output for the sequences but all necessary information
 * to visualize the alignment are contained in this class.
 * 
 * @author Andreas Dr&auml;ger
 * 
 */
public class AlignmentPair extends SimpleAlignment {

	/**
	 * Generated serial version identifier
	 */
	private static final long serialVersionUID = -8834131912021612261L;

	/**
	 * 
	 * @param query
	 * @param subject
	 * @param algorithm
	 * @return
	 * @throws Exception
	 */
	public static AlignmentPair align(Sequence query, Sequence subject,
			AlignmentAlgorithm algorithm) throws Exception {
		AlignmentPair pair = algorithm.pairwiseAlignment(query, subject);
		return pair;
	}

	/**
	 * Creates the Map required by the super class.
	 * 
	 * @param s1
	 * @param s2
	 * @return
	 */
	private static Map<String, SymbolList> createHashMap(Sequence s1,
			Sequence s2) {
		Map<String, SymbolList> m = new HashMap<String, SymbolList>();
		m.put(s1.getName(), s1);
		m.put(s2.getName(), s2);
		return m;
	}

	/**
	 * Number of identical elements in both sequences
	 */
	private int identicals;

	/**
	 * Start position in the query
	 */
	private int queryStart;

	/**
	 * @return the queryStart
	 */
	public int getQueryStart() {
		return queryStart;
	}

	/**
	 * @param queryStart
	 *            the queryStart to set
	 * @throws BioException
	 */
	void setQueryStart(int queryStart) throws BioException {
		this.queryStart = queryStart;
		init();
	}

	/**
	 * End position in the query
	 */
	private int queryEnd;

	/**
	 * @return the queryEnd
	 */
	public int getQueryEnd() {
		return queryEnd;
	}

	/**
	 * @param queryEnd
	 *            the queryEnd to set
	 * @throws BioException
	 */
	void setQueryEnd(int queryEnd) throws BioException {
		this.queryEnd = queryEnd;
		init();
	}

	/**
	 * Number of gaps in query
	 */
	private int nGapsQ;

	/**
	 * Number of gaps in subject
	 */
	private int nGapsS;

	/**
	 * Length of the query sequence
	 */
	private Sequence query;

	/**
	 * Start position in the subject
	 */
	private int subjectStart;

	/**
	 * @return the subjectStart
	 */
	public int getSubjectStart() {
		return subjectStart;
	}

	/**
	 * @param subjectStart
	 *            the subjectStart to set
	 * @throws BioException
	 */
	void setSubjectStart(int subjectStart) throws BioException {
		this.subjectStart = subjectStart;
		init();
	}

	/**
	 * End position in the subject
	 */
	private int subjectEnd;

	/**
	 * @return the subjectEnd
	 */
	public int getSubjectEnd() {
		return subjectEnd;
	}

	/**
	 * @param subjectEnd
	 *            the subjectEnd to set
	 * @throws BioException
	 */
	void setSubjectEnd(int subjectEnd) throws BioException {
		this.subjectEnd = subjectEnd;
		init();
	}

	/**
	 * Number of similar symbols
	 */
	private int similars;

	/**
	 * The subject sequence.
	 */
	private Sequence subject;

	/**
	 * Reference to the underlying substitution matrix of this alignment.
	 */
	private SubstitutionMatrix subMatrix;

	/**
	 * Time consumption to create this alignment.
	 */
	private long time;

	/**
	 * 
	 * @param queryStart
	 *            the start position in the query, where the alignment starts.
	 *            For example zero for normal Needleman-Wunsch-Alignments.
	 * @param queryEnd
	 *            the end position, that means the sequence coordinate, which is
	 *            the last symbol of the query sequence. Counting starts at
	 *            zero!
	 * @param queryLength
	 *            The length of the query sequence without gaps.
	 * @param subjectStart
	 *            These are all the same for the target. Have a look at these
	 *            above.
	 * @param subjectEnd
	 * @param subMatrix
	 *            the subsitution Matrix used for calculating the alignment
	 * @throws IllegalArgumentException
	 * @throws BioException
	 */
	public AlignmentPair(Sequence query, Sequence subject, int queryStart,
			int queryEnd, int subjectStart, int subjectEnd,
			SubstitutionMatrix subMatrix) throws IllegalArgumentException,
			BioException {
		super(createHashMap(query, subject));
		this.subMatrix = subMatrix;
		this.query = query;
		this.subject = subject;
		this.queryStart = queryStart;
		this.queryEnd = queryEnd;
		this.subjectStart = subjectStart;
		this.subjectEnd = subjectEnd;
		init();
	}

	/**
	 * 
	 * @param query
	 * @param subject
	 * @param subMatrix
	 * @throws IllegalArgumentException
	 * @throws BioException
	 */
	public AlignmentPair(Sequence query, Sequence subject,
			SubstitutionMatrix subMatrix) throws IllegalArgumentException,
			BioException {
		this(query, subject, 1, query.length(), 1, subject.length(), subMatrix);
	}

	/**
	 * 
	 * @throws BioException
	 */
	private void init() throws BioException {
		similars = 0;
		identicals = 0;
		nGapsQ = 0;
		nGapsS = 0;
		for (int i = 0; i < Math.min(queryEnd - queryStart , subjectEnd
				- subjectStart)  ; i++) {
		   
			Symbol a = query.symbolAt(i + queryStart);
			Symbol b = subject.symbolAt(i + subjectStart);
			boolean gap = false;
			if (a.equals(b))
				identicals++;

			// get score for this pair. if it is positive, they are similar...
			if (a.equals(query.getAlphabet().getGapSymbol())) {
				nGapsQ++;
				gap = true;
			}
			if (b.equals(subject.getAlphabet().getGapSymbol())) {
				nGapsS++;
				gap = true;
			}
			if (!gap && subMatrix != null && subMatrix.getValueAt(a, b) > 0)
				similars++;
		}
	}

	/**
	 * @return the time
	 */
	public long getComputationTime() {
		return time;
	}

	/**
	 * @return the ngapsq
	 */
	public int getNumGapsInQuery() {
		return nGapsQ;
	}

	/**
	 * @return the ngapst
	 */
	public int getNumGapsInSubject() {
		return nGapsS;
	}

	/**
	 * @return the identicals
	 */
	public int getNumIdenticals() {
		return identicals;
	}

	/**
	 * @return the similars
	 */
	public int getNumSimilars() {
		return similars;
	}

	/**
	 * 
	 * @return
	 */
	public float getPercentIdentityQuery() {
		return identicals / (float) (query.length() - nGapsQ) * 100;
	}

	/**
	 * 
	 * @return
	 */
	public float getPercentIdentitySubject() {
		return identicals / (float) (subject.length() - nGapsS) * 100;
	}

	/**
	 * 
	 * @return
	 */
	public float getPercentSimilarityQuery() {
		return similars / (float) query.length() * 100;
	}

	/**
	 * 
	 * @return
	 */
	public float getPercentSimilaritySubject() {
		return similars / (float) subject.length() * 100;
	}

	/**
	 * 
	 * @return
	 */
	public int getQueryLength() {
		return query.length();
	}

	/**
	 * 
	 * @return
	 */
	public int getSubjectLength() {
		return subject.length();
	}

	/**
	 * @return the subMatrix
	 */
	public SubstitutionMatrix getSubstitutionMatrix() {
		return subMatrix;
	}

	/**
	 * @param time
	 *            The time in milliseconds, which was needed to generate the
	 *            alignment.
	 */
	void setComputationTime(long time) {
		this.time = time;
	}

	/**
	 * 
	 * @return
	 * @throws BioException
	 */
	public String formatOutput() throws BioException {
		return formatOutput(60);
	}

	/**
	 * This method provides a BLAST-like formated alignment from the given
	 * <code>String</code>s, in which the sequence coordinates and the
	 * information "Query" or "Sbjct", respectively is added to each line. Each
	 * line contains <code>width</code> sequence characters including the gap
	 * symbols plus the meta information. There is one white line between two
	 * pairs of sequences.
	 * 
	 * @param width
	 *            the number of symbols to be displayed per line.
	 * @return formated String.
	 * @throws BioException
	 */
	public String formatOutput(int width) throws BioException {
		int i, j;
		/*
		 * Highlights equal symbols within the alignment, String match/missmatch
		 * representation
		 */
		StringBuffer path = new StringBuffer();
		for (i = 0; i < Math.min(queryEnd - queryStart, subjectEnd
				- subjectStart) + 1; i++) {
			Symbol a = query.symbolAt(i + queryStart);
			Symbol b = subject.symbolAt(i + subjectStart);
			if (!a.equals(query.getAlphabet().getGapSymbol())
					&& !b.equals(subject.getAlphabet().getGapSymbol())
					&& (subMatrix.getValueAt(a, b) >= 0 || a.equals(b)))
				path.append('|');
			else
				path.append(' ');
		}

		int maxLength = path.length();
		/*
		 * Math.max(queryEnd - queryStart, subjectEnd - subjectStart) + 1;
		 */
		Formatter output = new Formatter();
		output.format("%n Time (ms):  %s%n", time);
		output.format(" Length:     %d%n", maxLength);
		output.format("  Score:     %d%n", getScore());
		output.format("  Query:     %s, Length: %d%n", query.getName(), query
				.length()
				- nGapsQ);
		output.format("  Sbjct:     %s, Length: %d%n", subject.getName(),
				subject.length() - nGapsS);
		output.format(
				" Identities: %d/%d, i.e., %d %% (query) and %d %% (sbjct)%n",
				identicals, maxLength, Math.round(getPercentIdentityQuery()),
				Math.round(getPercentIdentitySubject()));
		output.format(
				" Similars:   %d/%d, i.e., %d %% (query) and %d %% (sbjct)%n",
				similars, maxLength, Math.round(getPercentSimilarityQuery()),
				Math.round(getPercentSimilaritySubject()));
		output.format(
				" No. gaps:   %d (%d %%) in query and %d (%d %%) in sbjct%n",
				nGapsQ, Math.round(getPercentGapsQuery()), nGapsS, Math
						.round(getPercentGapsTarget()));

		int queryLPos = queryStart, queryRPos, pathLPos = 0, pathRPos;
		int subjectLPos = subjectStart, subjectRPos;
		int ql = queryLPos - 1, qr = queryLPos - 1, qgaps;
		int sl = subjectLPos - 1, sr = subjectLPos - 1, sgaps;

		int widthLeft = String.valueOf((int) Math.max(queryStart, queryEnd))
				.length();
		int widthRight = String.valueOf((int) Math.max(queryEnd, subjectEnd))
				.length() + 1;

		// Take width of the meta information into account.
		width = Math.max(width - widthLeft - widthRight - 12, 2);

		for (i = 1; i <= Math.ceil((double) maxLength / width); i++) {

			// Query
			queryRPos = Math.min(queryStart + i * width - 1, Math.min(queryEnd,
					subjectEnd - subjectStart + queryStart));
			qgaps = 0;
			for (j = queryLPos; j <= queryRPos; j++)
				if (!query.symbolAt(j).equals(
						query.getAlphabet().getGapSymbol()))
					qr++;
				else
					qgaps++;
			if (qgaps <= queryRPos - queryLPos)
				ql++;
			output.format("%nQuery:   %" + widthLeft + "d ", ql);
			output.format("%s ", query.subStr(queryLPos, queryRPos));
			output.format("%-" + widthRight + "d%n", qr);
			queryLPos = queryRPos + 1;
			ql = qr;

			// Path
			pathRPos = Math.min(i * width, path.length());
			output.format("%-" + (widthLeft + 10) + "c%s", Character
					.valueOf(' '), path.substring(pathLPos, pathRPos));
			pathLPos = pathRPos;

			// Sbjct
			subjectRPos = Math.min(subjectStart + i * width - 1, Math.min(
					queryEnd - queryStart + subjectStart, subjectEnd));
			sgaps = 0;
			for (j = subjectLPos; j <= subjectRPos; j++)
				if (!subject.symbolAt(j).equals(
						subject.getAlphabet().getGapSymbol()))
					sr++;
				else
					sgaps++;
			if (sgaps <= subjectRPos - subjectLPos)
				sl++;
			output.format("%nSbjct:   %" + widthLeft + "d ", sl);
			output.format("%s ", subject.subStr(subjectLPos, subjectRPos));
			output.format("%-" + widthRight + "d%n", sr);
			subjectLPos = subjectRPos + 1;
			sl = sr;
		}

		return output.toString();
	}

	/**
	 * 
	 * @return
	 */
	public float getPercentGapsTarget() {
		return nGapsS / (float) subject.length() * 100;
	}

	/**
	 * 
	 * @return
	 */
	public float getPercentGapsQuery() {
		return nGapsQ / (float) query.length() * 100;
	}

}
