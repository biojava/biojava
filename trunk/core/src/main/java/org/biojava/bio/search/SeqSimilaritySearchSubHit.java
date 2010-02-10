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

package org.biojava.bio.search;

import java.util.Comparator;

import org.biojava.bio.Annotatable;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.seq.StrandedFeature.Strand;

/**
 * Objects of this type represent one particular sub-hit (one concrete
 * sequence stretch within a sequence and associated information) from
 * a sequence similarity search hit.
 *
 * @author Gerald Loeffler
 * @author Keith James
 */
public interface SeqSimilaritySearchSubHit extends Annotatable
{
    /**
     * This object is used as the label for the query sequence in the
     * alignment of the query sequence with this sub-hit sequence.
     */
    public static final String QUERY_LABEL = "Query";
  
    /**
     * Return the score of this sub-hit in the units defined by the
     * search algorithm.
     *
     * @return the score of this sub-hit. This is a mandatory piece of
     * information and hence may not be NaN.
     */
    public double getScore();

    /**
     * Return the P-value of this sub-hit.
     *
     * @return the P-value of this sub-hit. This is an optional (but
     * desired) piece of information and implementations of this
     * interface may return NaN if a P-value is not available for this
     * hit.
     */
    public double getPValue();
  
    /**
     * Return the E-value of this sub-hit.
     *
     * @return the E-value of this sub-hit. This is an optional (but
     * desired) piece of information and implementations of this
     * interface may return NaN if an E-value is not available for
     * this hit.
     */
    public double getEValue();
  
    /**
     * Return the start position of the sub-hit in the query sequence.
     *
     * @return an <code>int</code>.
     */
    public int getQueryStart();

    /**
     * Return the end position of the sub-hit in the query sequence.
     *
     * @return an <code>int</code>.
     */
    public int getQueryEnd();

    /**
     * Return the strand of the sub-hit with respect to the query
     * sequence. This may be null for protein sequences.
     *
     * @return a <code>Strand</code>.
     */
    public Strand getQueryStrand();

    /**
     * Return the start position of the sub-hit in the subject
     * sequence.
     *
     * @return an <code>int</code>.
     */
    public int getSubjectStart();

    /**
     * Return the start position of the sub-hit in the subject
     * sequence.
     *
     * @return an <code>int</code>.
     */
    public int getSubjectEnd();

    /**
     * Return the strand of the sub-hit with respect to the subject
     * sequence. This may be null for protein sequences.
     *
     * @return a <code>Strand</code>.
     */
    public Strand getSubjectStrand();

    /**
     * Return an alignment of (possibly part of) the query sequence
     * against (possibly part of) this hit sequence. In this
     * alignment, the query is identified by the label given by the
     * static field QUERY_LABEL.
     *
     * @return the alignment of the query sequence against this hit
     * sequence.
     */
    public Alignment getAlignment();

    /**
     * <code>byScore</code> contains a
     * <code>SeqSimilaritySearchSubHit</code> comparator which
     * compares by the score of the sub-hit.
     */
    public static final ByScoreComparator byScore = new ByScoreComparator();

    /**
     * <code>bySubjectStart</code> contains a
     * <code>SeqSimilaritySearchSubHit</code> comparator which
     * compares by the start position of the sub-hit on the subject
     * sequence.
     */
    public static final BySubjectStartComparator bySubjectStart
        = new BySubjectStartComparator();

    /**
     * <code>ByScoreComparator</code> compares
     * <code>SeqSimilaritySearchSubHit</code>s by their score.
     */
    public static final class ByScoreComparator implements Comparator
    {
        public int compare(Object o1, Object o2)
        {
            SeqSimilaritySearchSubHit h1 = (SeqSimilaritySearchSubHit) o1;
            SeqSimilaritySearchSubHit h2 = (SeqSimilaritySearchSubHit) o2;

            if (h1.getScore() > h2.getScore())
                return 1;
            else if (h1.getScore() < h2.getScore())
                return -1;
            else
                return 0;
        }
    }

    /**
     * <code>BySubjectStartComparator</code> compares
     * <code>SeqSimilaritySearchSubHit</code>s by their start position
     * on the subject sequence.
     */
    public static final class BySubjectStartComparator implements Comparator
    {
        public int compare(Object o1, Object o2)
        {
            SeqSimilaritySearchSubHit h1 = (SeqSimilaritySearchSubHit) o1;
            SeqSimilaritySearchSubHit h2 = (SeqSimilaritySearchSubHit) o2;

            if (h1.getSubjectStart() > h2.getSubjectStart())
                return 1;
            else if (h1.getSubjectStart() < h2.getSubjectStart())
                return -1;
            else
                return 0;
        }
    }
}
