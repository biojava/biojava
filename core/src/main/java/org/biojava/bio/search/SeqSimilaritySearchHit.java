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
import java.util.List;

import org.biojava.bio.Annotatable;
import org.biojava.bio.seq.StrandedFeature.Strand;

/**
 * Objects of this type represent one particular hit (sequence and
 * associated information) from a sequence similarity search.
 *
 * @author Gerald Loeffler
 * @author Keith James
 */
public interface SeqSimilaritySearchHit extends Annotatable
{
    /**
     * Return the overall score of this hit in the units defined by the
     * search algorithm.
     *
     * @return the overall score of this hit. This is a mandatory piece
     * of information and hence may not be NaN.
     */
    public double getScore();

    /**
     * Return the overall P-value of this hit.
     *
     * @return the overall P-value of this hit. This is an optional
     * (but desired) piece of information and implementations of this
     * interface may return NaN if a P-value is not available for this
     * hit.
     */
    public double getPValue();
  
    /**
     * Return the overall E-value of this hit.
     *
     * @return the overall E-value of this hit. This is an optional
     * (but desired) piece of information and implementations of this
     * interface may return NaN if an E-value is not available for
     * this hit.
     */
    public double getEValue();

    /**
     * Return the start position of the first sub-hit in the query
     * sequence.
     *
     * @return an <code>int</code>.
     */
    public int getQueryStart();

    /**
     * Return the end position of the last sub-hit in the query
     * sequence.
     *
     * @return an <code>int</code>.
     */
    public int getQueryEnd();

    /**
     * Return the strand of the hit with respect to the query
     * sequence. If the sub-hits are not all on the same strand this
     * should return the unknown strand. This may be null for protein
     * sequences.
     *
     * @return a <code>Strand</code>.
     */
    public Strand getQueryStrand();

    /**
     * Return the start position of the first sub-hit in the subject
     * sequence.
     *
     * @return an <code>int</code>.
     */
    public int getSubjectStart();

    /**
     * Return the end position of the last sub-hit in the subject
     * sequence.
     *
     * @return an <code>int</code>.
     */
    public int getSubjectEnd();

    /**
     * Return the strand of the sub-hit with respect to the subject
     * sequence. If the sub-hits are not all on the same strand this
     * should return the unknown strand. This may be null for protein
     * sequences.
     *
     * @return a <code>Strand</code>.
     */
    public Strand getSubjectStrand();

    /**
     * The sequence identifier of this hit within the sequence
     * database against which the search was performed.
     *
     * @return the (unique) sequence identifier for this hit, valid
     * within the sequence database against which this search was
     * performed. Never returns null.
     */
    public String getSubjectID();

    /**
     * Return all sub-hits for this sequence similarity search
     * hit. The sub-hits contain concrete alignments (and scores) for
     * sequence stretches from the sequence of this hit. The sub-hits
     * in the list returned by this method are sorted from best to
     * worst.
     * @return a List of SeqSimilaritySearchSubHit objects containing
     * all sub-hits for this hit.  Never returns null and the List is
     * guaranteed to contain at least one entry.
     */
    public List getSubHits();

    /**
     * <code>byScore</code> contains a
     * <code>SeqSimilaritySearchHit</code> comparator which
     * compares by their score.
     */
    public static final ByScoreComparator byScore = new ByScoreComparator();

    /**
     * <code>bySubHitCount</code> contains a
     * <code>SeqSimilaritySearchHit</code> comparator which
     * compares by their number of sub-hits.
     */
    public static final BySubHitCountComparator bySubHitCount =
        new BySubHitCountComparator();

    /**
     * <code>ByScoreComparator</code> compares
     * <code>SeqSimilaritySearchHit</code>s by their score.
     */
    public static final class ByScoreComparator implements Comparator
    {
        public int compare(Object o1, Object o2)
        {
            SeqSimilaritySearchHit h1 = (SeqSimilaritySearchHit) o1;
            SeqSimilaritySearchHit h2 = (SeqSimilaritySearchHit) o2;

            if (h1.getScore() > h2.getScore())
                return 1;
            else if (h1.getScore() < h2.getScore())
                return -1;
            else
                return 0;
        }
    }

    /**
     * <code>BySubHitCountComparator</code> compares
     * <code>SeqSimilaritySearchHit</code>s by their number of
     * sub-hits.
     */
    public static final class BySubHitCountComparator implements Comparator
    {
        public int compare(Object o1, Object o2)
        {
            SeqSimilaritySearchHit h1 = (SeqSimilaritySearchHit) o1;
            SeqSimilaritySearchHit h2 = (SeqSimilaritySearchHit) o2;

            if (h1.getSubHits().size() > h2.getSubHits().size())
                return 1;
            else if (h1.getSubHits().size() < h2.getSubHits().size())
                return -1;
            else
                return 0;
        }
    }
}
