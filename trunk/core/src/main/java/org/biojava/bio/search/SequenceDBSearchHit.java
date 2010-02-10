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

import java.util.Collections;
import java.util.List;

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ObjectUtil;

/**
 * <p><code>SequenceDBSearchHit</code> objects represent a similarity
 * search hit of a query sequence to a sequence referenced in a
 * SequenceDB object. The core data (score, E-value, P-value) have
 * accessors, while supplementary data are stored in the Annotation
 * object. Supplementary data are typically the more loosely formatted
 * details which vary from one search program to another (and between
 * versions of those programs).</p>
 *
 * <p>It is up to the user to define the meaning of the hit's
 * query/subject start/end/strand with respect to its constituent
 * sub-hits. One approach could be:</p>
 *
 * <ul>
 * <li>Hit query/subject start == start of first sub-hit</li>
 * <li>Hit query/subject   end == end of last sub-hit</li>
 * <li>Hit strand == POSITIVE if all sub-hits have strand POSITIVE</li>
 * <li>Hit strand == NEGATIVE if all sub-hits have strand NEGATIVE</li>
 * <li>Hit strand == UNKNOWN if sub-hits have either mixed or any UNKNOWN
 *     strands</li>
 * <li>Hit strand == null if the concept of strandedness is inappropriate
 *     for the sequence type i.e. for protein</li>
 * </ul>
 *
 * @author Keith James
 * @author Matthew Pocock
 * @since 1.1
 * @deprecated SimpleSeqSimilaritySearchHit has been made Annotatable
 * and is now functionally identical.
 * @see AbstractChangeable
 * @see SeqSimilaritySearchHit
 * @see Annotatable
 */
public class SequenceDBSearchHit extends AbstractChangeable
    implements SeqSimilaritySearchHit, Annotatable
{
    protected transient ChangeForwarder annotationForwarder;
    private double     score;
    private double     pValue;
    private double     eValue;
    private int        queryStart;
    private int        queryEnd;
    private Strand     queryStrand;
    private int        subjectStart;
    private int        subjectEnd;
    private Strand     subjectStrand;
    private String     subjectID;
    private Annotation annotation;
    private List       subHits;

    // Hashcode is cached after first calculation because the data on
    // which is is based do not change
    private int hc;
    private boolean hcCalc;

    /**
     * Creates a new <code>SequenceDBSearchHit</code> object.
     *
     * @param score a <code>double</code> value; the score of the hit,
     * which may not be NaN.
     * @param eValue a <code>double</code> value; the E-value of the
     * hit, which may be NaN.
     * @param pValue a <code>double</code> value; the P-value of the
     * hit, which may be NaN.
     * @param queryStart the start of the first sub-hit on the query
     * sequence.
     * @param queryEnd the end of the last sub-hit on the query
     * sequence.
     * @param queryStrand the strand of the sub-hits on the query
     * sequence, which may be null for protein similarities. If they
     * are not all positive or all negative, then this should be the
     * unknown strand.
     * @param subjectStart the start of the first sub-hit on the subject
     * sequence.
     * @param subjectEnd the end of the last sub-hit on the subject
     * sequence.
     * @param subjectStrand the strand of the sub-hits on the subject
     * sequence, which may be null for protein similarities. If they
     * are not all positive or all negative, then this should be the
     * unknown strand.
     * @param subjectID a <code>String</code> representing the ID in
     * the SequenceDB of the sequence which was hit, which may not be
     * null.
     * @param annotation an <code>Annotation</code> object, which may
     * not be null.
     * @param subHits a <code>List</code> object containing the
     * subhits, which may not be null. They should be sorted in the
     * order specified by the search program.
     */
    public SequenceDBSearchHit(double     score,
                               double     eValue,
                               double     pValue,
                               int        queryStart,
                               int        queryEnd,
                               Strand     queryStrand,
                               int        subjectStart,
                               int        subjectEnd,
                               Strand     subjectStrand,
                               String     subjectID,
                               Annotation annotation,
                               List       subHits)
    {
        if (Double.isNaN(score))
        {
            throw new IllegalArgumentException("score was NaN");
        }

        // pValue may be NaN
        // eValue may be NaN
        if (subjectID == null)
        {
            throw new IllegalArgumentException("subjectID was null");
        }

        if (annotation == null)
        {
            throw new IllegalArgumentException("annotation was null");
        }

        if (subHits == null)
        {
            throw new IllegalArgumentException("subHits was null");
        }

        this.score         = score;
        this.eValue        = eValue;
        this.pValue        = pValue;
        this.queryStart    = queryStart;
        this.queryEnd      = queryEnd;
        this.queryStrand   = queryStrand;
        this.subjectStart  = subjectStart;
        this.subjectEnd    = subjectEnd;
        this.subjectStrand = subjectStrand;
        this.subjectID     = subjectID;
        this.annotation    = annotation;
        this.subHits       = Collections.unmodifiableList(subHits);

        // Lock the annotation by vetoing all changes
        this.annotation.addChangeListener(ChangeListener.ALWAYS_VETO);

        hcCalc = false;
    }

    public double getScore()
    {
        return score;
    }

    public double getPValue()
    {
        return pValue;
    }

    public double getEValue()
    {
        return eValue;
    }

    public int getQueryStart()
    {
        return queryStart;
    }

    public int getQueryEnd()
    {
        return queryEnd;
    }

    public Strand getQueryStrand()
    {
        return queryStrand;
    }

    public int getSubjectStart()
    {
        return subjectStart;
    }

    public int getSubjectEnd()
    {
        return subjectEnd;
    }

    public Strand getSubjectStrand()
    {
        return subjectStrand;
    }

    public String getSubjectID()
    {
        return subjectID;
    }

    public List getSubHits()
    {
        return subHits;
    }

    /**
     * <code>getAnnotation</code> returns the Annotation associated
     * with this hit.
     *
     * @return an <code>Annotation</code>.
     */
    public Annotation getAnnotation()
    {
        return annotation;
    }

    public boolean equals(Object other)
    {
        if (other == this) return true;
        if (other == null) return false;

        if (! other.getClass().equals(this.getClass())) return false;

        SequenceDBSearchHit that = (SequenceDBSearchHit) other;

        if (! ObjectUtil.equals(this.score, that.score))
            return false;
        if (! ObjectUtil.equals(this.pValue, that.pValue))
            return false;
        if (! ObjectUtil.equals(this.eValue, that.eValue))
            return false;
        if (! ObjectUtil.equals(this.subjectID, that.subjectID))
            return false;
        if (! ObjectUtil.equals(this.subHits, that.subHits))
            return false;

        return true;
    }

    public int hashCode()
    {
        if (! hcCalc)
        {
            hc = ObjectUtil.hashCode(hc, score);
            hc = ObjectUtil.hashCode(hc, pValue);
            hc = ObjectUtil.hashCode(hc, eValue);
            hc = ObjectUtil.hashCode(hc, subjectID);
            hc = ObjectUtil.hashCode(hc, subHits);
            hcCalc = true;
        }

        return hc;
    }

    public String toString()
    {
        return "SequenceDBSearchHit to " + getSubjectID()
            + " with score " + getScore();
    }

    protected ChangeSupport getChangeSupport(ChangeType ct)
    {
        ChangeSupport cs = super.getChangeSupport(ct);

        if (annotationForwarder == null &&
            (ct.isMatchingType(Annotatable.ANNOTATION) || Annotatable.ANNOTATION.isMatchingType(ct)))
        {
            annotationForwarder =
                new ChangeForwarder.Retyper(this, cs, Annotation.PROPERTY);
            getAnnotation().addChangeListener(annotationForwarder,
                                              Annotatable.ANNOTATION);
        }

        return cs;
    }
}
