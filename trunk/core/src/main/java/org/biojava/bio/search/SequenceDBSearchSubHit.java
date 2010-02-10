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

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ObjectUtil;

/**
 * <p><code>SequenceDBSearchSubHit</code> objects represent sub-hits
 * which make up a hit. In the case of Blast these correspond to
 * HSPs.</p>
 *
 * @author Keith James
 * @since 1.1
 * @deprecated SimpleSeqSimilaritySearchSubHit has been made
 * Annotatable and is now functionally identical.
 * @see SeqSimilaritySearchSubHit
 */
public class SequenceDBSearchSubHit  extends AbstractChangeable
    implements SeqSimilaritySearchSubHit
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
    private Alignment  alignment;
    private Annotation annotation;

    // Hashcode is cached after first calculation because the data on
    // which is is based do not change
    private int hc;
    private boolean hcCalc;

    /**
     * Creates a new <code>SequenceDBSearchSubHit</code> object.
     *
     * @param queryStart an <code>int</code> value indicating the
     * start coordinate of the hit on the query sequence.
     * @param queryEnd an <code>int</code> value indicating the end
     * coordinate of the hit on the query sequence.
     * @param queryStrand a <code>Strand</code> object indicating the
     * strand of the hit with respect to the query sequence, which may
     * be null for protein similarities.
     * @param subjectStart an <code>int</code> value indicating the
     * start coordinate of the hit on the subject sequence.
     * @param subjectEnd an <code>int</code> value indicating the end
     * coordinate of the hit on the query sequence.
     * @param subjectStrand a <code>Strand</code> object indicating
     * the strand of the hit with respect to the query sequence, which
     * may be null for protein similarities.
     * @param score a <code>double</code> value; the score of the
     * subhit, which may not be NaN.
     * @param eValue a <code>double</code> the E-value of the
     * subhit, which may be NaN.
     * @param pValue a <code>double</code> value; the P-value of the
     * hit, which may be NaN.
     * @param alignment an <code>Alignment</code> object containing
     * the alignment described by the subhit region, which may not be
     * null.
     * @param annotation an <code>Annotation</code> object, which may
     * not be null.
     */
    public SequenceDBSearchSubHit(double    score,
                                  double    eValue,
                                  double    pValue,
                                  int       queryStart,
                                  int       queryEnd,
                                  Strand    queryStrand,
                                  int       subjectStart,
                                  int       subjectEnd,
                                  Strand    subjectStrand,
                                  Alignment alignment,
                                  Annotation annotation)
    {
        if (Double.isNaN(score))
        {
            throw new IllegalArgumentException("score was NaN");
        }

        // pValue may be NaN
        // eValue may be NaN
        if (alignment == null)
        {
            throw new IllegalArgumentException("alignment was null");
        }

        if (annotation == null)
        {
            throw new IllegalArgumentException("annotation was null");
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
        this.alignment     = alignment;
        this.annotation    = annotation;

        // Lock alignment by vetoing all changes
        this.alignment.addChangeListener(ChangeListener.ALWAYS_VETO);

        // Lock the annotation by vetoing all changes to properties
        annotation.addChangeListener(ChangeListener.ALWAYS_VETO);

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

    public Alignment getAlignment()
    {
        return alignment;
    }

    public Annotation getAnnotation()
    {
        return annotation;
    }

    public boolean equals(Object other)
    {
        if (other == this) return true;
        if (other == null) return false;

        if (! other.getClass().equals(this.getClass())) return false;

        SequenceDBSearchSubHit that = (SequenceDBSearchSubHit) other;

        if (! ObjectUtil.equals(this.score, that.score))
            return false;
        if (! ObjectUtil.equals(this.pValue, that.pValue))
            return false;
        if (! ObjectUtil.equals(this.eValue, that.eValue))
            return false;
        if (! ObjectUtil.equals(this.queryStart, that.queryStart))
            return false;
        if (! ObjectUtil.equals(this.queryEnd, that.queryEnd))
            return false;
        if (! ObjectUtil.equals(this.queryStrand, that.queryStrand))
            return false;
        if (! ObjectUtil.equals(this.subjectStart, that.subjectStart))
            return false;
        if (! ObjectUtil.equals(this.subjectEnd, that.subjectEnd))
            return false;
        if (! ObjectUtil.equals(this.subjectStrand, that.subjectStrand))
            return false;
        if (! ObjectUtil.equals(this.annotation, that.annotation))
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
            hc = ObjectUtil.hashCode(hc, queryStart);
            hc = ObjectUtil.hashCode(hc, queryEnd);
            hc = ObjectUtil.hashCode(hc, queryStrand);
            hc = ObjectUtil.hashCode(hc, subjectStart);
            hc = ObjectUtil.hashCode(hc, subjectEnd);
            hc = ObjectUtil.hashCode(hc, subjectStrand);
            hc = ObjectUtil.hashCode(hc, annotation);
            hcCalc = true;
        }

        return hc;
    }

    public String toString()
    {
        return "SequenceDBSearchSubHit with score " + getScore();
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
