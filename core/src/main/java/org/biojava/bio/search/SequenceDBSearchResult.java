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
import java.util.Map;

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ObjectUtil;

/**
 * <code>SequenceDBSearchResult</code> objects represent a result of a
 * search of a <code>SymbolList</code> against the sequences within a
 * <code>SequenceDB</code> object. The core data (query sequence,
 * database, search parameters, hits) have accessors, while
 * supplementary data are stored in the <code>Annotation</code>
 * object. Supplementary data are typically the more loosely formatted
 * details which vary from one search program to another (and between
 * versions of those programs).
 *
 * @author Keith James
 * @author Matthew Pocock
 * @since 1.1
 * @deprecated SimpleSeqSimilaritySearchResult has been made
 * Annotatable and is now functionally identical.
 * @see AbstractChangeable
 * @see SeqSimilaritySearchResult
 * @see Annotatable
 */
public class SequenceDBSearchResult extends AbstractChangeable
    implements SeqSimilaritySearchResult, Annotatable
{
    protected transient ChangeForwarder annotationForwarder;

    private Sequence   querySequence;
    private SequenceDB sequenceDB;
    private Map        searchParameters;
    private Annotation annotation;
    private List       hits;

    // Hashcode is cached after first calculation because the data on
    // which is is based do not change
    private int hc;
    private boolean hcCalc;

    /**
     * Creates a new <code>SequenceDBSearchResult</code>.
     *
     * @param querySequence a <code>Sequence</code>.
     * @param sequenceDB a <code>SequenceDB</code>.
     * @param searchParameters a <code>Map</code>.
     * @param annotation an <code>Annotation</code>.
     * @param hits a <code>List</code>.
     */
    public SequenceDBSearchResult(Sequence   querySequence,
                                  SequenceDB sequenceDB,
                                  Map        searchParameters,
                                  List       hits,
                                  Annotation annotation)
    {
        if (querySequence == null)
        {
            throw new IllegalArgumentException("querySequence was null");
        }

        if (sequenceDB == null)
        {
            throw new IllegalArgumentException("sequenceDB was null");
        }

        if (searchParameters != null)
        {
            this.searchParameters =
                Collections.unmodifiableMap(searchParameters);
        }

        if (annotation == null)
        {
            throw new IllegalArgumentException("annotation was null");
        }

        if (hits == null)
        {
            throw new IllegalArgumentException("hits was null");
        }

        // Lock the sequenceDB by vetoing all changes
        sequenceDB.addChangeListener(ChangeListener.ALWAYS_VETO);

        // Lock the querySeq by vetoing all changes
        querySequence.addChangeListener(ChangeListener.ALWAYS_VETO);

        // Lock the annotation by vetoing all changes to properties
        annotation.addChangeListener(ChangeListener.ALWAYS_VETO);

        this.querySequence = querySequence;
        this.sequenceDB    = sequenceDB;
        this.annotation    = annotation;
        this.hits          = Collections.unmodifiableList(hits);

        hcCalc = false;
    }

    public Sequence getQuerySequence()
    {
        return querySequence;
    }

    public SequenceDB getSequenceDB()
    {
        return sequenceDB;
    }

    public Map getSearchParameters()
    {
        return searchParameters;
    }

    public List getHits()
    {
        return hits;
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

        SequenceDBSearchResult that = (SequenceDBSearchResult) other;

        if (! ObjectUtil.equals(this.querySequence, that.querySequence))
            return false;
        if (! ObjectUtil.equals(this.sequenceDB, that.sequenceDB))
            return false;
        if (! ObjectUtil.equals(this.searchParameters, that.searchParameters))
            return false;
        if (! ObjectUtil.equals(this.annotation, that.annotation))
            return false;
        if (! ObjectUtil.equals(this.hits, that.hits))
            return false;

        return true;
    }

    public int hashCode()
    {
        if (! hcCalc)
        {
            hc = ObjectUtil.hashCode(hc, querySequence);
            hc = ObjectUtil.hashCode(hc, sequenceDB);
            hc = ObjectUtil.hashCode(hc, searchParameters);
            hc = ObjectUtil.hashCode(hc, hits);
            hc = ObjectUtil.hashCode(hc, annotation);
            hcCalc = true;
        }

        return hc;
    }

    public String toString()
    {
        return "SequenceDBSearchResult of " + getQuerySequence()
            + " against " + getSequenceDB().getName();
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
