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

import java.util.List;
import java.util.Map;

import org.biojava.bio.Annotatable;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.db.SequenceDB;

/**
 * Objects of this type represent one particular result of a sequence
 * similarity search.
 *
 * @author Gerald Loeffler
 * @author Keith James
 */
public interface SeqSimilaritySearchResult extends Annotatable
{
    /**
     * Returns the query sequence which was used to perform the search.
     *
     * @return the <code>Sequence</code> object used to search the
     * <code>SequenceDB</code>. Never returns null.
     */
    public Sequence getQuerySequence();

    /**
     * Returns the sequence database against which the search was
     * performed.
     *
     * @return the <code>SequenceDB object</code> against which the
     * search was carried out. Never returns null.
     */
    public SequenceDB getSequenceDB();

    /**
     * Returns the search parameters used in the search that produced
     * this search result.
     *
     * @return the (immutable) search parameter <code>Map
     * object</code>. May return null.
     */
    public Map getSearchParameters();

    /**
     * Return all hits in this sequence similarity search result. The
     * hits are sorted from best to worst.
     *
     * @return an (immutable) <code>List</code> of
     * <code>SeqSimilaritySearchHit</code> objects containing all hits
     * in the search result. Never returns null but may return an
     * empty list.
     */
    public List getHits();
}
