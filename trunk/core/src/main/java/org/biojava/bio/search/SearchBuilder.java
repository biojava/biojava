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

import org.biojava.bio.BioException;

/**
 * The <code>SearchBuilder</code> interface is to be used by objects
 * which accumulate state via a <code>SearchContentHandler</code> and
 * then construct a <code>SeqSimilaritySearchResult</code> object.
 *
 * @author Keith James
 * @since 1.1
 * @see SearchContentHandler
 */
public interface SearchBuilder extends SearchContentHandler
{
    /**
     * The <code>makeSearchResult</code> method returns a
     * <code>SeqSimilaritySearchResult</code> instance created from
     * accumulated data.
     *
     * @return a <code>SeqSimilaritySearchResult</code>.
     *
     * @exception BioException if an error occurs.
     */
    public SeqSimilaritySearchResult makeSearchResult()
	throws BioException;
}
