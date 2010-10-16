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

import java.util.Map;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.symbol.SymbolList;

/**
 * Objects of this type represent one particular installation (not
 * just implementation) of a sequence similarity searcher such as
 * BLASTP.  Objects of this type must be immutable such that all
 * parameters affecting the outcome of a search must be passed as an
 * argument to the search() method.
 *
 * @author <A href="mailto:Gerald.Loeffler@vienna.at">Gerald
 * Loeffler</A> for the <A href="http://www.imp.univie.ac.at">IMP</A>
 *
 */
public interface SeqSimilaritySearcher
{
    /**
     * Return a set of all databases that can be searched with this
     * sequence similarity searcher. This method reflects the fact
     * that objects of this type represent a particular installation
     * of a searcher which has access to only a limited number of
     * databases.
     *
     * @return an unmodifiable Set of SequenceDB objects. Never
     * returns null but may return an empty set.
     */
    Set getSearchableDBs();
  
    /**
     * Using this sequence similarity searcher, search with the given
     * sequence against the given sequence database. It is a
     * precondition that the sequence and the sequence database be
     * compatible with each other and with this sequence similarity
     * searcher, otherwise an IllegalArgumentException is thrown.
     *
     * <p> Particular implementations of a searcher will differ in the
     * number and kind of parameters that can be used to customize a
     * search. All these parameters must be passed as an argument to
     * this method, i.e. these parameters are not part of the state of
     * this searcher.
     *
     * <p> This method performs a synchronous search, i.e. it will
     * block until the search result becomes available.
     *
     * @param querySeq the sequence with which to search.  May not be
     * null otherwise an IllegalArgumentException is thrown.
     * @param db the sequence database against which the similarity
     * search will be performed.  May not be null otherwise an
     * IllegalArgumentException is thrown. Must also be an element of
     * the set of searchable dbs returned by getSearchableDBs().
     * @param searchParameters parameters that customize the
     * search. Null must always be a legal value for this argument and
     * results in a default search being performed. If this map
     * contains keys and/or values that are not supported by a
     * particular implementation of this interface, an
     * IllegalArgumentException is thrown.
     *
     * @return the sequence similarity search result that encapsulates
     * the result of the search. Never returns null.
     *
     * @exception BioException if the actual search fails for some
     * reason. If, however, the search can not even be started,
     * because a precondition is not met, an IllegalArgumentException
     * is thrown.
     */
    SeqSimilaritySearchResult search(SymbolList querySeq,
				     SequenceDB db,
				     Map        searchParameters) 
	throws BioException;
}
