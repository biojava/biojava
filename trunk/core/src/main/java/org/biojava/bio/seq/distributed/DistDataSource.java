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

package org.biojava.bio.seq.distributed;

import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;

/**
 * <p>Object which contributes data to a DistributedSequenceDB.</p>
 *
 * <p>DistDataSource is responsible for providing some information about what sequences exist,
 * what the SymbolList associated with it and what features are here. Typically, the objects
 * returned from DistributedSequenceDB will be composed from information from multiple
 * DistDataSource instances.</p>
 *
 * @author Thomas Down
 * @since 1.2
 *
 * Take an instance of this interface and add it to a DistributedSequenceDB.
 *
 * Implement this if you have information about some seqeunces but do not wish to or can
 * not integrate this with the main sequence at source. For example, if you have some
 * locally annotated features for SwissProt entries, you could create a DistDataSource
 * providing just your features and let the DistDataSource API integrate these in software
 * with a SwissProt sequence db provider.
 *
 * DistDataSource instances can provided sequence information and feature information. These
 * are integrated seperately. To provide sequences, implement hasSequence(), getSequence() and
 * ids(). ids(false).contains(id) should equal hasSequence(id). Features are provided by implementing
 * hasFeatures(), and the two getFeatures() methods. In the case where hasFeatures() returns true,
 * getFeatures() should return a FeatureHolder. If it is false, getFeatures() may raise a
 * BioException. If these rules are not followed, the results are undefined and may not be
 * consistent.
 */

public interface DistDataSource {
  /**
   * Find out if this DistDataSource provides the sequence information for a sequence ID.
   *
   * @param id  the String id of a sequence
   * @return true if this DistDataSource provides the primary sequence, false otherwise
   */
  public boolean hasSequence(String id) throws BioException;
  
  /**
   * Find out if this DistDataSource can provide features on a sequence with a particular ID.
   *
   * @param id  the String id of a sequence
   * @return true if this DistDataSource provides features for the sequence, false otherwise
   */
  public boolean hasFeatures(String id) throws BioException;

  /**
   * Get all features matching a FeatureFilter provided by this DistDataSource.
   * You can simulate getFeatures(id, ff, recurse) by using the advanced FeatureFilter
   * implementations.
   * @param ff  the FeatureFilter to search with
   * @return a FeatureHolder with all matching filters
   *
   **/
  public FeatureHolder getFeatures(FeatureFilter ff) throws BioException;
  
  /**
   * Get all features matching a FeatureFilter on a Sequence with an ID and recurse flats.
   * You can simulate getFeatures(ff) by adding the apropreate FeatureFilter implementations.
   *
   * @param id  the ID of the Sequence
   * @param ff  the FeatureFilter to search with
   * @param recurse true if we are to recurse the feature hierachy, false otherwise
   * @return a FeatureHolder containing all feature matching
   * @throws BioException if the features could not be fetched
   *
   */
  public FeatureHolder getFeatures(String id, FeatureFilter ff, boolean recurse) throws BioException;
  
  /**
   * Get a Sequence object for an ID.
   *
   * @param id  the ID of the Sequence to fetch
   * @return a Seqeunce if hasSequence(id) would return true
   * @throws BioException if either the ID could not be resolved or if the
   *         sequence could not be fetched
   */
  public Sequence getSequence(String id) throws BioException;
  
  /**
   * <p>Get the complete set of sequence IDs provided by this DistDataSource.</p>
   *
   * <p>If the recurse flat is true, the IDs associated with the top level will be returned.
   * However, if it is false, then IDs should be returned for all levels of an assembly
   * hierachy including the top level IDs.</p>
   *
   * @param topLevel  if true, return top level IDs, otherwise all IDs
   * @return a Set of String IDs
   * @throws BioException  if the IDs could not be fetched
   */
  public Set ids(boolean topLevel) throws BioException;
}
