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


package org.biojava.bio.dp;

import org.biojava.bio.dist.Distribution;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * <p>
 * A state in a markov process that has an emission spectrum.
 * </p>
 *
 * <p>
 * These states have an associated Distribution. Within an HMM, these are the
 * states that actualy make your observed sequence. They also must supply
 * training behaviour to set the emission spectrum up.
 * </p>
 *
 * @author Matthew Pocock
 */
public interface EmissionState extends State, Trainable, Changeable {
  /**
   * <p>
   * This signals that the distribution associate with an EmissionState has
   * been altered.
   * </p>
   *
   * <p>
   * If the distribution has changed its weights, then the event'e
   * getChainedEvent method will return the event fired by the distribution. If
   * one distribution has been replaced by another, then the new and old
   * Distributions will be in current and previous, respectively.
   * </p>
   */
  public static final ChangeType DISTRIBUTION = new ChangeType(
    "The associated ditribution has changed",
    "org.biojava.bio.dp.EmissionState",
    "DISTRIBUTION"
  );
  
  /**
   * <p>
   * This signals that the advance array has been altered.
   * </p>
   *
   * <p>
   * current and previous should hold the current and previous advances,
   * respectively.
   * </p>
   */
  public static final ChangeType ADVANCE = new ChangeType(
    "The associated advance array has changed",
    "org.biojava.bio.dp.EmissionState",
    "ADVANCE"
  );
  
  /**
   * Determine the number of symbols this state advances along
   * one or more symbol lists.  In the simple case, this method
   * should almost always return {1} if it is a true `emmision'
   * state, or {0} if it is a dot state which only emits a gap
   * character.  For pairwise HMMs, it will normally return {1, 1}
   * for match state, and {0, 1} or {1, 0} for a gap state.  Under
   * some circumstances it may be valid to return values other
   * than 1 or 0, but you should consider the consequences for
   * HMM architecture very carefully, and contact the authors.
   *
   * Developers may wish to return a copy of some underlying array from this method
   * as code outside could modify the array you give
   */
  public int[] getAdvance();
  
  /**
   * Set the advance array.
   *
   * @param advance  an array of ints, specifying how many symbols are consumed
   *        from each sequence
   * @throws ChangeVetoException  if the implementation doesn't support setting
   *         advance, or if the change is vetoed
   */
  public void setAdvance(int[] advance) throws ChangeVetoException;
  
  /**
   * <p>
   * Get the Distribution associated with this state.
   * </p>
   *
   * <p>
   * If the state is to be added to an HMM, then the state's emission spectrum
   * must be compatible with the HMM - that is, their emission alphabets must
   * match.
   * </p>
   *
   * @return the current Distribution object used by this state
   */
  public Distribution getDistribution();
  
  /**
   * Set the Distribution associated with this state.
   * 
   * @param dis the new Distribution to use
   * @throws ChangeVetoException  if the implementation doesn't support setting
   *         the distribution, or if the change is vetoed
   */
  public void setDistribution(Distribution dis)
  throws ChangeVetoException;
}
