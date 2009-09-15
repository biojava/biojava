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
package org.biojavax.bio.phylo.io.phylip;

import org.biojava.bio.seq.io.ParseException;

/**
 * Listens to events fired by the PHYLIP parser. Use these events to handle data
 * directly or construct objects.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */

public interface PHYLIPFileListener {
  
  /**
   * About to start a new file.
   */
  public void startFile();

  /**
   * Finished reading a file.
   */
  public void endFile() throws ParseException;

  /**
   * Set the number of sequences in the alignment.
   * 
   * @param count
   *        the expected number of sequences
   */
  public void setSequenceCount(int count);
  
  /**
   * Set the number of sites in the alignment
   * 
   * @param count
   *        the expected number of sites
   */
  public void setSitesCount(int count);
  
  /**
   * Set the name of the sequence which is about to be received.  If the name has already
   * been seen, the sequence should be appended.
   * 
   * @param name
   *        the label for the current sequence
   */
  public void setCurrentSequenceName(String name);
  
  /**
   * Receive sequence data for the current sequence.
   * 
   * @param sequence
   *        sequence text for the current sequence
   */
  public void receiveSequence(String sequence);

}
