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

package org.biojava.bio.seq.db;

import java.io.File;

/**
 * This defines an index entry for an individual sequence within a set of
 * indexed files.
 *
 * @author Matthew Pocock
 */
public interface Index {
  /**
   * The file to retrieve from.
   *
   * @return the File containing this Sequence
   */
  File getFile();
  
  /**
   * Skipping this number of bytes through the file should put the file pointer
   * to the first byte of the sequence.
   *
   * @return the offset within the file
   */
  long getStart();
  
  /**
   * The entry can be slurped out of the file by grabbing length bytes from
   * start. If the length can't be read from a store then this method should
   * return -1.
   * 
   * @return the length in bytes of this indexed entry
   */
  int getLength();
  
  /**
   * The ID of the sequence at this position in this file.
   *
   * @return the ID of the indexed Sequence
   */
  String getID();
}
