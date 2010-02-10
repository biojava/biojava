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
 * This is a no-frills implementation of the Index interface.
 * <p>
 * The file, start and ID are explicitly maintained as immutable properties of
 * the index. This implementation should be appropriate for many indexing
 * schemes. However, some schemes may wish to implement this interface as a
 * wrapper around a simple file offset, or an array index.
 *
 * @author Matthew Pocock
 */
public class SimpleIndex implements Index {
  private final File file;
  private final long start;
  private final int length;
  private final String id;
  
  /**
   * Build the index using the given file, start and id
   *
   * @param file the File this sequence is in
   * @param start how many bytes to skip to reach the first byte of the sequence
   * @param length how many bytes can be pulled out of the file to grab the record
   * @param id the ID of the sequence
   */
  public SimpleIndex(File file, long start, int length, String id) {
    this.file = file;
    this.start = start;
    this.length = length;
    this.id = id;
  }
  
  public File getFile() {
    return file;
  }
    
  public long getStart() {
    return start;
  }
  
  public int getLength() {
    return length;
  }

  public String getID() {
    return id;
  }
}
