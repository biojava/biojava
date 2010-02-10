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


package org.biojava.bio.seq.io;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.SequenceIterator;

/**
 * Writes all of the sequences from a SequenceIterator to a stream with a
 * particular format.
 * <p>
 * This can be wired from a StreamReader to make a simple file-format conversion
 * utility, or can be used to write out the sequences in a database to disk.
 *
 * <p>More functionality is offered by {@link org.biojavax.bio.seq.io.RichStreamWriter RichStreamWriter},
 * Use of this interface is prefered.</p>
 *
 * @author Matthew Pocock
 * @see org.biojavax.bio.seq.io.RichStreamWriter
 */
public class StreamWriter {
  /**
   * The format to write in.
   */
  private SequenceFormat format;
  
  /**
   * The stream to write to.
   */
  private PrintStream os;

  /**
   * Write each of the sequences in ss to the stream in the given format.
   *
   * @param ss  the SequenceIterator to loop over
   * @throws IOException if the stream has any problems
   */
  public void writeStream(SequenceIterator ss)
              throws IOException {
    while(ss.hasNext()) {
      try {
        format.writeSequence(ss.nextSequence(), os);
      } catch (BioException se) {
        se.printStackTrace();
      }
    }
  }

  /**
   * Generate a new StreamWriter to the stream os and using format.
   *
   * @param os  the OutputStream to write to
   * @param format the SequenceFormat to write with
   */
  public StreamWriter(OutputStream os, SequenceFormat format) {
    this.os = new PrintStream(os);
    this.format = format;
  }
}
