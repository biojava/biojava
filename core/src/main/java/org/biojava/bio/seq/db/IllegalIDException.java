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

import org.biojava.bio.BioException;

/**
 * @author Matthew Pocock
 * @author Keith James
 */
public class IllegalIDException extends BioException {
  public IllegalIDException() {
    super();
  }

  public IllegalIDException(Throwable t) {
    super(t);
  }

  public IllegalIDException(String message) {
    super(message);
  }

  public IllegalIDException(Throwable t, String message) {
    super(message, t);
  }
}
