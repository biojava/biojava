package org.biojava.nbio.core.exceptions;
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
 * Created on DATE
 *
 */
/**
 * General abstraction of different parsing errors
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class ParserException extends RuntimeException {

  private static final long serialVersionUID = -4101924035353204493L;

  public ParserException(String message) {
    super(message);
  }

  public ParserException(Exception e) {
    super(e);
  }

  public ParserException(String message, Exception e) {
    super(message, e);
  }

}
