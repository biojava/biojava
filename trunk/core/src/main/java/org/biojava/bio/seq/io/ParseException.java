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

import org.biojava.bio.BioException;

/**
 * ParseException should be thrown to indicate that there was a problem with
 * parsing sequence information.
 *
 * @author Matthew Pocock
 */
public class ParseException extends BioException {
  public ParseException() {
    super();
  }

  public ParseException(String message) {
    super(message);
  }

  public ParseException(Throwable nested) {
    super(nested);
  }

  public ParseException(Throwable nested, String message) {
    super(message, nested);
  }
  
  private static StringBuffer message;
  private static String intro = "\n\nA Exception Has Occurred During Parsing. \n"+
          "Please submit the details that follow to biojava-l@biojava.org or "+
          "post a bug report to http://bugzilla.open-bio.org/ \n\n";
  /**
   * Make a new error message.
   * @param format the format object that was doing the parsing
   * @param accession the accession number of the record that failed
   * @param identifier the identifier of the sequence that failed (eg the GI number for genbank)
   * @param comments any additional information
   * @param parseBlock the chunk of the file the parser was trying to parse when the error occured
   * @return the formatted error message
   */
  public static String newMessage(Class format, String accession, String identifier, String comments, String parseBlock){
      message = new StringBuffer();
      message.append(intro);
      message.append("Format_object=").append(format.getName()).append('\n');
      message.append("Accession=").append(accession).append('\n');
      message.append("Id=").append(identifier).append('\n');
      message.append("Comments=").append(comments).append('\n');
      message.append("Parse_block=").append(parseBlock).append('\n');
      message.append("Stack trace follows ....\n\n");
      return message.toString();
  }
}
