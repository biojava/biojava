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

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.GenbankFormat;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.seq.io.SequenceFormat;
import org.biojava.bio.symbol.Alphabet;

/**
 * @author Lei Lai
 * @author Matthew Pocock
 * @author Shuvankar Mukherjee
 */
public class GenpeptSequenceDB {
  private static SequenceFormat format = new GenbankFormat();
  private static String DBName = "Genpept";
  private boolean IOExceptionFound = false;
  private boolean ExceptionFound = false;

  protected SequenceFormat getSequenceFormat() {
    return format;
  }

  protected Alphabet getAlphabet() {
    return ProteinTools.getTAlphabet();
  }

  protected URL getAddress(String id) throws MalformedURLException {
    String defaultReturnFormat = "text";
    FetchURL seqURL = new FetchURL(DBName, defaultReturnFormat);
    String baseurl = seqURL.getbaseURL();
    String db = seqURL.getDB();

    String url = baseurl + db + "&id=" + id + "&rettype=gb";

    return new URL(url);
  }

  //ask user to input id and return format
  protected URL getAddress(String id, String format) throws
      MalformedURLException {
    FetchURL seqURL = new FetchURL(DBName, format);
    String baseurl = seqURL.getbaseURL();
    if (! (baseurl.equalsIgnoreCase("")))
      baseurl = seqURL.getbaseURL();
    String db = seqURL.getDB();
//	String returnFormat = seqURL.getReturnFormat();

    String url = baseurl + db + "&id=" + id + "&rettype=gb";

    return new URL(url);
  }

  public String getName() {
    return DBName;
  }

  public Sequence getSequence(String id) throws BioException {
    try {
      IOExceptionFound = false;
      ExceptionFound = false;
      URL queryURL = getAddress(id); //achieve URL based on ID


      //System.err.println("got data from " + queryURL);
      DataInputStream in = new DataInputStream(queryURL.openStream());
      BufferedReader reader = new BufferedReader(new InputStreamReader(in));
      SequenceIterator seqI = SeqIOTools.readGenpept(reader);

      return seqI.nextSequence();

    } catch (Exception e) {
      System.out.println(e.toString());
      IOExceptionFound = true;
      ExceptionFound = true;
      return null;
    }
  }

  public boolean checkIOException() {
    return IOExceptionFound;
  }

  public boolean checkException() {
    return ExceptionFound;
  }
}
