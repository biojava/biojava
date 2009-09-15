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
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.Socket;
import java.net.URL;
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.GenbankFormat;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.seq.io.SequenceFormat;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.utils.ChangeVetoException;

/**
 * This class contains functions accessing DNA sequences in Genbank format.
 *
 * @author Lei Lai
 * @author Matthew Pocock
 * @author Laurent Jourdren
 * @author Shuvankar Mukherjee
 * @author Mark Schreiber
 * @author Richard Holland
 */
public class GenbankSequenceDB
{
  private static SequenceFormat format = new GenbankFormat();;//return format of the sequence
  private static String DBName="Genbank";//predefined the database name -- Genbank
  protected boolean IOExceptionFound=false;//check if IOException is found
  protected boolean ExceptionFound=false;//check if any exception is found
  protected static final String urlBatchSequences =
    "http://www.ncbi.nlm.nih.gov:80/entrez/eutils/efetch.fcgi";

  public GenbankSequenceDB() {
  }

  protected SequenceFormat getSequenceFormat()
  {
    return format;
  }

  protected Alphabet getAlphabet()
  {
    return DNATools.getDNA();
  }

  /**
   * Get the URL object for locating sequence object using eutils.
   * The default value of the return format of the sequence object is text.
   **/
  protected URL getAddress (String id) throws MalformedURLException
  {
        String defaultReturnFormat="text";
        FetchURL seqURL = new FetchURL(DBName, defaultReturnFormat);
        String baseurl = seqURL.getbaseURL();
        String db = seqURL.getDB();
        //String returnFormat = seqURL.getReturnFormat();

        String url = baseurl+db+"&id="+id+"&rettype=gb";

    return new URL (url);
  }

  /**
   * Get the URL object for locating sequence object using eutils.
   * User could specify the return format of the sequence object.
   */
  protected URL getAddress(String id, String format) throws MalformedURLException
  {
        FetchURL seqURL = new FetchURL(DBName, format);
        String baseurl = seqURL.getbaseURL();
        if (!(baseurl.equalsIgnoreCase("")))
                baseurl = seqURL.getbaseURL();
        String db = seqURL.getDB();
//	String returnFormat = seqURL.getReturnFormat();
//	String url = baseurl+db+"&"+returnFormat+"&id="+id;
        String url = baseurl+db+"&id="+id+"&rettype=gb";

    return new URL (url);
  }

  public String getName()
  {
    return DBName;
  }

  public Sequence getSequence(String id) throws Exception {
    try {
      IOExceptionFound = false;
      ExceptionFound = false;
      URL queryURL = getAddress(id); //get URL based on ID

      //  System.err.println("query is "+ queryURL.toString());

      //System.err.println("got data from " + queryURL);

      DataInputStream in = new DataInputStream(queryURL.openStream());
      BufferedReader reader = new BufferedReader(new InputStreamReader(in));
      SequenceIterator seqI = SeqIOTools.readGenbank(reader);

      return seqI.nextSequence();

    } catch (Exception e) {
      System.out.println("Exception found in GenbankSequenceDB -- getSequence");
      System.out.println(e.toString());
      ExceptionFound = true;
      IOExceptionFound = true;
      return null;
    }
  }

  public boolean checkIOException()
  {
        return IOExceptionFound;
  }

  public boolean checkException()
  {
        return ExceptionFound;
  }

  /**
   * Create the Http Post Request to fetch (in batch mode) a list of sequence
   * with Genbank.
   * @param url URL of the request
   * @param list List of sequence identifier
   * @return The Post request.
   */
  protected String makeBatchRequest(URL url, Set list) {

    StringBuffer params = new StringBuffer();
    params.append("db=nucleotide&rettype=gb&id=");

    boolean b = true;
    for (Iterator i = list.iterator(); b;) {
      String idSequence = (String) i.next();
      params.append(idSequence);
      if(i.hasNext()){
        params.append(",");
      }else{
        b =false;
        //params.append("\r\n");
      }
    }

    StringBuffer header = new StringBuffer();
    header.append("POST ");
    header.append(url.getPath());
    header.append(
      " HTTP/1.0\r\n"
        + "Connection: close\r\n"
        + "Accept: text/html, text/plain\r\n"
        + "Host: ");

    header.append(url.getHost());
    header.append(
      "\r\n"
        + "User-Agent: Biojava/GenbankSequenceDB\r\n"
        + "Content-Type: application/x-www-form-urlencoded\r\n"
        + "Content-Length: ");
    header.append(params.length());
    header.append("\r\n\r\n");

    StringBuffer request = new StringBuffer();
    request.append(header);
    request.append(params);

    return request.toString();
  }

  /**
   * Retrieve sequences from a Genbank
   *
   * @param list List of NCBI sequence number (GI), accession, accession.version,
   * fasta or seqid.
   * @return The database object (HashSequenceDB) with downloaded sequences.
   */
  public SequenceDB getSequences(Set list) throws BioException {

    return getSequences(list, null);
  }

  /**
   * Retrieve sequences from a Genbank
   *
   * @param list List of NCBI sequence number (GI), accession, accession.version,
   * fasta or seqid.
   * @param database Where to store sequences. if database is null, use an
   * HashSequenceDB Objet.
   * @return The database object with downloaded sequences.
   */
  public SequenceDB getSequences(Set list, SequenceDB database)
    throws BioException {

    if (database == null)
      database = new HashSequenceDB();

    try {

      URL url = new URL(urlBatchSequences);
      int port = url.getPort();
      String hostname = url.getHost();

      //Open the connection and the streams
      Socket s = new Socket(hostname, port);

      InputStream sin = s.getInputStream();
      BufferedReader fromServer =
        new BufferedReader(new InputStreamReader(sin));
      OutputStream sout = s.getOutputStream();
      PrintWriter toServer = new PrintWriter(new OutputStreamWriter(sout));

      // Put the Post request to the server
      toServer.print(makeBatchRequest(url, list));
      toServer.flush();

      // Delete response headers
      boolean finEntete = false;
      for (String l = null;
        ((l = fromServer.readLine()) != null) && (!finEntete);
        )
        if (l.equals(""))
          finEntete = true;

      SequenceIterator seqI = SeqIOTools.readGenbank(fromServer);

      while (seqI.hasNext())
        database.addSequence(seqI.nextSequence());

    } catch (MalformedURLException e) {
      throw new BioException(e,"Exception found in GenbankSequenceDB -- getSequences");
    } catch (IOException e) {
      throw new BioException(e,"Exception found in GenbankSequenceDB -- getSequences");
    } catch (BioException e) {
      throw new BioException(e,"Exception found in GenbankSequenceDB -- getSequences");
    } catch (ChangeVetoException e) {
      throw new BioException(e,"Exception found in GenbankSequenceDB -- getSequences");
    }

    return database;
  }
}
