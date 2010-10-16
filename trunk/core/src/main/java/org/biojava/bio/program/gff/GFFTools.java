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

package org.biojava.bio.program.gff;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ParserException;

/**
 * @author Mark Schreiber
 * @author Matthew Pocock
 * @since 1.2
 */

public class GFFTools {

  /**
   * Flag to indicate that there is no score info.
   */
  public static double NO_SCORE = Double.NEGATIVE_INFINITY;

  /**
   * Flag to indicate that there is no frame info.
   */
  public static int NO_FRAME = -1;

  /**
   * Reads a <code>GFFEntrySet</code> from a file with no filtering.
   *
   * @param fileName the file containing the GFF
   * @throws FileNotFoundException if file is not found
   * @throws ParserException if format is wrong
   * @throws BioException if format is wrong
   * @throws IOException if file reading error occurs
   * @return a <code>GFFEntrySet</code> encapsulating the records read from the file
   * @deprecated use: readGff(File)
   */
  public static GFFEntrySet readGFF(String fileName)
    throws FileNotFoundException, ParserException, BioException, IOException
  {
    return readGFF(fileName, GFFRecordFilter.ACCEPT_ALL);
  }

  /**
   * Reads a GFFEntrySet from a file with the specified filter.
   *
   * @param fileName the file containing the GFF
   * @param recFilt the filter to use
   * @throws FileNotFoundException if file is not found
   * @throws ParserException if format is wrong
   * @throws BioException if format is wrong
   * @throws IOException if file reading error occurs
   * @return a <code>GFFEntrySet</code> encapsulating the records read from the file
   * @deprecated use: readGff(File,GFFRecordFilter)
   */
  public static GFFEntrySet readGFF(String fileName, GFFRecordFilter recFilt)
    throws FileNotFoundException, ParserException, BioException, IOException
  {
    GFFEntrySet gffEntries = new GFFEntrySet();
    GFFFilterer filterer = new GFFFilterer(gffEntries.getAddHandler(),recFilt);
    GFFParser parser = new GFFParser();
    parser.parse(new BufferedReader(new FileReader(fileName)),filterer);
    return gffEntries;
  }
  
 /**
  * Reads a <code>GFFEntrySet</code> from a file with no filtering.
  *
  * @param inFile the File containing the GFF
  * @throws FileNotFoundException if file is not found
  * @throws ParserException if format is wrong
  * @throws BioException if format is wrong
  * @throws IOException if file reading error occurs
  * @return a <code>GFFEntrySet</code> encapsulating the records read from the file
  */
  public static GFFEntrySet readGFF(File inFile)
    throws FileNotFoundException, ParserException, BioException, IOException
  {
    return readGFF(inFile, GFFRecordFilter.ACCEPT_ALL);
  }

  /**
   * Reads a GFFEntrySet from a file with the specified filter.
   *
   * @param inFile the File containing the GFF
   * @param recFilt the filter to use
   * @throws FileNotFoundException if file is not found
   * @throws ParserException if format is wrong
   * @throws BioException if format is wrong
   * @throws IOException if file reading error occurs
   * @return a <code>GFFEntrySet</code> encapsulating the records read from the file
   */
  public static GFFEntrySet readGFF(File inFile, GFFRecordFilter recFilt)
    throws FileNotFoundException, ParserException, BioException, IOException
  {
    GFFEntrySet gffEntries = new GFFEntrySet();
    GFFFilterer filterer = new GFFFilterer(gffEntries.getAddHandler(),recFilt);
    GFFParser parser = new GFFParser();
    parser.parse(new BufferedReader(new FileReader(inFile)),filterer);
    return gffEntries;
  }

  /**
   * Read all GFF entries from a buffered reader.
   *
   * This will read up untill the end of the reader.
   *
   * @param gffIn  the BufferedReader to read text from
   * @return a GFFEntrySet containing all of the GFF that could be read
   * @throws parserException  if the text could not be parsed as GFF
   * @throws BioException if there was some error reading the GFF
   * @throws IOException if there was an error with the reader
   */
  public static GFFEntrySet readGFF(BufferedReader gffIn)
    throws ParserException, BioException, IOException
  {
    return readGFF(gffIn, GFFRecordFilter.ACCEPT_ALL);
  }

  /**
   * Read all GFF entries matching a filter from a buffered reader.
   *
   * This will read up untill the end of the reader.
   *
   * @param gffIn  the BufferedReader to read text from
   * @return a GFFEntrySet containing all of the GFF that could be read
   * @throws parserException  if the text could not be parsed as GFF
   * @throws BioException if there was some error reading the GFF
   * @throws IOException if there was an error with the reader
   */
  public static GFFEntrySet readGFF(BufferedReader gffIn, GFFRecordFilter recFilt)
    throws ParserException, BioException, IOException
  {
    GFFEntrySet gffEntries = new GFFEntrySet();
    GFFFilterer filterer = new GFFFilterer(gffEntries.getAddHandler(),recFilt);
    GFFParser parser = new GFFParser();
    parser.parse(gffIn, filterer);
    return gffEntries;
  }

  /**
   * Writes a GFFEntrySet to a file.
   *
   * @param fileName the file to write to
   * @param ents the entries to write
   * @throws IOException if file writing fails
   */
  public static void writeGFF(String fileName, GFFEntrySet ents)
    throws IOException
  {
    PrintWriter pw = new PrintWriter(new FileWriter(fileName));
    writeGFF(pw, ents);
    pw.close();
  }
  
  /**
   * Writes a GFFEntrySet to a file.
   *
   * @param outFile  the file to write to
   * @param ents  the entry set to write
   * @throws IOException if writing to the file fails
   */
  public static void writeGFF(File outFile, GFFEntrySet ents)
    throws IOException
  {
    PrintWriter pw = new PrintWriter(new FileWriter(outFile));
    writeGFF(pw, ents);
    pw.close();
  }

  /**
   * Writes a GFFEntrySet to a PrintWriter.
   *
   * @param pw  the PrintWriter to write to
   * @param ents the entries to write
   * @throws IOException if file writing fails
   */
  public static void writeGFF(PrintWriter pw, GFFEntrySet ents)
    throws IOException
  {
    GFFWriter writer = new GFFWriter(pw);
    ents.streamRecords(writer);
  }

  /**
   * Annotates a sequence with the features from a GFF entry set with sequence
   * name matching this sequence.
   *
   * @param seq the <code>Sequence</code> to annotate.
   * @param ents the the GFF features to annotate it with.
   * @return a reference to a newly annotated sequence.
   */
  public static Sequence annotateSequence(Sequence seq, GFFEntrySet ents){
    Sequence annotated;
    try {
      annotated = ents.getAnnotator().annotate(seq);
    }
    catch (ChangeVetoException ex) {
      throw new BioError("Assertion Error: Unable to annotate sequence",ex);
    }catch (BioException ex) {
      throw new BioError("Assertion Error: Unable to annotate sequence",ex);
    }
    return annotated;
  }

  /**
   * Annotates a sequence with the features from a GFF entry set.
   *
   * @param seq the <code>Sequence</code> to annotate.
   * @param ents the the GFF features to annotate it with.
   * @param checkSeqName  boolean flat, if true only annotate sequence with
   *        features that have matching sequence names, otherwise annotate
   *        all features
   * @return a reference to a newly annotated sequence.
   */
  public static Sequence annotateSequence(
    Sequence seq,
    GFFEntrySet ents,
    boolean checkSeqName
  ) {
    Sequence annotated;
    try {
      annotated = ents.getAnnotator(checkSeqName).annotate(seq);
    }
    catch (ChangeVetoException ex) {
      throw new BioError("Assertion Error: Unable to annotate sequence",ex);
    }catch (BioException ex) {
      throw new BioError("Assertion Error: Unable to annotate sequence",ex);
    }
    return annotated;
  }

  /**
   * Annotates all sequences in a sequence DB with features from a GFF entry set.
   *
   * @param seqs  the SequenceDB to annotate
   * @param ents  the GFFEntrySet to annote with
   * @return a SequenceDB with all the annotations on
   */
  public static SequenceDB annotateSequences(SequenceDB seqs, GFFEntrySet ents)
    throws IllegalIDException, BioException{
    Set names = new HashSet();

    //get the list of names for each sequence
    for (Iterator i = ents.lineIterator(); i.hasNext(); ) {
      Object o = i.next();
      if(o instanceof GFFRecord){//only process GFFRecords not comments
        GFFRecord record = (GFFRecord)o;
        if(! names.contains(record.getSeqName())){
          names.add(record.getSeqName());
        }
      }
    }

    //filter entry set into subsets with same names, use that subset to annotate
    //the correct sequence.
    for (Iterator i = names.iterator(); i.hasNext(); ) {
      final String name = (String)i.next();
      GFFRecordFilter filt = new GFFRecordFilter(){
        public boolean accept(GFFRecord rec){
          return rec.getSeqName().equals(name);
        }
      };

      GFFEntrySet filtered = ents.filter(filt);
      Sequence seq = seqs.getSequence(name);
      seq = GFFTools.annotateSequence(seq, filtered);
    }

    return seqs;
  }

  /**
   * Creates a GFFEntrySet containing one entry for each feature on a sequence.
   *
   * @param seq  the Sequence to create features for
   * @return a new GFFEntrySet with gff records for each featre on the sequence
   * @throws BioException if something went wrong GFF-ifying the sequences
   *         features
   */
  public static GFFEntrySet gffFromSequence(Sequence seq)
  throws BioException {
    SequencesAsGFF sagff = new SequencesAsGFF();
    GFFEntrySet gffES = new GFFEntrySet();
    sagff.processSequence(seq, gffES.getAddHandler());
    return gffES;
  }
  
  /**
   * Creates a GFFEntrySet containing one entry for each feature on each
   * sequence of a SequenceDB.
   *
   * <p><em>Note:</em> This converts all features in the whole database to
   * in-memorey GFFRecord instances. This will take up considerable memory for
   * large databases.</p>
   *
   * @param seqDB  the SequenceDB to create features for
   * @return  a new GFFEntrySet with gff records for each feature on the database
   * @throws BioException if something went wrong GFF-ifying the sequences
   *         features
   */
public static GFFEntrySet gffFromSeqDB(SequenceDB seqDB)
  throws BioException {
    GFFEntrySet gffES = new GFFEntrySet();
    for(SequenceIterator si = seqDB.sequenceIterator(); si.hasNext(); ) {
      Sequence seq = si.nextSequence();
      SequencesAsGFF sagff = new SequencesAsGFF();
      sagff.processSequence(seq, gffES.getAddHandler());
    }
    return gffES;
  }
}
