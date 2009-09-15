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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Vector;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ParseErrorEvent;
import org.biojava.utils.ParseErrorListener;
import org.biojava.utils.ParseErrorSource;

/**
 * <p>
 * Format processor for handling EMBL records and similar files.  This
 * takes a very simple approach: all `normal' attribute lines are
 * passed to the listener as a tag (first two characters) and a value
 * (the rest of the line from the 6th character onwards).  Any data
 * between the special `SQ' line and the "//" entry terminator is
 * passed as a SymbolReader.
 * </p>
 *
 * <p>
 * This low-level format processor should normally be used in
 * conjunction with one or more `filter' objects, such as
 * EmblProcessor.
 * </p>
 *
 * <p>
 * Many ideas borrowed from the old EmblFormat processor by Thomas
 * Down and Thad Welch.
 * </p>
 *
 * @author Thomas Down
 * @author Greg Cox
 * @author Keith James
 * @author Len Trigg
 * @author Lorna Morris
 * @since 1.1
 * @deprecated Use org.biojavax.bio.seq.io.EMBLFormat instead
 */

public class EmblLikeFormat implements
    SequenceFormat,
    Serializable,
    ParseErrorSource,
    ParseErrorListener {
  public static final String DEFAULT = "EMBL";

  protected static final String ID_TAG = "ID";
  protected static final String SIZE_TAG = "SIZE";
  protected static final String STRAND_NUMBER_TAG = "STRANDS";
  protected static final String TYPE_TAG = "TYPE";
  protected static final String CIRCULAR_TAG = "CIRCULAR";
  protected static final String DIVISION_TAG = "DIVISION";
  protected static final String DR_TAG = "DR"; //Lorna: new tag

  protected static final String ACCESSION_TAG = "AC";
  protected static final String VERSION_TAG = "SV";
  protected static final String DATE_TAG = "DT";
  protected static final String DEFINITION_TAG = "DE";
  protected static final String KEYWORDS_TAG = "KW";
  protected static final String SOURCE_TAG = "OS";
  protected static final String ORGANISM_TAG = "OC";
  protected static final String ORGANISM_XREF_TAG = "OX";
  protected static final String REFERENCE_TAG = "RN";
  protected static final String COORDINATE_TAG = "RP";
  protected static final String REF_ACCESSION_TAG = "RX";
  protected static final String AUTHORS_TAG = "RA";
  protected static final String REF_XREF_TAG = "RX";
  protected static final String TITLE_TAG = "RT";
  protected static final String JOURNAL_TAG = "RL";
  protected static final String COMMENT_TAG = "CC";
  protected static final String FEATURE_TAG = "FH";
  protected static final String SEPARATOR_TAG = "XX";
  protected static final String FEATURE_TABLE_TAG = "FT";
  protected static final String START_SEQUENCE_TAG = "SQ";
  protected static final String END_SEQUENCE_TAG = "//";

  private boolean elideSymbols = false;
  private Vector mListeners = new Vector();

  /**
   * <p>Specifies whether the symbols (SQ) part of the entry should
   * be ignored. If this property is set to <code>true</code>, the
   * parser will never call addSymbols on the
   * <code>SeqIOListener</code>, but parsing will be faster if
   * you're only interested in header information.</p>
   *
   * <p> This property also allows the header to be parsed for files
   * which have invalid sequence data.</p>
   */
  public void setElideSymbols(boolean b) {
    elideSymbols = b;
  }

  /**
   * Return a flag indicating if symbol data will be skipped
   * when parsing streams.
   */
  public boolean getElideSymbols() {
    return elideSymbols;
  }

  public boolean readSequence(BufferedReader reader,
                              SymbolTokenization symParser,
                              SeqIOListener listener) throws
      IllegalSymbolException, IOException, ParseException {

    //EmblReferenceProperty reference = null; //lorna

    if (listener instanceof ParseErrorSource) {
      ( (ParseErrorSource) (listener)).addParseErrorListener(this);
    }

    String line;
    StreamParser sparser = null;
    boolean hasMoreSequence = true;
    boolean hasInternalWhitespace = false;

    listener.startSequence();

    while ( (line = reader.readLine()) != null) {
      if (line.startsWith(END_SEQUENCE_TAG)) {
        if (sparser != null) {
          // End of symbol data
          sparser.close();
          sparser = null;
        }

        // Allows us to tolerate trailing whitespace without
        // thinking that there is another Sequence to follow
        while (true) {
          reader.mark(1);
          int c = reader.read();

          if (c == -1) {
            hasMoreSequence = false;
            break;
          }

          if (Character.isWhitespace( (char) c)) {
            hasInternalWhitespace = true;
            continue;
          }

          if (hasInternalWhitespace)
            System.err.println(
                "Warning: whitespace found between sequence entries");

          reader.reset();
          break;
        }

        listener.endSequence();
        return hasMoreSequence;
      }
      else if (line.startsWith(START_SEQUENCE_TAG)) {
        // Adding a null property to flush the last feature;
        // Needed for Swissprot files because there is no gap
        // between the feature table and the sequence data
        listener.addSequenceProperty(SEPARATOR_TAG, "");

        sparser = symParser.parseStream(listener);
      }
      else {
        if (sparser == null) {
          // Normal attribute line
          String tag = line.substring(0, 2);
          String rest = null;
          if (line.length() > 5) {
            rest = line.substring(5);
          }

          if (tag.equals(REFERENCE_TAG)) { //only 1 reference_tag!

            try {
              //lorna added, tags read in order, when a complete set goes through,
              //spit out a single annotation event
              ReferenceAnnotation refAnnot = new ReferenceAnnotation();

              refAnnot.setProperty(tag, rest);
              while (! (tag.equals(SEPARATOR_TAG))) {
                // Normal attribute line

                line = reader.readLine();

                tag = line.substring(0, 2);

                if (line.length() > 5) {
                  rest = line.substring(5);
                }
                else {
                  rest = null; //for XX lines
                }

                if (refAnnot.containsProperty(tag)) {

                  Object property = refAnnot.getProperty(tag);
                  ArrayList properties;

                  if (property instanceof String) {
                    properties = new ArrayList();
                    properties.add(property);
                    properties.add(rest);
                    refAnnot.setProperty(tag, properties);
                  }
                  if (property instanceof ArrayList) {
                    ( (ArrayList) property).add(rest);
                  }
                }
                else {
                  refAnnot.setProperty(tag, rest);
                }
                //mark_s: required for parsing swissprot
                //fixme: it is actually possible to have more than one JOURNAL_TAG
                //so should really only break after the last one.
                if(tag.equals(JOURNAL_TAG))
                  break;
              }
              listener.addSequenceProperty(ReferenceAnnotation.class, refAnnot);

            } catch (ChangeVetoException cve) {
              cve.printStackTrace();
            }

          }
          // lorna, end
          else { //lorna
            listener.addSequenceProperty(tag, rest);
          } //lorna
        }
        else {
          // Sequence line
          if (!elideSymbols)
            processSequenceLine(line, sparser);
        }
      }
    }

    if (sparser != null)
      sparser.close();

    throw new IOException(
        "Premature end of stream or missing end tag '//' for EMBL");
  }

  /**
   * Dispatch symbol data from SQ-block line of an EMBL-like file.
   */
  protected void processSequenceLine(String line, StreamParser parser) throws
      IllegalSymbolException, ParseException {
    char[] cline = line.toCharArray();
    int parseStart = 0;
    int parseEnd = 0;

    while (parseStart < cline.length) {
      while (parseStart < cline.length && cline[parseStart] == ' ')
        ++
          parseStart;
      if (parseStart >= cline.length)
        break;

      if (Character.isDigit(cline[parseStart]))
        return;

      parseEnd = parseStart + 1;
      while (parseEnd < cline.length && cline[parseEnd] != ' ') {
        if (cline[parseEnd] == '.' || cline[parseEnd] == '~') {
          cline[parseEnd] = '-';
        }
        ++parseEnd;
      }

      // Got a segment of read sequence data
      parser.characters(cline, parseStart, parseEnd - parseStart);

      parseStart = parseEnd;
    }
  }

  public void writeSequence(Sequence seq, PrintStream os) throws IOException {
    writeSequence(seq, getDefaultFormat(), os);
  }

  /**
   * <code>writeSequence</code> writes a sequence to the specified
   * <code>PrintStream</code>, using the specified format.
   *
   * @param seq a <code>Sequence</code> to write out.
   * @param format a <code>String</code> indicating which sub-format
   * of those available from a particular
   * <code>SequenceFormat</code> implemention to use when
   * writing.
   * @param os a <code>PrintStream</code> object.
   *
   * @exception IOException if an error occurs.
   * @deprecated use writeSequence(Sequence seq, PrintStream os)
   */
  public void writeSequence(Sequence seq, String format, PrintStream os) throws
      IOException {
    SeqFileFormer former;

    if (format.equalsIgnoreCase("EMBL"))
      former = new EmblFileFormer();
    else if (format.equalsIgnoreCase("SWISSPROT"))
      former = new SwissprotFileFormer();
    else
      throw new IllegalArgumentException("Unknown format '"
                                         + format
                                         + "'");
    former.setPrintStream(os);

    SeqIOEventEmitter emitter =
        new SeqIOEventEmitter(GenEmblPropertyComparator.INSTANCE,
                              GenEmblFeatureComparator.INSTANCE);

    emitter.getSeqIOEvents(seq, former);
  }

  /**
   * <code>getDefaultFormat</code> returns the String identifier for
   * the default format written by a <code>SequenceFormat</code>
   * implementation.
   *
   * @return a <code>String</code>.
   * @deprecated
   */
  public String getDefaultFormat() {
    return DEFAULT;
  }

  /**
   * <p>
   * This method determines the behaviour when a bad line is processed.
   * Some options are to log the error, throw an exception, ignore it
   * completely, or pass the event through.
   * </p>
   *
   * <p>
   * This method should be overwritten when different behavior is desired.
   * </p>
   *
   * @param theEvent The event that contains the bad line and token.
   */
  public void BadLineParsed(ParseErrorEvent theEvent) {
    notifyParseErrorEvent(theEvent);
  }

  /**
   * Adds a parse error listener to the list of listeners if it isn't already
   * included.
   *
   * @param theListener Listener to be added.
   */
  public synchronized void addParseErrorListener(ParseErrorListener theListener) {
    if (mListeners.contains(theListener) == false) {
      mListeners.addElement(theListener);
    }
  }

  /**
   * Removes a parse error listener from the list of listeners if it is
   * included.
   *
   * @param theListener Listener to be removed.
   */
  public synchronized void removeParseErrorListener(ParseErrorListener
      theListener) {
    if (mListeners.contains(theListener) == true) {
      mListeners.removeElement(theListener);
    }
  }

  // Protected methods
  /**
   * Passes the event on to all the listeners registered for ParseErrorEvents.
   *
   * @param theEvent The event to be handed to the listeners.
   */
  protected void notifyParseErrorEvent(ParseErrorEvent theEvent) {
    Vector listeners;
    synchronized (this) {
      listeners = (Vector) mListeners.clone();
    }

    for (int index = 0; index < listeners.size(); index++) {
      ParseErrorListener client = (ParseErrorListener) listeners.elementAt(
          index);
      client.BadLineParsed(theEvent);
    }
  }
}