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

package	org.biojava.bio.seq.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.Vector;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ParseErrorEvent;
import org.biojava.utils.ParseErrorListener;


/**
 * Format reader for GenBank files. Converted from the old style io to
 * the new by working from <code>EmblLikeFormat</code>.
 *
 * @author Thomas Down
 * @author Thad	Welch
 * Added GenBank header	info to	the sequence annotation. The ACCESSION header
 * tag is not included.	Stored in sequence.getName().
 * @author Greg	Cox
 * @author Keith James
 * @author Matthew Pocock
 * @author Ron Kuhn
 * @deprecated Use org.biojavax.bio.seq.io.GenbankFormat
 */
public class GenbankFormat
        implements SequenceFormat,
        Serializable,
        org.biojava.utils.ParseErrorListener,
        org.biojava.utils.ParseErrorSource {
    public static final String DEFAULT = "GENBANK";
    
    protected static final String LOCUS_TAG = "LOCUS";
    protected static final String SIZE_TAG = "SIZE";
    protected static final String STRAND_NUMBER_TAG = "STRANDS";
    protected static final String TYPE_TAG = "TYPE";
    protected static final String CIRCULAR_TAG = "CIRCULAR";
    protected static final String DIVISION_TAG = "DIVISION";
    protected static final String DATE_TAG = "MDAT";
    
    protected static final String ACCESSION_TAG = "ACCESSION";
    protected static final String VERSION_TAG = "VERSION";
    protected static final String GI_TAG = "GI";
    protected static final String KEYWORDS_TAG = "KW";
    protected static final String DEFINITION_TAG = "DEFINITION";
    protected static final String SOURCE_TAG = "SOURCE";
    protected static final String ORGANISM_TAG = "ORGANISM";
    protected static final String REFERENCE_TAG = "REFERENCE";
    protected static final String COORDINATE_TAG = "COORDINATE";
    protected static final String REF_ACCESSION_TAG = "";
    protected static final String AUTHORS_TAG = "AUTHORS";
    protected static final String TITLE_TAG = "TITLE";
    protected static final String JOURNAL_TAG = "JOURNAL";
    protected static final String PUBMED_TAG = "PUBMED";
    protected static final String MEDLINE_TAG = "MEDLINE";
    protected static final String COMMENT_TAG = "COMMENT";
    protected static final String FEATURE_TAG = "FEATURES";
    protected static final String BASE_COUNT_TAG = "BASE";
    protected static final String FEATURE_FLAG = "FT";
    protected static final String START_SEQUENCE_TAG = "ORIGIN";
    protected static final String END_SEQUENCE_TAG = "//";
    
    protected static final String FEATURE_LINE_PREFIX = "     ";
    
    private Vector mListeners = new Vector();
    private boolean elideSymbols = false;
    
    /**
     * Reads a sequence from the specified reader using the Symbol
     * parser and Sequence Factory provided. The sequence read in must
     * be in Genbank format.
     *
     * @return boolean True if there is another sequence in the file; false
     * otherwise
     */
    public boolean readSequence(BufferedReader reader,
            SymbolTokenization symParser,
            SeqIOListener listener)
            throws IllegalSymbolException, IOException, ParseException {
        String line;
        boolean hasAnotherSequence    = true;
        boolean hasInternalWhitespace = false;
        
        GenbankContext ctx = new GenbankContext(symParser, listener);
        ctx.addParseErrorListener(this);
        ctx.setElideSymbols(this.getElideSymbols());
        
        listener.startSequence();
        
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(END_SEQUENCE_TAG)) {
                // To close the StreamParser encapsulated in the
                // GenbankContext object
                ctx.processLine(line);
                
                // Allows us to tolerate trailing whitespace without
                // thinking that there is another Sequence to follow
                while (true) {
                    reader.mark(1);
                    int c = reader.read();
                    
                    if (c == -1) {
                        hasAnotherSequence = false;
                        break;
                    }
                    
                    if (Character.isWhitespace((char) c)) {
                        hasInternalWhitespace = true;
                        continue;
                    }
                    
                    if (hasInternalWhitespace)
                        System.err.println("Warning: whitespace found between sequence entries");
                    
                    reader.reset();
                    break;
                }
                
                listener.endSequence();
                return hasAnotherSequence;
            }
            ctx.processLine(line);
        }
        
        throw new IOException("Premature end of stream for GENBANK");
    }
    
    public void	writeSequence(Sequence seq, PrintStream os)
    throws IOException {
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
    public void writeSequence(Sequence seq, String format, PrintStream os) throws IOException {
        SeqFileFormer former;
        
        if (format.equalsIgnoreCase("GENBANK"))
            former = new GenbankFileFormer();
        else if (format.equalsIgnoreCase("GENPEPT"))
            former = new GenpeptFileFormer();
        else if (format.equalsIgnoreCase("REFSEQ:PROTEIN"))
            former = new ProteinRefSeqFileFormer();
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
     * the default format.
     *
     * @return a <code>String</code>.
     * @deprecated
     */
    public String getDefaultFormat() {
        return DEFAULT;
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
    public synchronized void removeParseErrorListener(
            ParseErrorListener theListener) {
        if (mListeners.contains(theListener) == true) {
            mListeners.removeElement(theListener);
        }
    }
    
    /**
     * This method determines the behaviour when a bad line is processed.
     * Some options are to log the error, throw an exception, ignore it
     * completely, or pass the event through.
     * <P>
     * This method should be overwritten when different behavior is desired.
     *
     * @param theEvent The event that contains the bad line and token.
     */
    public void BadLineParsed(org.biojava.utils.ParseErrorEvent theEvent) {
        notifyParseErrorEvent(theEvent);
    }
    
// Protected methods
    /**
     * Passes the event on to all the listeners registered for ParseErrorEvents.
     *
     * @param theEvent The event to be handed to the listeners.
     */
    protected void notifyParseErrorEvent(ParseErrorEvent theEvent) {
        Vector listeners;
        synchronized(this) {
            listeners = (Vector)mListeners.clone();
        }
        
        int lnrCount = listeners.size();
        for (int index = 0; index < lnrCount; index++) {
            ParseErrorListener client = (ParseErrorListener)listeners.elementAt(index);
            client.BadLineParsed(theEvent);
        }
    }
    
    public boolean getElideSymbols() {
        return elideSymbols;
    }
    
    /**
     * Use this method to toggle reading of sequence data. If you're only
     * interested in header data set to true.
     * @param elideSymbols set to true if you don't want the sequence data.
     */
    public void setElideSymbols(boolean elideSymbols) {
        this.elideSymbols = elideSymbols;
    }
}

