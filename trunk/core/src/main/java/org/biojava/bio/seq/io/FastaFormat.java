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
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ParseErrorEvent;
import org.biojava.utils.ParseErrorListener;


/**
 * Format object representing FASTA files. These files are almost pure
 * sequence data. The only `sequence property' reported by this parser
 * is PROPERTY_DESCRIPTIONLINE, which is the contents of the
 * sequence's description line (the line starting with a '>'
 * character). Normally, the first word of this is a sequence ID. If
 * you wish it to be interpreted as such, you should use
 * FastaDescriptionLineParser as a SeqIO filter.
 *
 * If you pass it a RichSeqIOListener, you'll get RichSequence objects
 * in return. Likewise, if you write RichSequence objects, you'll get
 * absolutely correct FASTA formatted output.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Greg Cox
 * @author Lukas Kall
 * @author Richard Holland
 * @author Mark Schreiber
 * @deprecated Use org.biojavax.bio.seq.io.FastaFormat
 */

public class FastaFormat implements SequenceFormat,
        Serializable,
        org.biojava.utils.ParseErrorListener,
        org.biojava.utils.ParseErrorSource {
    public static final String DEFAULT = "FASTA";
    
    /**
     * Constant string which is the property key used to notify
     * listeners of the description lines of FASTA sequences.
     */
    public final static String PROPERTY_DESCRIPTIONLINE = "description_line";
    
    protected Vector mListeners = new Vector();
    
    /**
     * The line width for output.
     */
    protected int lineWidth = 60;
    
    /**
     * Retrive the current line width.
     *
     * @return the line width
     */
    public int getLineWidth() {
        return lineWidth;
    }
    
    /**
     * Set the line width.
     * <p>
     * When writing, the lines of sequence will never be longer than the line
     * width.
     *
     * @param width the new line width
     */
    public void setLineWidth(int width) {
        this.lineWidth = width;
    }
    
    /**
     * Reads information from a flatfile to a <code>SeqIOListener</code>
     * using a <code>SymbolTokenizer</code> to convert sequence strings
     * to <code>Symbol</code> objects.
     * @param reader The reader that is the source of the information
     * @param symParser converts text seqeunce to biojava objects
     * @param siol The listener that listens for event callbacks from this class.
     * The listener can be a <code>RichSeqIOListener</code>.
     * @throws org.biojava.bio.symbol.IllegalSymbolException if <code>symParser</code>
     * doesn't know how to convert the text sequence into biojava <code>Symbol</code>s 
     * @throws java.io.IOException if there is a problem reading.
     * @throws org.biojava.bio.seq.io.ParseException if the source cannot be parsed.
     * @return true if there is another unread sequence in the source.
     */
    public boolean readSequence(
            BufferedReader reader,
            SymbolTokenization symParser,
            SeqIOListener siol
            )	throws
            IllegalSymbolException,
            IOException,
            ParseException {
        String line = reader.readLine();
        if (line == null) {
            throw new IOException("Premature stream end");
        }
        while(line.length() == 0) {
            line = reader.readLine();
            if (line == null) {
                throw new IOException("Premature stream end");
            }
        }
        if (!line.startsWith(">")) {
            throw new IOException("Stream does not appear to contain FASTA formatted data: " + line);
        }
        
        siol.startSequence();
        
        String description = line.substring(1).trim();
        
        String regex = "(\\S+)(\\s+(.*))*";
        Pattern p = Pattern.compile(regex);
        Matcher m = p.matcher(description);
        if (!m.matches()) {
            throw new IOException("Stream does not appear to contain FASTA formatted data: " + line);
        }
        
        String name = m.group(1);
        
        siol.setName(name);
        siol.addSequenceProperty(PROPERTY_DESCRIPTIONLINE, description);
        
        boolean seenEOF = readSequenceData(reader, symParser, siol);
        siol.endSequence();
        
        return !seenEOF;
    }
    
    private boolean readSequenceData(
            BufferedReader r,
            SymbolTokenization parser,
            SeqIOListener listener
            ) throws
            IOException,
            IllegalSymbolException {
        char[] cache = new char[512];
        boolean reachedEnd = false, seenEOF = false;
        StreamParser sparser = parser.parseStream(listener);
        
        while (!reachedEnd) {
            r.mark(cache.length + 1);
            int bytesRead = r.read(cache, 0, cache.length);
            if (bytesRead < 0) {
                reachedEnd = seenEOF = true;
            } else {
                int parseStart = 0;
                int parseEnd = 0;
                while (!reachedEnd && parseStart < bytesRead && cache[parseStart] != '>') {
                    parseEnd = parseStart;
                    
                    while (parseEnd < bytesRead &&
                            cache[parseEnd] != '\n' &&
                            cache[parseEnd] != '\r'
                            ) {
                        ++parseEnd;
                    }
                    
                    sparser.characters(cache, parseStart, parseEnd - parseStart);
                    
                    parseStart = parseEnd + 1;
                    while (parseStart < bytesRead &&
                            (cache[parseStart] == '\n' ||
                            cache[parseStart] == '\r') ) {
                        ++parseStart;
                    }
                }
                if (parseStart < bytesRead && cache[parseStart] == '>') {
                    try {
                        r.reset();
                    } catch (IOException ioe) {
                        throw new IOException(
                                "Can't reset: " +
                                ioe.getMessage() +
                                " parseStart=" + parseStart +
                                " bytesRead=" + bytesRead
                                );
                    }
                    if (r.skip(parseStart) != parseStart) {
                        throw new IOException("Couldn't reset to start of next sequence");
                    }
                    reachedEnd = true;
                }
            }
        }
        
        sparser.close();
        return seenEOF;
    }
    
    /**
     * Return a suitable description line for a Sequence. If the
     * sequence's annotation bundle contains PROPERTY_DESCRIPTIONLINE,
     * this is used verbatim.  Otherwise, the sequence's name is used.
     */
    protected String describeSequence(Sequence seq) {
        String description = null;
        
        Annotation seqAnn = seq.getAnnotation();
        
        if(seqAnn.containsProperty(PROPERTY_DESCRIPTIONLINE)) {
            description = (String) seqAnn.getProperty(PROPERTY_DESCRIPTIONLINE);
        } else {
            description = seq.getName();
        }
        
        return description;
    }
    
    /**
     * Writes a <code>Sequence</code> or <code>RichSequence</code> to a 
     * <code>PrintStream</code> in FASTA format. If the sequence is a 
     * <code>RichSequence</code> the format of the header will be in line with
     * the NCBI standard.
     * @param seq the sequence to format
     * @param os the stream to write the sequence to. To print to screen use
     * <code>System.out</code>
     * @throws java.io.IOException if data cannot be written to <code>os</code>
     */
    public void writeSequence(Sequence seq, PrintStream os)
    throws IOException {
        os.print(">");
        os.println(describeSequence(seq));
        
        int length = seq.length();
        
        for (int pos = 1; pos <= length; pos += lineWidth) {
            int end = Math.min(pos + lineWidth - 1, length);
            os.println(seq.subStr(pos, end));
        }
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
    public void writeSequence(Sequence seq, String format, PrintStream os)
    throws IOException {
        if (! format.equalsIgnoreCase(getDefaultFormat()))
            throw new IllegalArgumentException("Unknown format '"
                    + format
                    + "'");
        writeSequence(seq, os);
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
    public synchronized void removeParseErrorListener(ParseErrorListener theListener) {
        if (mListeners.contains(theListener) == true) {
            mListeners.removeElement(theListener);
        }
    }
    
    /**
     * This method determines the behaviour when a bad line is processed.
     * Some options are to log the error, throw an exception, ignore it
     * completely, or pass the event through.
     * <p>
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
        
        for (int index = 0; index < listeners.size(); index++) {
            ParseErrorListener client = (ParseErrorListener)listeners.elementAt(index);
            client.BadLineParsed(theEvent);
        }
    }
}
