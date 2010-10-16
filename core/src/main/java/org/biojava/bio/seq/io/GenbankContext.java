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

import java.util.StringTokenizer;
import java.util.Vector;

import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ParseErrorEvent;
import org.biojava.utils.ParseErrorListener;
import org.biojava.utils.ParseErrorSource;

/**
 * Encapsulate state used while	reading	data from a specific
 * Genbank file.
 *
 * @author Thomas Down
 * @author Greg Cox
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 */
class GenbankContext implements org.biojava.utils.ParseErrorListener, org.biojava.utils.ParseErrorSource
{
    private final static int HEADER = 1;
    private final static int FEATURES = 2;
    private final static int SEQUENCE = 3;
    private final static int TAG_LENGTH = 12;

    private int status;
    private SymbolTokenization symParser;
    private StreamParser streamParser;
    private String headerTag = "";
    private StringBuffer headerTagText = new StringBuffer();
    private SeqIOListener listener;
    private Vector mListeners = new Vector();
    private boolean elideSymbols;

    /**
     * Constructor that takes the listener and the Symbol parser from the
     * GenbankFormat.
     *
     * @param theSymbolParser Symbol parser to use in processing the file
     * @param theListener Listener to notify when field has been processed
     */
    protected GenbankContext(SymbolTokenization theSymbolParser,
			     SeqIOListener theListener)
    {
	this.status = HEADER;
	this.listener = theListener;

	this.symParser = theSymbolParser;
	this.streamParser = symParser.parseStream(listener);
         if (this.listener instanceof ParseErrorSource)
           ((ParseErrorSource)(this.listener)).addParseErrorListener(this);
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
    public void BadLineParsed(org.biojava.utils.ParseErrorEvent theEvent)
    {
        notifyParseErrorEvent(theEvent);
    }

    /**
     * Processes the line passed in in the context of previous lines processed.
     * This is the method that does the real work of the class.
     *
     * @throws ParseException Thrown when an error occurs parsing the file
     * @throws IllegalSymbolException Thrown when there is an illegal symbol
     * in the file; i.e. there is a symbol in the file that is not in the
     * alphabet being used
     * @param line The line to process
     */
    protected void processLine(String line)
	throws ParseException, IllegalSymbolException
    {
	if (line.startsWith(GenbankFormat.FEATURE_TAG))
	{
	    status = FEATURES;
	    this.saveSeqAnno();
	}
	else if (line.startsWith(GenbankFormat.START_SEQUENCE_TAG))
	{
	    status = SEQUENCE;
	    this.saveSeqAnno();
            // Additional commit to push the final feature off the stack
            headerTag = line;
            this.saveSeqAnno();
	}
	else if (line.startsWith(GenbankFormat.END_SEQUENCE_TAG))
	{
	    streamParser.close();
	}
	else if (status == FEATURES)
	{
	    processFeatureLine(line);
	}
	else if (status == SEQUENCE && elideSymbols == false)
	{
	    processSeqLine(line, streamParser);
	}
	else if (status == HEADER)
	{
	    processHeaderLine(line);
	}
    }

    /**
     * Adds a parse error listener to the list of listeners if it isn't already
     * included.
     *
     * @param theListener Listener to be added.
     */
    public synchronized void addParseErrorListener(
            ParseErrorListener theListener)
    {
        if (mListeners.contains(theListener) == false)
        {
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
            ParseErrorListener theListener)
    {
        if (mListeners.contains(theListener) == true)
        {
            mListeners.removeElement(theListener);
        }
    }

    /**
     * Passes the event on to all the listeners registered for ParseErrorEvents.
     *
     * @param theEvent The event to be handed to the listeners.
     */
    protected void notifyParseErrorEvent(ParseErrorEvent theEvent)
    {
        Vector listeners;
        synchronized(this)
        {
            listeners = (Vector) mListeners.clone();
        }

        int lnrCount = listeners.size();
        for (int index = 0; index < lnrCount; index++)
        {
            ParseErrorListener client =
                    (ParseErrorListener)listeners.elementAt(index);
            client.BadLineParsed(theEvent);
        }
    }


    /**
     * Private method that processes a line, assuming that it is a sequence
     * line instead of a feature line or header line.  Removes whitespace and
     * location numbers and hands the sequence data to theParser to convert to
     * biojava symbols
     *
     * @throws IllegalSymbolException Thrown when there is an illegal symbol
     * in the file; i.e. there is a symbol in the file that is not in the
     * alphabet being used
     * @param line The line to process
     * @param theParser Parser to convert the string data to biojava symbols
     */
    private void processSeqLine(String line, StreamParser theParser)
	throws IllegalSymbolException
    {
	char[] cline = line.toCharArray();
	int parseStart = 0;
	int parseEnd = 0;

	while (parseStart < cline.length)
	{
	    while ((parseStart < cline.length) &&
		  ((cline[parseStart] == ' ') ||
		   (Character.isDigit(cline[parseStart]))))
	    {
		// Read past leading spaces and numbers
		++parseStart;
	    }
	    if (parseStart >= cline.length)
	    {
		break;
	    }

	    parseEnd = parseStart + 1;
	    while ((parseEnd < cline.length && cline[parseEnd] != ' '))
	    {
                if (cline[parseEnd] == '.' || cline[parseEnd] == '~') {
                    cline[parseEnd] = '-';
                }
		++parseEnd;
	    }

	    // Got a segment of read sequence data
	    theParser.characters(cline, parseStart, parseEnd - parseStart);
	    parseStart = parseEnd;
	}
    }

    /**
     * Private method to process a line assuming it's a feature line. A
     * feature line is defined to be a line between the FEATURE tag and the
     * ORIGIN tag. The BASE COUNT line is processed here.
     *
     * @throws ParseException Thrown when an error occurs parsing the file
     * @param line The line to be processed
     */
    private void processFeatureLine(String line)
	throws ParseException
    {
	// Check the line is really a feature line
	if (line.startsWith(GenbankFormat.FEATURE_LINE_PREFIX))
	{
	    this.saveSeqAnno();
	    // Flag value as a feature line for GenbankProcessor. By a
	    // strange coincidence, this happens to be the same as the EMBL
	    // value
	    headerTag = GenbankFormat.FEATURE_FLAG;
	    headerTagText = new StringBuffer(line.substring(5));
	}
	else
	{
	    // Otherwise, process it as a header line
	    processHeaderLine(line);
	}
    }

    /**
     * Private method to process a line assuming it's a header line.  A header
     * line is defined to be a line before the FEATURE tag appears in the file.
     *
     * @throws ParseException Thrown when an error occurs parsing the file
     * @param line The line to be processed
     */
    private void processHeaderLine(String line)
        throws ParseException
    {
        if (line.startsWith(GenbankFormat.LOCUS_TAG))
        {
        	// Genbank changed the format of the Locus line for release 127.
        	// The new format is incompatible with the old.
        	if (this.isLocusLinePre127(line))
        	{
                    this.parseLocusLinePre127(line);
        	}
        	else
        	{
                    this.parseLocusLinePost127(line);
        	}
        }
        else if (line.startsWith(GenbankFormat.VERSION_TAG))
        {
            // VERSION line is a special case because it contains both
            // the VERSION field and the GI number
            this.saveSeqAnno();
            StringTokenizer lineTokens = new StringTokenizer(line);
            headerTag = lineTokens.nextToken();
            headerTagText = new StringBuffer(lineTokens.nextToken());

            if (lineTokens.hasMoreTokens()) {
                String nextToken = lineTokens.nextToken();
                if (nextToken.startsWith(GenbankFormat.GI_TAG))
                {
                    this.saveSeqAnno();
                    headerTag = GenbankFormat.GI_TAG; // Possibly should be UID?
                    headerTagText =
                        new StringBuffer(nextToken.substring(3));
                }
            }
        }
        else if (hasHeaderTag(line))
        {	// line	has a header tag
            this.saveSeqAnno();
            headerTag =	line.substring(0, TAG_LENGTH).trim();
            headerTagText = new StringBuffer(line.substring(TAG_LENGTH));
        }
        // gbpri1.seq (Release 125.0) has a line which is not
        // TAG_LENGTH long. Patch offered by Ron Kuhn (rkuhn@cellomics.com)
        else if (line.length() >= TAG_LENGTH)
        {	// keep	appending tag text value
            headerTagText.append(" " + line.substring(TAG_LENGTH));
        }
    }

    /**
     * Checks which version of the locus line format is used.  The algorithm
     * switches on the size of the line; <75 means pre-127, otherwise it's 127.
     *
     * @param theLine the line to check the format of.
     * @return TRUE if the line is in Genbank release 126 or earlier format.
     * FALSE otherwise
     */
    private boolean isLocusLinePre127(String theLine)
    {
    	return (theLine.length() < 75);
    }

    /**
     * Parses the locus line assuming it is in pre release 127 format.
     *
     * @param theLine Locus line to parse.
     * @throws ParseException If the line is too short.
     */
    private void parseLocusLinePre127(String theLine)
    	throws ParseException
    {
        if (theLine.length() < 73)
        {
            StringTokenizer locusTokens = new StringTokenizer(theLine);
            //Locus Tag
            locusTokens.nextToken();
            if (locusTokens.hasMoreTokens()) {
                saveSeqAnno2(GenbankFormat.LOCUS_TAG, locusTokens.nextToken());
            }
            else {
                throw new ParseException("LOCUS line too short [" + theLine + "]");
            }
        }
        else {
            saveSeqAnno2(GenbankFormat.LOCUS_TAG, theLine.substring(12, 22));
            saveSeqAnno2(GenbankFormat.SIZE_TAG, theLine.substring(22, 29));
            saveSeqAnno2(GenbankFormat.STRAND_NUMBER_TAG, theLine.substring(33, 35));
            saveSeqAnno2(GenbankFormat.TYPE_TAG, theLine.substring(36, 41));
            saveSeqAnno2(GenbankFormat.CIRCULAR_TAG, theLine.substring(42, 52));
            saveSeqAnno2(GenbankFormat.DIVISION_TAG, theLine.substring(52, 55));
            saveSeqAnno2(GenbankFormat.DATE_TAG, theLine.substring(62, 73));
        }
    }

    /**
     * Parses the locus line assuming it is in post release 127 format.
     * Will also handle the case where the strand tag is optional.  That
     * is the only tag that is supported as optional.  Awaiting a response from
     * NCBI if it is or if their data is incorrectly formatted.
     *
     * @param theLine Locus line to parse.
     * @throws ParseException If the line is too short.
     */
    private void parseLocusLinePost127(String theLine)
    	throws ParseException
    {
        if (theLine.length() < 79)
        {
            throw new ParseException("LOCUS line too short [" + theLine + "]");
        }

        StringTokenizer locusTokens = new StringTokenizer(theLine);
        if ((locusTokens.countTokens() == 8) || (locusTokens.countTokens() == 7))
        {
            // While Genbank 127 documentation doesn't allow the strand tag to
            // be optional, some files don't have it.  A key assumption here is
            // that this is the only tag that is optional.  The parser will
            // generate incorrect data if this assumption is violated.
            boolean includedStrandTag = (locusTokens.countTokens() == 8);

            // LOCUS tag; not stored
            locusTokens.nextToken();
            // Locus name
            saveSeqAnno2(GenbankFormat.LOCUS_TAG, locusTokens.nextToken());
            // Sequence length
            saveSeqAnno2(GenbankFormat.SIZE_TAG, locusTokens.nextToken());
            // "bp"; not stored
            locusTokens.nextToken();
            // Strand information
            // Both the strand and type are in the same token.  The strand
            // information is an optional part, so this is a bit hairy
            // Some files do not have a strand token.  While this is not allowed
            // by Genbank's documentation, we will treat this as an optional
            // token.  It is the only optional token this parser will allow.
            if (includedStrandTag)
            {
                String strandString = locusTokens.nextToken();
                StringTokenizer strandTokens = new StringTokenizer(strandString, "-");
                if (strandTokens.countTokens() > 1)
                {
                    saveSeqAnno2(GenbankFormat.STRAND_NUMBER_TAG, strandTokens.nextToken());
                }
                saveSeqAnno2(GenbankFormat.TYPE_TAG, strandTokens.nextToken());
            }

            // Circularity
            saveSeqAnno2(GenbankFormat.CIRCULAR_TAG, locusTokens.nextToken());
            // Division code
            saveSeqAnno2(GenbankFormat.DIVISION_TAG, locusTokens.nextToken());
            // Date in dd-MMM-yyyy format
            saveSeqAnno2(GenbankFormat.DATE_TAG, locusTokens.nextToken());
        }
        else
        {
            //Try using currently specified positions for tokens
            saveSeqAnno2(GenbankFormat.LOCUS_TAG, theLine.substring(12, 28));
            saveSeqAnno2(GenbankFormat.SIZE_TAG, theLine.substring(29, 40));
            saveSeqAnno2(GenbankFormat.STRAND_NUMBER_TAG, theLine.substring(44, 46));
            saveSeqAnno2(GenbankFormat.TYPE_TAG, theLine.substring(47, 53));
            saveSeqAnno2(GenbankFormat.CIRCULAR_TAG, theLine.substring(55, 63));
            saveSeqAnno2(GenbankFormat.DIVISION_TAG, theLine.substring(64, 67));
            saveSeqAnno2(GenbankFormat.DATE_TAG, theLine.substring(68, 79));
        }
    }

    /**
     * Passes the tag and the text to the listener.
     *
     * @throws ParseException Thrown when an error occurs parsing the file
     */
    private void saveSeqAnno()
	throws ParseException
    {
	if (! headerTag.equals(""))
	{ // save tag and its text
	    listener.addSequenceProperty(headerTag, headerTagText.substring(0));
	    headerTag = "";
	    headerTagText = new StringBuffer("");
	}
    }

    /**
     * Private method to process a header tag and associated value.
     *
     * @param tag The tag to add
     * @param value The value of the associated tag
     * @throws ParseException Thrown when an error occurs parsing the file
     */
    private void saveSeqAnno2(String tag, String value)
        throws ParseException
    {
        value = value.trim();	// strip whitespace
        if (value.length() > 0) {
            this.saveSeqAnno();
            headerTag = tag;
            headerTagText = new StringBuffer(value);
        }
    }

    /**
     * @return does the line contain a header tag.
     * Yes, if any of the leading TAG_LENGTH characters aren't a space
     */
    private boolean hasHeaderTag(String line)
    {
	boolean isHeaderTag = false;

	char [] lar = line.toCharArray();
        int len = Math.min(lar.length, TAG_LENGTH); // handles empty lines better
	for (int i = 0; i < len && i < lar.length; i++)
	{
            if (lar[i] != ' ')
            {
                isHeaderTag = true;
                break;
            }
        }

	return isHeaderTag;
    }

    public boolean getElideSymbols()
    {
        return elideSymbols;
    }

    public void setElideSymbols(boolean elideSymbols)
    {
        this.elideSymbols = elideSymbols;
    }
}
