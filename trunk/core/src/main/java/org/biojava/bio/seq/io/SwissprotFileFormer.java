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

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;

import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.FuzzyLocation;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.Symbol;

/**
 * Formats a sequence into Swissprot/TrEMBL format.  Modeled after
 * EmblFileFormer.
 *
 * @author Greg Cox
 * @since 1.2
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 */
public class SwissprotFileFormer extends AbstractGenEmblFileFormer
    implements SeqFileFormer
{
    // Main qualifier formatting buffer
    private StringBuffer qb = new StringBuffer();
    // Utility formatting buffer
    private StringBuffer ub = new StringBuffer();

    // Buffers for each possible sequence property line
    private StringBuffer idb = null;
    private StringBuffer acb = null;
    private StringBuffer dtb = null;
    private StringBuffer deb = null;
    private StringBuffer svb = null;
    private StringBuffer kwb = null;
    private StringBuffer osb = null;
    private StringBuffer ocb = null;
    private StringBuffer ccb = null;
    private StringBuffer ftb = new StringBuffer();

    // Static variables
    static int LOCATION_WIDTH = 6;

    // Member variables
    PrintStream mStream;

    // Constructors and initialization

    /**
     * Creates a new <code>SwissprotFileFormer</code> using
     * <code>System.out</code> stream.
     */
    protected SwissprotFileFormer()
    {
        super();
        this.mStream = System.out;
    }

    /**
     * Creates a new <code>SwissprotFileFormer</code> using the
     * specified stream.
     *
     * @param theStream a <code>PrintStream</code> object.
     */
    protected SwissprotFileFormer(PrintStream theStream)
    {
        super();
        this.mStream = theStream;
    }

    // Interface implementations
    // SeqIOListener methods

    /**
     * Start the processing of a sequence.  This method exists primarily
     * to enforce the life-cycles of SeqIOListener objects.
     */
    public void startSequence() throws ParseException
    {
    }

    /**
     * Notify the listener that processing of the sequence is complete.
     */
    public void endSequence() throws ParseException
    {
    }

    /**
     * The name is printed out as part of the identifier line.  It will be
     * replaced if an ID keyword exists in the annotations.
     *
     * @param theName the String that should be returned by getName for the sequence
     * being parsed
     */
    public void setName(String theName) throws ParseException
    {
        idb = new StringBuffer("ID   " + theName);
    }

    /**
     * Null implementation.  This object formats and prints a sequence.  The
     * URI alone cannot be printed in Swissprot format.  Therefore, it's
     * easiest to ignore it.	 
     * @param theURI the new URI of the sequence
     */
    public void setURI(String theURI) throws ParseException
    {
    }

    /**
     * Prints out the sequences properties in order.
     * Prints out the symbol array passed in in lines of 60, blocks of 10
     *
     * @param theAlphabet The alphabet of the symbol data
     * @param theSymbols An array containing symbols
     * @param theStart The start offset of valid data within the array
     * @param theLength The number of valid symbols in the array
     *
     * @throws IllegalAlphabetException if we can't cope with this
     *                                  alphabet.
     */
    public void addSymbols(Alphabet theAlphabet,
                           Symbol[] theSymbols,
                           int theStart,
                           int theLength)
        throws IllegalAlphabetException
    {

        PrintStream stream = this.getPrintStream();

        // Print out all of the sequence properties in order
        if (idb != null) {stream.println(idb); stream.println("XX");}
        if (acb != null) {stream.println(acb); stream.println("XX");}
        if (svb != null) {stream.println(svb); stream.println("XX");}
        if (dtb != null) {stream.println(dtb); stream.println("XX");}
        if (deb != null) {stream.println(deb); stream.println("XX");}
        if (kwb != null) {stream.println(kwb); stream.println("XX");}
        if (osb != null) {stream.println(osb);}
        if (ocb != null) {stream.println(ocb); stream.println("XX");}
        if (ccb != null) {stream.println(ccb); stream.println("XX");}
        if (ftb.length() != 0) {
            stream.print(ftb);
        }

        this.printOutSequenceHeaderLine(theAlphabet, theSymbols, theStart, theLength);

        List brokenLines = this.breakSymbolArray(theAlphabet, theSymbols,
                                                 theStart, theLength);

        java.util.Iterator iterator = brokenLines.iterator();
        String leader = "     ";
        while(iterator.hasNext())
        {
            stream.print(leader + iterator.next() + nl);
        }
        stream.println("//");
    }

    /**
     * Formats sequence properties into form suitable for printing to
     * file.
     *
     * @param key    The key of the sequence property
     * @param value  The value of the sequence property
     *
     * @returns      Properly formated string
     */
    private String sequenceBufferCreator(Object key, Object value) {
        StringBuffer temp = new StringBuffer();

        if (value == null) {
            temp.append((String) key);
        }
        else if (value instanceof ArrayList) {
            Iterator iter = ((ArrayList) value).iterator();
            while (iter.hasNext()) {
                temp.append((String) key + "   " + iter.next());
                if (iter.hasNext())
                    temp.append(nl);
            }
        }
        else {
            StringTokenizer valueToke = new StringTokenizer((String) value, " ");
            int fullline = 80;
            int length = 0;
            String token = valueToke.nextToken();

            while (true) {
                temp.append((String) key + "  ");
                length = (temp.length() % (fullline + 1)) + token.length() + 1;
                if (temp.length() % (fullline + 1) == 0) length = 81 + token.length();
                while (length <= fullline && valueToke.hasMoreTokens()) {
                    temp.append(" " + token);
                    token = valueToke.nextToken();
                    length = (temp.length() % (fullline + 1)) + token.length() + 1;
                    if (temp.length() % (fullline + 1) == 0) length = 81 + token.length();
                }
                if (valueToke.hasMoreTokens()) {
                    for(int i = length-token.length(); i < fullline; i++) {
                        temp.append(" ");
                    }
                    temp.append(nl);
                }
                else if (length <= fullline) {
                    temp.append(" " + token);
                    break;
                }
                else {
                    temp.append(nl);
                    temp.append((String) key + "   " + token);
                    break;
                }
            }
        }

        return temp.toString();
    }

    /**
     * Notify the listener of a sequence-wide property.  This might
     * be stored as an entry in the sequence's annotation bundle.
     * Checks for possible known properties to be shown in the file.
     *
     * @param key Key the property will be stored under
     * @param value Value stored under the key
     */
    public void addSequenceProperty(Object key, Object value) throws ParseException
    {
        if (key.equals("ID")) {
            idb.setLength(0);
            idb.append("ID   " + (String) value);
        }
        else if (key.equals("DT") || key.equals("MDAT")) {
            dtb = new StringBuffer(sequenceBufferCreator("DT", value));
        }
        else if (key.equals("DE") || key.equals("DEFINITION")) {
            deb = new StringBuffer(sequenceBufferCreator("DE", value));
        }
        else if (key.equals("SV") || key.equals("VERSION")) {
            svb = new StringBuffer(sequenceBufferCreator("SV", value));
        }
        else if (key.equals("KW") || key.equals("KEYWORDS")) {
            kwb = new StringBuffer(sequenceBufferCreator("KW", value));
        }
        else if (key.equals("OS") || key.equals("SOURCE")) {
            osb = new StringBuffer(sequenceBufferCreator("OS", value));
        }
        else if (key.equals("OC") || key.equals("ORGANISM")) {
            ocb = new StringBuffer(sequenceBufferCreator("OC", value));
        }
        else if (key.equals("CC") || key.equals("COMMENT")) {
            ccb = new StringBuffer(sequenceBufferCreator("CC", value));
        }
        else if (key.equals(SwissprotProcessor.PROPERTY_SWISSPROT_ACCESSIONS))
        {
            acb = new StringBuffer();
            acb.append("AC   ");
            for (Iterator ai = ((List) value).iterator(); ai.hasNext();)
            {
                acb.append((String) ai.next());
                acb.append(";");
            }
        }
    }

    /**
     * Null implementation.
     *
     * @param templ The template for this new feature object
     */
    public void startFeature(Feature.Template templ) throws ParseException
    {
        // There are 19 spaces in the leader
        String leader = "FT                   ";

        ub.setLength(0);
        ub.append(leader);

        StringBuffer lb = formatLocation(ub, templ.location);

        lb.replace(5, 5 + templ.type.length(), templ.type);
        ftb.append(lb + nl);
    }

    /**
     * Null implementation.
     */
    public void endFeature() throws ParseException
    {
    }

    /**
     * Null implementation
     *
     * @param key Key the property will be stored under
     * @param value Value stored under the key
     */

    public void addFeatureProperty(Object key, Object value) throws ParseException
    {
        // There are 19 spaces in the leader
        String leader = "FT                   ";

        // Don't print internal data structures
        if (key.equals(Feature.PROPERTY_DATA_KEY))
            return;

        // The value may be a collection if several qualifiers of the
        // same type are present in a feature
        if (Collection.class.isInstance(value))
        {
            for (Iterator vi = ((Collection) value).iterator(); vi.hasNext();)
            {
                qb.setLength(0);
                ub.setLength(0);
                StringBuffer fb = formatQualifierBlock(qb,
                                                       formatQualifier(ub, key, vi.next()).toString(),
                                                       leader,
                                                       80);
                ftb.append(fb + nl);
            }
        }
        else
        {
            qb.setLength(0);
            ub.setLength(0);
            StringBuffer fb = formatQualifierBlock(qb,
                                                   formatQualifier(ub, key, value).toString(),
                                                   leader,
                                                   80);
            ftb.append(fb + nl);
        }
    }

    // SeqFileFormer methods
    /**
     * <code>getPrintStream</code> returns the
     * <code>PrintStream</code> to which an instance of SwissprotFileFormer
     * will write the formatted data. The default is System.out
     *
     * @return the <code>PrintStream</code> which will be written to.
     */
    public PrintStream getPrintStream()
    {
        return(this.mStream);
    }

    /**
     * <code>setPrintStream</code> informs an instance which
     * <code>PrintStream</code> to use.
     *
     * @param theStream a <code>PrintStream</code> to write to.
     */
    public void setPrintStream(PrintStream theStream)
    {
        this.mStream = theStream;
    }

    /**
     * <code>formatLocation</code> creates a String representation of
     * a <code>Location</code>. Strand information is ignored, as Swissprot
     * files represent proteins. An alternative form of this function does not
     * take a Strand; that form is available only on SwissprotFileFormer; it
     * is not part of the SeqFileFormer interface.
     *
     * @param theBuffer a <code>StringBuffer</code> to append the location
     * to.
     * @param theLocation a <code>Location</code> to format.
     * @param theStrand a <code>StrandedFeature.Strand</code> indicating nothing
     * of relevance
     *
     * @return a <code>StringBuffer</code> with the location appended.
     */
    public StringBuffer formatLocation(StringBuffer theBuffer,
                                       Location theLocation,
                                       StrandedFeature.Strand theStrand)
    {
        return(this.formatLocation(theBuffer, theLocation));
    }

    /**
     * Creates a string representation of the location of a feature
     *
     * @param theFeature The feature with the location to format
     * @return String The formatted location
     */
    public String formatLocation(Feature theFeature)
    {
        StringBuffer toReturn = this.formatLocation(new StringBuffer(), theFeature.getLocation());
        return toReturn.toString();
    }

    // Public methods
    /**
     * <code>formatLocation</code> creates a String representation of
     * a <code>Location</code>. The stringbuffer returned represents columns
     * 15-27 of the Swissprot feature table entry. An alternative form of this
     * function takes a Strand; that form is part of the SeqFileFormer
     * interface.
     *
     * @param theBuffer a <code>StringBuffer</code> to append the location
     * to.
     * @param theLocation a <code>Location</code> to format.
     *
     * @return a <code>StringBuffer</code> with the location appended.
     */
    public StringBuffer formatLocation(StringBuffer theBuffer,
                                       Location theLocation)
        {
            // Five Location cases, each treated seperately:
            //   Point Location: "     5      5"
            //   Range Location: "     5     10"
            //   Fuzzy Location: "    <5     10"
            //   Fuzzy Location: "     ?     10"
            //   Fuzzy Location: "   ?24     35" (Not in the current
            //       specification, but used anyways
            StringBuffer startPoint = new StringBuffer(LOCATION_WIDTH);
            StringBuffer endPoint   = new StringBuffer(LOCATION_WIDTH);
            if((theLocation instanceof PointLocation) ||
               (theLocation instanceof RangeLocation))
            {
                //   Point Location: "     5      5"
                //   Range Location: "     5     10"
                startPoint = formatPoint(theLocation.getMin(), theLocation.getMin(), false);
                endPoint = formatPoint(theLocation.getMax(), theLocation.getMax(), false);
            }
            else if(theLocation instanceof FuzzyLocation)
            {
                // Handle all fuzzy location types through the magic of delegation.
                // If you pass things around long enough, someone's bound to do it
                // for you
                FuzzyLocation tempLocation = (FuzzyLocation)theLocation;
                //System.out.println("OuterMin: " + tempLocation.getOuterMin());
                //System.out.println("InnerMin: " + tempLocation.getInnerMin());
                //System.out.println("InnerMax: " + tempLocation.getInnerMax());
                //System.out.println("OuterMax: " + tempLocation.getOuterMax());
                startPoint = this.formatPoint(tempLocation.getOuterMin(),
                                              tempLocation.getInnerMin(), tempLocation.isMinFuzzy());
                endPoint = this.formatPoint(tempLocation.getInnerMax(),
                                            tempLocation.getOuterMax(), tempLocation.isMaxFuzzy());
            }

            return new StringBuffer(startPoint.toString() + " " + endPoint.toString());
        }

    // Protected methods
    /**
     * Prints out sequence header with only length data.
     *
     * @param theAlphabet The alphabet of the symbol data
     * @param theSymbols An array containing symbols
     * @param theStart The start offset of valid data within the array
     * @param theLength The number of valid symbols in the array
     *
     * @throws IllegalAlphabetException if we can't cope with this
     *                                  alphabet.
     */
    protected void printOutSequenceHeaderLine(Alphabet theAlphabet,
                                              Symbol[] theSymbols,
                                              int theStart,
                                              int theLength)
        throws IllegalAlphabetException
    {
        this.getPrintStream().println("SQ   SEQUENCE   " + theLength + " AA;   ");
    }

    /**
     * Converts the symbol list passed in into an array of strings.  The
     * strings will be blocks of ten, with six blocks on a line.
     *
     * @param theAlphabet The alphabet of the symbol data
     * @param theSymbols An array containing symbols
     * @param theStart The start offset of valid data within the array
     * @param theLength The number of valid symbols in the array
     * @return The symbol list passed in broken into blocks of ten
     * characters, six to a string.
     *
     * @throws IllegalAlphabetException if we can't cope with this
     *                                  alphabet.
     */
    protected List breakSymbolArray(Alphabet theAlphabet,
                                    Symbol[] theSymbols,
                                    int theStart,
                                    int theLength)
        throws IllegalAlphabetException
    {
        List returnList = new ArrayList(theLength / 60 + 1);
        int blockCount = 0;
        int blockIndex = 0;
        StringBuffer tempString = new StringBuffer();
        SymbolTokenization tokenization;
        try {
            tokenization = theAlphabet.getTokenization("token");
        } catch (Exception ex) {
            throw new IllegalAlphabetException(ex, "Couldn't get tokenization for this alphabet");
        }
        for(int i = theStart; i < theStart + theLength; i++)
        {
            try
            {
                theAlphabet.validate(theSymbols[i]);
            }
            catch (IllegalSymbolException e)
            {
                throw new IllegalAlphabetException(e);
            }

            // Every six completed blocks, put on the stack to return
            if(blockIndex == 10)
            {
                tempString.append(' ');
                blockIndex = 0;
                blockCount++;
            }

            if(blockCount == 6)
            {
                returnList.add(tempString.substring(0));
                tempString.setLength(0);
                blockCount = 0;
                blockIndex = 0;
            }
            try {
                tempString.append(tokenization.tokenizeSymbol(theSymbols[i]));
            } catch (IllegalSymbolException ex) {
                throw new IllegalAlphabetException(ex, "Couldn't tokenize symbols");
            }
            blockIndex++;
        }

        // Add the last line on
        if(tempString.length() != 0)
        {
            returnList.add(tempString.substring(0));
        }
        return returnList;
    }

    /**
     * Simple method that adds spaces onto the buffer passed in.  This method
     * exists to refactor some code used in location formatting.  It isn't
     * intended to be generally used.
     *
     * @param theBuffer Buffer to append whitespace to.
     * @param theLength Ammount of whitespace to append.
     */
    protected void fillBuffer(StringBuffer theBuffer, int theLength)
    {
        for(int i = 0; i < theLength; i++)
        {
            theBuffer.append(' ');
        }
    }

    /**
     * Formats the points from fuzzy locations.  This is called easily with
     * this.formatPoint(FuzzyLocation.getInnerMax(), FuzzyLocation.getOuterMax(), FuzzyLocation.isFuzzyMax())
     *
     * @param theMaxIndex Inner index of the fuzzy point
     * @param theMinIndex Outer index of the fuzzy point
     * @param isFuzzy Indicates if this point is fuzzy
     */
    protected StringBuffer formatPoint(int theMinIndex, int theMaxIndex, boolean isFuzzy)
    {
        StringBuffer bufferToReturn = new StringBuffer(LOCATION_WIDTH);
        if(isFuzzy == false)
        {
            String tempString = Integer.toString(theMinIndex);
            int offset = LOCATION_WIDTH - tempString.length();
            this.fillBuffer(bufferToReturn, offset);
            bufferToReturn.append(tempString);
        }
        else
        {
            // MIN_VALUE to MAX_VALUE is the ? location regardless of which end is which
            if((theMinIndex == Integer.MIN_VALUE) && (theMaxIndex == Integer.MAX_VALUE))
            {
                int offset = LOCATION_WIDTH - 1;
                this.fillBuffer(bufferToReturn, offset);
                bufferToReturn.append('?');
            }
            // If the outer index is MIN_VALUE, that's <n
            else if(theMinIndex == Integer.MIN_VALUE)
            {
                String tempString = Integer.toString(theMaxIndex);
                int offset = LOCATION_WIDTH - tempString.length() - 1;
                this.fillBuffer(bufferToReturn, offset);
                bufferToReturn.append('<');
                bufferToReturn.append(tempString);
            }
            // If the outer index is MAX_VALUE, that's >n
            else if(theMaxIndex == Integer.MAX_VALUE)
            {
                String tempString = Integer.toString(theMinIndex);
                int offset = LOCATION_WIDTH - tempString.length() - 1;
                this.fillBuffer(bufferToReturn, offset);
                bufferToReturn.append('>');
                bufferToReturn.append(tempString);
            }
            // The only swissprot location left is ?nn
            else if(theMinIndex == theMaxIndex)
            {
                String tempString = Integer.toString(theMinIndex);
                int offset = LOCATION_WIDTH - tempString.length() - 1;
                this.fillBuffer(bufferToReturn, offset);
                bufferToReturn.append('?');
                bufferToReturn.append(tempString);
            }
            else
            {
                // The location cannot be formatted in Swissprot format
                // Revisit
                System.out.println("Error in formatPoint");
                System.out.println("\tInner: " + theMinIndex);
                System.out.println("\tOuter: " + theMaxIndex);
                System.out.println("\tFuzzy: " + isFuzzy);
            }
        }
        return bufferToReturn;
    }

    // Private methods
}
