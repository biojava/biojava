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
import java.util.List;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * This class performs the detailed formatting of refseq protein entries.
 * Functionality is essentially identical to GenbankFileFormer except that
 * SimpleFeatures are created intead of StrandedFeatures
 *
 * @author Greg Cox
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 */
public class ProteinRefSeqFileFormer extends GenbankFileFormer
{
    // Constructors and initialization

    protected ProteinRefSeqFileFormer()
    {
        super(System.out);
    }

    protected ProteinRefSeqFileFormer(PrintStream theStream)
    {
        super(theStream);
    }

    // Interface implementations
    public void addSymbols(Alphabet  theAlphabet,
                           Symbol [] theSymbols,
                           int       theStart,
                           int       theLength)
        throws IllegalAlphabetException
    {
        // Get newline character
        String newLine = System.getProperty("line.separator");
        this.getPrintStream().print("ORIGIN" + newLine);
        List brokenLines = this.breakSymbolArray(theAlphabet, theSymbols,
                                                 theStart, theLength);

        java.util.Iterator iterator = brokenLines.iterator();
        String leader = "     ";
        while(iterator.hasNext())
        {
            this.getPrintStream().print(leader + iterator.next() + newLine);
        }
        this.getPrintStream().println("//" + newLine);
    }

    public void startFeature(Feature.Template templ)
        throws ParseException
    {
        // CAUTION: This hasn't been tested.  Use at your own risk.

        // There are 21 spaces in the leader
        String leader = "                     ";
        int    strand = 0;

        StringBuffer theBuffer = new StringBuffer();
        formatLocationBlock(theBuffer, templ.location, strand, leader, 80);
        theBuffer.replace(5, 5 + templ.type.length(), templ.type);
        this.getPrintStream().println(theBuffer);
    }

    // Public methods

    // Protected methods

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
        } catch (BioException ex) {
            throw new BioError("Expected tokenization",ex);
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

            try
            {
                tempString.append(tokenization.tokenizeSymbol(theSymbols[i]));
            }
            catch (IllegalSymbolException e)
            {
                throw new IllegalAlphabetException(e);
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
}
