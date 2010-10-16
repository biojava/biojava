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
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * <code>GenpeptFileFormer</code> performs the detailed formatting of
 * Genpept entries for writing to a <code>PrintStream</code>. It is based
 * on GenbankFileFormer.
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 *
 * @author Nimesh Singh
 */
class GenpeptFileFormer extends AbstractGenEmblFileFormer
    implements SeqFileFormer
{
    private PrintStream stream;

    // Main sequence formatting buffer
    private StringBuffer sq = new StringBuffer();
    // Utility formatting buffer
    private StringBuffer ub = new StringBuffer();

    // I am adding an incredible number of buffers here to hold everything.
    private StringBuffer idb = null;
    private StringBuffer acb = null;
    private StringBuffer deb = null;
    private StringBuffer svb = null;
    private StringBuffer kwb = null;
    private StringBuffer osb = null;
    private StringBuffer ocb = null;
    private StringBuffer ccb = null;
    private StringBuffer ftb = new StringBuffer();

    //locusline buffers
    private StringBuffer typeb = new StringBuffer();
    private StringBuffer strb = new StringBuffer();
    private StringBuffer sizeb = new StringBuffer();
    private StringBuffer circb = new StringBuffer();
    private StringBuffer mdatb = new StringBuffer();
    private StringBuffer divb = new StringBuffer();

    //End of my changes here

    private SymbolTokenization proteinTokenization;

    {
        try {
            proteinTokenization = ProteinTools.getTAlphabet().getTokenization("token");
        } catch (BioException ex) {
            throw new BioError("Couldn't initialize tokenizer for the PROTEIN alphabet",ex);
        }
    }

    /**
     * Creates a new <code>GenpeptFileFormer</code> using
     * <code>System.out</code> stream.
     */
    protected GenpeptFileFormer()
    {
        super();
        stream = System.out;
    }

    /**
     * Creates a new <code>GenpeptFileFormer</code> using the
     * specified stream.
     *
     * @param stream a <code>PrintStream</code>.
     */
    protected GenpeptFileFormer(PrintStream stream)
    {
        this.stream = stream;
    }

    public PrintStream getPrintStream()
    {
        return stream;
    }

    public void setPrintStream(PrintStream stream)
    {
        this.stream = stream;
    }

    public void setName(String id) throws ParseException {
        idb = new StringBuffer("LOCUS       " + id);
    }

    public void startSequence() throws ParseException { }

    public void endSequence() throws ParseException { }

    public void setURI(String uri) throws ParseException { }

    public void addSymbols(Alphabet  alpha,
                           Symbol [] syms,
                           int       start,
                           int       length)
        throws IllegalAlphabetException
    {
        try {

            // My Changes are here
            locusLineCreator(length);
            if (idb != null) {stream.println(idb); }
            if (acb != null) {stream.println(acb); }
            if (svb != null) {stream.println(svb); }
            if (deb != null) {stream.println(deb); }
            if (kwb != null) {stream.println(kwb); }
            if (osb != null) {stream.println(osb);}
            if (ocb != null) {stream.println(ocb); }
            if (ccb != null) {stream.println(ccb); }

            if (ftb.length() != 0) {
                ftb.insert(0, "FEATURES             Location/Qualifiers" + nl);
                stream.print(ftb);
            }

            sq.setLength(0);
            sq.append("ORIGIN");

            // Print sequence summary header
            stream.println(sq);

            int fullLine = length / 60;
            int partLine = length % 60;

            int lineCount = fullLine;
            if (partLine > 0)
                lineCount++;

            int lineLens [] = new int [lineCount];

            // All lines are 60, except last (if present)
            Arrays.fill(lineLens, 60);
            lineLens[lineCount - 1] = partLine;

            // Prepare line 80 characters wide, sequence is subset of this
            char [] emptyLine = new char [80];

            for (int i = 0; i < lineLens.length; i++)
            {
                sq.setLength(0);
                ub.setLength(0);

                // How long is this chunk?
                int len = lineLens[i];

                // Prep the whitespace
                Arrays.fill(emptyLine, ' ');
                sq.append(emptyLine);

                // Prepare a Symbol array same length as chunk
                Symbol [] sa = new Symbol [len];

                // Get symbols and format into blocks of tokens
                System.arraycopy(syms, start + (i * 60), sa, 0, len);

                String blocks = (formatTokenBlock(ub, sa, 10, proteinTokenization)).toString();

                sq.replace(10, blocks.length() + 10, blocks);

                // Calculate the running residue count and add to the line
                String count = Integer.toString((i * 60) + 1);
                sq.replace((9 - count.length()), 9, count);

                // Print formatted sequence line
                stream.println(sq);
            }

            // Print end of entry
            stream.println("//");
        } catch (IllegalSymbolException ex) {
            throw new IllegalAlphabetException(ex, "Protein not tokenizing");
        }
    }

    private String sequenceBufferCreator(Object key, Object value) {
        StringBuffer temp = new StringBuffer();

        if (value == null) {
            temp.append((String) key);
        }
        else if (value instanceof ArrayList) {
            Iterator iter = ((ArrayList) value).iterator();
            temp.append((String) key + " " + iter.next());
            while (iter.hasNext()) {
                temp.append(nl + "            " + iter.next());
            }
        }
        else {
            StringTokenizer valueToke = new StringTokenizer((String) value, " ");
            int fullline = 80;
            int length = 0;
            temp.append((String) key);
            if (valueToke.hasMoreTokens()) {
                String token = valueToke.nextToken();

                while (true) {
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
                        temp.append(nl + "           ");
                    }
                    else if (length <= fullline) {
                        temp.append(" " + token);
                        break;
                    }
                    else {
                        temp.append(nl);
                        temp.append("            " + token);
                        break;
                    }
                }
            }
            else {
                temp.append(" ");
            }
        }

        return temp.toString();
    }

    private StringBuffer fixLength(StringBuffer temp, int length) {
        while (temp.length() < length) {
            temp.append(" ");
        }
        return temp;
    }

    private void locusLineCreator(int size) {
        idb = fixLength(idb, 30);
        typeb = fixLength(typeb, 8);

        sizeb.insert(0, size);
        while(sizeb.length() < 12) {sizeb.insert(0, " ");}
        sizeb.append(" aa ");

        if (strb.length() > 0) {
            strb.append("-");
        }
        strb = fixLength(strb, 3);
        circb = fixLength(circb, 9);
        mdatb = fixLength(mdatb, 11);
        divb = fixLength(divb, 4);
        idb.insert(29, sizeb);
        idb.insert(44, strb);
        idb.insert(47, typeb);
        idb.insert(55, circb);
        idb.insert(64, divb);
        idb.insert(68, mdatb);
        idb.setLength(79);
        System.out.println("idb length: " + idb.length());
    }

    public void addSequenceProperty(Object key, Object value)
        throws ParseException
    {
        if (key.equals("LOCUS")) {
            idb.setLength(0);
            idb.append("LOCUS       " + (String) value);
        }
        else if (key.equals("TYPE")) {
            typeb.append(value);
        }
        else if (key.equals("DIVISION")) {
            divb.append(value);
        }
        else if (key.equals("CIRCULAR")) {
            circb.append(value);
        }
        else if (key.equals("DT") || key.equals("MDAT")) {
            if (value instanceof ArrayList) {
                mdatb.append(((ArrayList) value).get(0));
            }
            else {
                mdatb.append(value);
            }
        }
        else if (key.equals("DE") || key.equals("DEFINITION")) {
            deb = new StringBuffer(sequenceBufferCreator("DEFINITION ", value));
        }
        else if (key.equals("SV") || key.equals("VERSION")) {
            if (svb != null) {
                svb.insert(11, (String) value);
            }
            else {
                svb = new StringBuffer("VERSION     " + (String) value);
            }
        }
        else if (key.equals("GI")) {
            if (svb != null) {
                svb.append("  GI:" + (String) value);
            }
            else {
                svb = new StringBuffer("VERSION       GI:" + (String) value);
            }
        }
        else if (key.equals("KW") || key.equals("KEYWORDS")) {
            kwb = new StringBuffer(sequenceBufferCreator("KEYWORDS   ", value));
        }
        else if (key.equals("OS") || key.equals("SOURCE")) {
            osb = new StringBuffer(sequenceBufferCreator("SOURCE     ", value));
        }
        else if (key.equals("OC") || key.equals("ORGANISM")) {
            ocb = new StringBuffer(sequenceBufferCreator("  ORGANISM ", value));
        }
        else if (key.equals("CC") || key.equals("COMMENT")) {
            ccb = new StringBuffer(sequenceBufferCreator("COMMENT    ", value));
        }
        else if (key.equals(GenbankProcessor.PROPERTY_GENBANK_ACCESSIONS))
        {
            ub.setLength(0);
            ub.append("ACCESSION   ");
            List l;
            if(value instanceof List){
                l = (List)value;
            }else{
                l = new ArrayList();
                l.add(value);
            }
            for (Iterator ai = l.iterator(); ai.hasNext();)
            {
                ub.append((String) ai.next());
            }
            acb = new StringBuffer(ub.toString());
        }
    }

    //null implementation
    public void startFeature(Feature.Template templ)
        throws ParseException
    {
    }

    public void endFeature() throws ParseException { }

    //null implementation
    public void addFeatureProperty(Object key, Object value)
        throws ParseException
    {
    }
}
