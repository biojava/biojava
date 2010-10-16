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
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.taxa.EbiFormat;
import org.biojava.bio.taxa.Taxon;

/**
 * <p><code>EmblFileFormer</code> performs the detailed formatting of
 * EMBL entries for writing to a <code>PrintStream</code>. Currently
 * the formatting of the header is not correct. This really needs to
 * be addressed in the parser which is merging fields which should
 * remain separate.</p>
 *
 * <p>The event generator used to feed events to this class should
 * enforce ordering of those events. This class will stream data
 * directly to the <code>PrintStream</code></p>.
 *
 * <p>This implementation requires that all the symbols be added in
 * one block as is does not buffer the tokenized symbols between
 * calls.</p>
 *
 * @author Keith James
 * @author Len Trigg (Taxon output)
 * @author Lorna Morris
 * @since 1.2
 * @deprecated Use org.biojavax.bio.seq.io framework instead
 */
public class EmblFileFormer extends AbstractGenEmblFileFormer
    implements SeqFileFormer
{
    // Tags which are special cases, not having "XX" after them
    private static List NON_SEPARATED_TAGS = new ArrayList();

    static
    {
        NON_SEPARATED_TAGS.add(EmblLikeFormat.SOURCE_TAG);
        NON_SEPARATED_TAGS.add(EmblLikeFormat.REFERENCE_TAG);
        NON_SEPARATED_TAGS.add(EmblLikeFormat.COORDINATE_TAG);
        NON_SEPARATED_TAGS.add(EmblLikeFormat.REF_ACCESSION_TAG);
        NON_SEPARATED_TAGS.add(EmblLikeFormat.AUTHORS_TAG);
        NON_SEPARATED_TAGS.add(EmblLikeFormat.TITLE_TAG);
        NON_SEPARATED_TAGS.add(EmblLikeFormat.FEATURE_TAG);
        NON_SEPARATED_TAGS.add(EmblLikeFormat.JOURNAL_TAG);//Lorna: added
        NON_SEPARATED_TAGS.add(EmblLikeFormat.REF_XREF_TAG);//RichardH: added
        NON_SEPARATED_TAGS.add(EmblLikeFormat.SEPARATOR_TAG);//Lorna: added
    }

    // 19 spaces
    private static String FT_LEADER =
        EmblLikeFormat.FEATURE_TABLE_TAG + "                   ";

    // 3 spaces
    private static String SQ_LEADER = "   ";

    // 80 spaces
    private static String EMPTY_LINE =
        "                                        " +
        "                                        ";

    private PrintStream stream;

    private String accLine;

    /**
     * Creates a new <code>EmblFileFormer</code> using
     * <code>System.out</code> stream.
     */
    protected EmblFileFormer()
    {
        this(System.out);
    }

    /**
     * Creates a new <code>EmblFileFormer</code> using the specified
     * stream.
     *
     * @param stream a <code>PrintStream</code>.
     */
    protected EmblFileFormer(PrintStream stream)
    {
        super();
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

    public void setName(String id) throws ParseException
    {
    }

    public void startSequence() throws ParseException
    {
       aCount = 0;
       cCount = 0;
       gCount = 0;
       tCount = 0;
       oCount = 0;
    }

    public void endSequence() throws ParseException
    {
        stream.println(EmblLikeFormat.END_SEQUENCE_TAG);
    }

    public void setURI(String uri) throws ParseException { }

    public void addSymbols(Alphabet  alpha,
                           Symbol [] syms,
                           int       start,
                           int       length)
        throws IllegalAlphabetException
    {
        try
        {
            int end = start + length - 1;

            for (int i = start; i <= end; i++)
            {
                Symbol sym = syms[i];

                if (sym == a)
                    aCount++;
                else if (sym == c)
                    cCount++;
                else if (sym == g)
                    gCount++;
                else if (sym == t)
                    tCount++;
                else
                    oCount++;
            }

            StringBuffer sb = new StringBuffer(EmblLikeFormat.SEPARATOR_TAG);
            sb.append(nl);
            sb.append("SQ   Sequence ");
            sb.append(length + " BP; ");
            sb.append(aCount + " A; ");
            sb.append(cCount + " C; ");
            sb.append(gCount + " G; ");
            sb.append(tCount + " T; ");
            sb.append(oCount + " other;");

            // Print sequence summary header
            stream.println(sb);

            int fullLine = length / 60;
            int partLine = length % 60;

            int lineCount = fullLine;
            if (partLine > 0)
                lineCount++;

            int lineLens [] = new int [lineCount];

            // All lines are 60, except last (if present)
            Arrays.fill(lineLens, 60);

            if (partLine > 0)
                lineLens[lineCount - 1] = partLine;

            for (int i = 0; i < lineLens.length; i++)
            {
                // Prep the whitespace
                StringBuffer sq = new StringBuffer(EMPTY_LINE);

                // How long is this chunk?
                int len = lineLens[i];
                // Prepare a Symbol array same length as chunk
                Symbol [] sa = new Symbol [len];

                // Get symbols and format into blocks of tokens
                System.arraycopy(syms, start + (i * 60), sa, 0, len);

                sb = new StringBuffer();

                String blocks = (formatTokenBlock(sb, sa, 10,
                         alpha.getTokenization("token"))).toString();

                sq.replace(5, blocks.length() + 5, blocks);

                // Calculate the running residue count and add to the line
                String count = Integer.toString((i * 60) + len);
                sq.replace((80 - count.length()), 80, count);

                // Print formatted sequence line
                stream.println(sq);
            }
        }
        catch (BioException ex)
        {
            throw new IllegalAlphabetException(ex, "Alphabet not tokenizing");
        }
    }

        public void addSequenceProperty(Object key, Object value)
        throws ParseException
    {
        StringBuffer sb = new StringBuffer();

        // Ignore separators if they are sent to us. The parser should
        // be ignoring these really (lorna: I've changed this so they are ignored in SeqIOEventEmitter)
        //if (key.equals(EmblLikeFormat.SEPARATOR_TAG))
            //return;

        String tag = key.toString();
        String leader = tag + SQ_LEADER;
        String line = "";
        int wrapWidth = 85 - leader.length();

        // Special case: accession number
        if (key.equals(EmblProcessor.PROPERTY_EMBL_ACCESSIONS))
        {
            accLine = buildPropertyLine((Collection) value, ";", true);
            return;
        }
        else if (key.equals(EmblLikeFormat.ACCESSION_TAG))
        {
            line = accLine;
        } else if (key.equals(OrganismParser.PROPERTY_ORGANISM)) {
            Taxon taxon = (Taxon) value;
            addSequenceProperty(EmblLikeFormat.SOURCE_TAG, taxon);
            addSequenceProperty(EmblLikeFormat.ORGANISM_TAG, taxon.getParent());
            addSequenceProperty(EmblLikeFormat.ORGANISM_XREF_TAG, taxon);
            return;
        }
        if (value instanceof String)
        {
            line = (String) value;
        }
        else if (value instanceof Collection)
        {
            // Special case: date lines
            if (key.equals(EmblLikeFormat.DATE_TAG))
            {
                line = buildPropertyLine((Collection) value, nl + leader, false);
                wrapWidth = Integer.MAX_VALUE;
            }
            //lorna :added 21.08.03, DR lines are another special case. Each one goes onto a separate line.
            else if (key.equals(EmblLikeFormat.DR_TAG))
            {
                line = buildPropertyLine((Collection) value, nl + leader, false);
                wrapWidth = Integer.MAX_VALUE;
            }
            else if (key.equals(EmblLikeFormat.AUTHORS_TAG))
            {
                line = buildPropertyLine((Collection) value, nl + leader, false); //lorna: add space here?
                wrapWidth = Integer.MAX_VALUE;
            }
            else if (key.equals(EmblLikeFormat.REF_ACCESSION_TAG))
            {
                line = buildPropertyLine((Collection) value, nl + leader, false);
                wrapWidth = Integer.MAX_VALUE;
            }
            else
            {
                line = buildPropertyLine((Collection) value, " ", false);
            }
        } else if (value instanceof Taxon) {
            if (key.equals(EmblLikeFormat.ORGANISM_TAG)) {
                line = EbiFormat.getInstance().serialize((Taxon) value);
            } else if (key.equals(EmblLikeFormat.SOURCE_TAG)) {
                line = EbiFormat.getInstance().serializeSource((Taxon) value);
            } else if (key.equals(EmblLikeFormat.ORGANISM_XREF_TAG)) {
                line = EbiFormat.getInstance().serializeXRef((Taxon) value);
            }
        }

        if (line.length() == 0)
        {
            stream.println(tag);
        }
        else
        {
            sb = formatSequenceProperty(sb, line, leader, wrapWidth);
            stream.println(sb);
        }
        // Special case: those which don't get separated
        if (! NON_SEPARATED_TAGS.contains(key))
            stream.println(EmblLikeFormat.SEPARATOR_TAG);
        // Special case: feature header
        if (key.equals(EmblLikeFormat.FEATURE_TAG))
            stream.println(EmblLikeFormat.FEATURE_TAG);
    }


    public void startFeature(Feature.Template templ)
        throws ParseException
    {
        int strand = 0;

        if (templ instanceof StrandedFeature.Template)
            strand = ((StrandedFeature.Template) templ).strand.getValue();

        StringBuffer sb = new StringBuffer(FT_LEADER);
        sb = formatLocationBlock(sb, templ.location, strand, FT_LEADER, 80);
        sb.replace(5, 5 + templ.type.length(), templ.type);
        stream.println(sb);
    }

    public void endFeature() throws ParseException { }

    public void addFeatureProperty(Object key, Object value)
    {
        // Don't print internal data structures
        if (key.equals(Feature.PROPERTY_DATA_KEY))
            return;

        StringBuffer fb;
        StringBuffer sb;

        // The value may be a collection if several qualifiers of the
        // same type are present in a feature
        if (value instanceof Collection)
        {
            for (Iterator vi = ((Collection) value).iterator(); vi.hasNext();)
            {
                fb = new StringBuffer();
                sb = new StringBuffer();

                fb = formatQualifierBlock(fb,
                                          formatQualifier(sb, key, vi.next()).substring(0),
                                          FT_LEADER,
                                          80);
                stream.println(fb);
            }
        }
        else
        {
            fb = new StringBuffer();
            sb = new StringBuffer();

            fb = formatQualifierBlock(fb,
                                      formatQualifier(sb, key, value).substring(0),
                                      FT_LEADER,
                                      80);
            stream.println(fb);
        }
    }

    private String buildPropertyLine(Collection property,
                                     String separator,
                                     boolean terminate)
    {
        StringBuffer sb = new StringBuffer();

        for (Iterator pi = property.iterator(); pi.hasNext();)
        {
            sb.append(pi.next().toString());
            sb.append(separator);
        }

        if (terminate)
        {
            return sb.substring(0);
        }
        else
        {
            return sb.substring(0, sb.length() - separator.length());
        }
    }
}
