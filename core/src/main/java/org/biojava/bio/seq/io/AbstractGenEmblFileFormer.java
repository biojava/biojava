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

import java.io.IOException;
import java.io.InputStream;
import java.text.BreakIterator;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.RemoteFeature;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.BetweenLocation;
import org.biojava.bio.symbol.FuzzyLocation;
import org.biojava.bio.symbol.FuzzyPointLocation;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ClassTools;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

/**
 * <code>AbstractGenEmblFileFormer</code> contain file formatting code
 * common to both GenBank and EMBL file formats.
 *
 * @author Keith James
 * @author Greg Cox
 * @since 1.2
 */
class AbstractGenEmblFileFormer
{
    private static final String FEATURE_DATA_FILE =
        "org/biojava/bio/seq/io/FeatureQualifier.xml";

    private static final Map   FEATURE_DATA = new HashMap();
    private static final Map QUALIFIER_DATA = new HashMap();

    static
    {
        /* This loads an XML file containing information on which
         * qualifiers are valid (or even mandatory) for a particular
         * feature key. It also indicates whether the value should be
         * contained within quotes.
         */
        loadFeatureData(FEATURE_DATA_FILE, FEATURE_DATA, QUALIFIER_DATA);
    }

    /* Defines types of qualifier lines encountered: FIRST - the first
     * line of the qualifier, OVERWIDE - contains a very wide token
     * which can not be wrapped without breaking it, NOFIT - a line
     * which is too wide, but may be wrapped without breaking tokens,
     * FIT - a line which is short enough to add without wrapping it.
     */
    static final int FIRST    = 0;
    static final int OVERWIDE = 1;
    static final int NOFIT    = 2;
    static final int FIT      = 3;

    /* Defines the various types of Location available. Each is
     * represented differently within an EMBL/Genbank flatfile.
     */
    static final int            RANGE = 0;
    static final int            POINT = 1;
    static final int      FUZZY_RANGE = 2;
    static final int      FUZZY_POINT = 3;
    static final int BETWEEN_LOCATION = 4;

    // Get separator for system
    String nl = System.getProperty("line.separator");

    SymbolTokenization dnaTokenization;
    AlphabetIndex dnaIndex;

    int aCount, cCount, gCount, tCount, oCount;
    Symbol a, c , g, t;

    {
        dnaIndex = AlphabetManager.getAlphabetIndex(DNATools.getDNA());
        a = DNATools.a();
        c = DNATools.c();
        g = DNATools.g();
        t = DNATools.t();

        try
        {
            dnaTokenization = DNATools.getDNA().getTokenization("token");
        }
        catch (Exception e)
        {
            throw new BioError("Couldn't initialize tokenizer for the DNA alphabet",e);
        }
    }

    /**
     * <code>formatSequenceProperty</code> formats text into
     * EMBL/Genbank style header lines.
     *
     * @param sb a <code>StringBuffer</code>.
     * @param text a <code>String</code> to format.
     * @param leader a <code>String</code> to append to the start of
     * each line.
     * @param wrapWidth an <code>int</code> indicating the number of
     * columns per line.
     *
     * @return a <code>StringBuffer</code>.
     */
    StringBuffer formatSequenceProperty(StringBuffer sb,
                                        String       text,
                                        String       leader,
                                        int          wrapWidth)
    {
        BreakIterator boundary = BreakIterator.getLineInstance();
        boundary.setText(text);
        int start = boundary.first();
        int end = boundary.next();
        int lineLength = 0;

        sb.append(leader);
        while (end != BreakIterator.DONE)
        {
            String word = text.substring(start, end);
            lineLength = lineLength + word.length();

            if (lineLength >= wrapWidth)
            {
                sb.append(nl);
                sb.append(leader);
                lineLength = word.length();
            }

            sb.append(word);
            start = end;
            end = boundary.next();
        }

        return sb;
    }

    /**
     * <code>formatQualifierBlock</code> formats text into
     * EMBL/Genbank style qualifiers.
     *
     * @param text a <code>String</code> to format.
     * @param leader a <code>String</code> to append to the start of
     * each line.
     * @param wrapWidth an <code>int</code> indicating the number of
     * columns per line.
     *
     * @return a <code>StringBuffer</code>.
     */
    StringBuffer formatQualifierBlock(StringBuffer sb,
                                      String       text,
                                      String       leader,
                                      int          wrapWidth)
    {
        int tokenType = FIRST;
        int  position = leader.length();

        sb.append(leader);

        StringTokenizer t = new StringTokenizer(text);

    TOKEN:
        while (t.hasMoreTokens())
        {
            String s = t.nextToken();
            String separator = "";
            int tokenLen = s.length();

            // The first token has to be treated differently. It
            // always starts immediately after the '=' character
            if (! (tokenType == FIRST))
            {
                separator = " ";

                if (tokenLen + 1 > wrapWidth)
                    tokenType = OVERWIDE;
                else if (position + tokenLen + 1 > wrapWidth)
                    tokenType = NOFIT;
                else
                    tokenType = FIT;
            }

            switch (tokenType)
            {
                case FIRST:
                    // The first line always always starts immediately
                    // after the '=' character, even if it means
                    // forcing a break
                    if (! (position + tokenLen > wrapWidth))
                    {
                        sb.append(s);
                        position += s.length();
                        tokenType = FIT;
                        continue TOKEN;
                    }
                    separator = " ";

                case OVERWIDE:
                    // Force breaks in the token until the end is
                    // reached
                    for (int i = 0; i < tokenLen; i++)
                    {
                        if (position == wrapWidth)
                        {
                            sb.append(nl);
                            sb.append(leader);
                            position = leader.length();
                        }
                        sb.append(s.charAt(i));
                        position++;
                    }

                    position = tokenLen % wrapWidth;
                    break;

                case NOFIT:
                    // Token won't fit, so pass it to the next line
                    sb.append(nl);
                    sb.append(leader);
                    sb.append(s);
                    position = tokenLen + leader.length();
                    break;

                case FIT:
                    // Token fits on this line
                    sb.append(separator);
                    sb.append(s);
                    position += (tokenLen + 1);
                    break;

                default:
                    // Nothing
                    break;
            } // end switch
        } // end while

        return sb;
    }

    /**
     * <code>formatQualifier</code> creates a qualifier string and
     * adds quotes or parens to EMBL/Genbank qualifier values as
     * specified in the internally loaded XML feature table
     * description.
     *
     * @param key an <code>Object</code> (the qualifier key).
     * @param value an <code>Object</code> (the qualifier content).
     *
     * @return a <code>String</code> bounded by the correct tokens.
     */
    StringBuffer formatQualifier(StringBuffer sb, Object key, Object value)
    {
        sb.append('/');
        sb.append(key);

        // Default is to quote unknown qualifiers
        String form = "quoted";
        if (QUALIFIER_DATA.containsKey(key))
            form = (String) ((Map) QUALIFIER_DATA.get(key)).get("form");

        // This is a slight simplification. There are some types of
        // qualifier which are unquoted unless they contain
        // spaces. We all love special cases, don't we?
        if (form.equals("quoted"))
        {
            sb.append("=\"");
            sb.append(value);
            sb.append("\"");
        }
        else if (form.equals("bare"))
        {
            sb.append('=');
            sb.append(value);
        }
        else if (form.equals("paren"))
        {
            sb.append('(');
            sb.append(value);
            sb.append(')');
        }
        else if (! form.equals("empty"))
        {
            sb.append('=');
            sb.append(value);
        }

        return sb;
    }

    /**
     * <code>formatTokenBlock</code> divides up the tokens
     * representing the <code>Symbols</code> into blocks of the
     * specified length, with a single space delimeter.
     *
     * @param syms a <code>Symbol []</code> array whose tokens are to
     * be formatted.
     * @param blockSize an <code>int</code> indicating the size of
     * each block.
     *
     * @return a <code>StringBuffer</code>.
     */
    StringBuffer formatTokenBlock(StringBuffer       sb,
                                  Symbol []          syms,
                                  int                blockSize,
                                  SymbolTokenization tokenization)
        throws IllegalSymbolException
    {
        for (int i = 0; i < syms.length; i++)
        {
            sb.append(tokenization.tokenizeSymbol(syms[i]));
            if ((i + 1) % blockSize == 0)
                sb.append(' ');
        }
        return sb;
    }

    /**
     * Formats the location of a feature.  This version is required
     * when formatting remote locations, since the location field of a
     * remote feature is the projection of that feature on the
     * sequence.  When a distinction is made between 'order' and
     * 'join' this method will likely be extended for that also.
     *
     * @param theFeature The feature with the location to format
     *
     * @return String The formatted location
     */
    public String formatLocation(Feature theFeature)
    {
        String formattedLocation = null;
        StrandedFeature.Strand featureStrand = StrandedFeature.POSITIVE;

        String joinType = "join";
        Annotation ann = theFeature.getAnnotation();

        if (ann.containsProperty(Feature.PROPERTY_DATA_KEY))
        {
            Map dat = (Map) ann.getProperty(Feature.PROPERTY_DATA_KEY);

            if (dat.containsKey("JoinType"))
            {
                joinType = (String) dat.get("JoinType");
            }
        }

        if (theFeature instanceof RemoteFeature)
        {
            StringBuffer tempBuffer = new StringBuffer();
            List regionList = ((RemoteFeature)theFeature).getRegions();

            if (regionList.size() > 1)
            {
                tempBuffer.append(joinType);
                tempBuffer.append('(');
            }

            java.util.ListIterator tempIterator = regionList.listIterator();

            while (tempIterator.hasNext())
            {
                // Set up remote location
                RemoteFeature.Region tempRegion = (RemoteFeature.Region)(tempIterator.next());

                if (tempRegion.getSeqID() != null)
                {
                    tempBuffer.append(tempRegion.getSeqID());
                    tempBuffer.append(':');
                }

                // Add the actual location
                String tempLocation = null;
                if (theFeature instanceof StrandedFeature)
                {
                    featureStrand = ((StrandedFeature)theFeature).getStrand();
                }
                tempLocation = this.formatLocation(tempRegion.getLocation(), featureStrand);
                tempBuffer.append(tempLocation);

                // Only have commas between two subregions
                if (tempIterator.hasNext())
                {
                    tempBuffer.append(',');
                }
            }

            // Close out the join block if needed
            if (regionList.size() > 1)
            {
                tempBuffer.append(')');
            }
            formattedLocation = tempBuffer.substring(0);
        }
        else
        {
            if (theFeature instanceof StrandedFeature)
            {
                featureStrand = ((StrandedFeature)theFeature).getStrand();
            }

            StringBuffer tempBuffer = new StringBuffer();
            formattedLocation = this.formatLocationBlock(tempBuffer,
                                                         theFeature.getLocation(), featureStrand.getValue(), "",
                                                         Integer.MAX_VALUE, joinType).toString();
        }

        return formattedLocation;
    }

    /**
     * <code>formatLocation</code> creates an EMBL/Genbank style
     * representation of a <code>Location</code>. This is a
     * convenience method only. The version which has a
     * <code>StringBuffer</code> parameter (and returns the
     * <code>StringBuffer</code>) is preferred.  If a compound location is
     * formatted using this method, it is returned as a join-type location
     * rather than an order-type.
     *
     * @param loc a <code>Location</code> to format.
     * @param strand a <code>StrandedFeature.Strand</code>
     * indicating the <code>Location</code>'s strand.
     *
     * @return a <code>StringBuffer</code>.
     */
    public String formatLocation(Location               loc,
                                 StrandedFeature.Strand strand)
    {
        // Using arbitrary leader and wrapwidth wide enough to always
        // make one line.
        StringBuffer sb = formatLocationBlock(new StringBuffer(),
                                              loc,
                                              strand.getValue(),
                                              "",
                                              Integer.MAX_VALUE,
                                              "join");
        return sb.substring(0);
    }

    /**
     * <code>formatLocation</code> creates an EMBL/Genbank style
     * representation of a <code>Location</code>. Supported location
     * forms:
     *
     * <pre>
     *   123
     *  <123 or >123
     *  (123.567)
     *  (123.567)..789
     *   123..(567.789)
     *  (123.345)..(567.789)
     *   123..456
     *  <123..567 or 123..>567 or <123..>567
     *   123^567
     *   AL123465:(123..567)
     * </pre>
     *
     * If a compound location is formatted using this method, it is returned as
     * a join-type location rather than an order-type.
     *
     * To preserve the join/order distinction; and to format locations like
     * AL123465:(123..567), use the formatLocation(Feature) method.
     *
     * @param sb a <code>StringBuffer</code to which the location will
     * be appended.
     * @param loc a <code>Location</code> to format.
     * @param strand a <code>StrandedFeature.Strand</code>
     * indicating the <code>Location</code>'s strand.
     *
     * @return a <code>StringBuffer</code>.
     */
    public StringBuffer formatLocation(StringBuffer           sb,
                                       Location               loc,
                                       StrandedFeature.Strand strand)
    {
        // Using arbitrary leader and wrapwidth wide enough to always
        // make one line
        return formatLocationBlock(sb, loc, strand.getValue(), "", Integer.MAX_VALUE, "join");
    }

    /**
     * <code>formatLocationBlock</code> creates an EMBL/Genbank style
     * representation of a <code>Location</code> wrapped to a specific
     * width.
     *
     * If a compound location is formatted using this method, it is returned as
     * a join-type location rather than an order-type.
     *
     * To preserve the join/order distinction; and to format locations like
     * AL123465:(123..567), use the formatLocation(Feature) method.
     *
     * @param loc a <code>Location</code> to use as a template.
     * @param strand an <code>int</code> indicating the
     * <code>Location</code>'s strand.
     * @param leader a <code>String</code> to append to the start of
     * each line.
     * @param wrapWidth an <code>int</code> indicating the number of
     * columns per line.
     *
     * @return a <code>StringBuffer</code>.
     */
    StringBuffer formatLocationBlock(StringBuffer sb,
                                     Location     loc,
                                     int          strand,
                                     String       leader,
                                     int          wrapWidth)
    {
        return this.formatLocationBlock(sb, loc, strand, leader, wrapWidth, "join");
    }

    /**
     * <code>formatLocationBlock</code> creates an EMBL/Genbank style
     * representation of a <code>Location</code> wrapped to a specific
     * width.
     *
     * @param loc a <code>Location</code> to use as a template.
     * @param strand an <code>int</code> indicating the
     * <code>Location</code>'s strand.
     * @param leader a <code>String</code> to append to the start of
     * each line.
     * @param wrapWidth an <code>int</code> indicating the number of
     * columns per line.
     * @param joinType Only used if the location is a compound
     * location. It is prepended to the list of locations.
     *
     * @return a <code>StringBuffer</code>.
     */
    StringBuffer formatLocationBlock(StringBuffer sb,
                                     Location     loc,
                                     int          strand,
                                     String       leader,
                                     int          wrapWidth,
                                     String       joinType)
    {
        // Indicates how many characters have been added to the
        // current line
        int       position = leader.length();
        boolean       join = false;
        boolean complement = false;

        List locs = new ArrayList();
        for (Iterator li = loc.blockIterator(); li.hasNext();)
        {
            locs.add(li.next());
        }

        /* There are issues here about choosing various forms:
         * join(complement(...),complement(...))
         * complement(join(...,...))
         *
         * The former has the locations sorted in reverse order.
         */

        Collections.sort(locs, Location.naturalOrder);

        if (!loc.isContiguous())
        {
            join = true;
            sb.append(joinType);
            sb.append('(');
            position += 5;
        }

        if (strand == -1)
        {
            Collections.reverse(locs);
            complement = true;
        }

        int locType = 0;

        // Records the length of the String(s) added to the buffer to
        // determine whether we need to wrap the line
        int pre, post;
        int diff = 0;

        for (Iterator li = locs.iterator(); li.hasNext();)
        {
            Location thisLoc = (Location) li.next();

            pre = sb.length();

            if (PointLocation.class.isInstance(thisLoc))
                locType = POINT;
            else if (FuzzyLocation.class.isInstance(thisLoc))
                locType = FUZZY_RANGE;
            else if (FuzzyPointLocation.class.isInstance(thisLoc))
                locType = FUZZY_POINT;
            else if (BetweenLocation.class.isInstance(thisLoc))
                locType = BETWEEN_LOCATION;
            else
                locType = RANGE;

            StringBuffer lb = new StringBuffer();
            switch (locType)
            {
                case POINT:
                    PointLocation pl = (PointLocation) thisLoc;

                    sb.append(complement                                     ?
                              toComplement(formatPoint(lb, pl).substring(0)) :
                              formatPoint(lb, pl).substring(0));
                    break;

                case FUZZY_RANGE:
                    FuzzyLocation fl = (FuzzyLocation) thisLoc;

                    sb.append(complement                                          ?
                              toComplement(formatFuzzyRange(lb, fl).substring(0)) :
                              formatFuzzyRange(lb, fl).substring(0));
                    break;

                case FUZZY_POINT:
                    FuzzyPointLocation fpl = (FuzzyPointLocation) thisLoc;

                    sb.append(complement                                           ?
                              toComplement(formatFuzzyPoint(lb, fpl).substring(0)) :
                              formatFuzzyPoint(lb, fpl).substring(0));
                    break;

                case RANGE:
                    sb.append(complement                                     ?
                              toComplement(formatRange(lb, thisLoc).substring(0)) :
                              formatRange(lb, thisLoc).substring(0));
                    break;

                case BETWEEN_LOCATION:
                    BetweenLocation tempLocation = (BetweenLocation) thisLoc;
                    String formattedLocation = formatBetween(lb, tempLocation).toString();
                    if (complement)
                    {
                        formattedLocation = toComplement(formattedLocation);
                    }
                    sb.append(formattedLocation);
                    break;

                default:
                    // Maybe exception here?
                    break;
            }

            // If there is another location after this
            if ((locs.indexOf(thisLoc) + 1) < locs.size())
                sb.append(',');

            post = sb.length();

            // The number of characters just added
            diff = post - pre;

            // If we have exceeded the line length
            if ((position + diff) > wrapWidth)
            {
                // Insert a newline just prior to this location string
                sb.insert((sb.length() - diff), nl + leader);
                position = leader.length() + diff;
            }
            else
            {
                position += diff;
            }
        }

        if (join)
        {
            sb.append(')');
            // If adding the ")" has made the line too long, move the
            // last range to the next line
            if ((position + 1) > wrapWidth)
            {
                sb.insert((sb.length() - diff), nl + leader);
                position++;
                diff++;
            }
        }

        return sb;
    }

    /**
     * <code>formatFuzzyRange</code> creates an EMBL/Genbank style
     * String representation of a <code>FuzzyLocation</code>.
     *
     * @param sb a <code>StringBuffer</code> to format the location
     * into.
     * @param fl a <code>FuzzyLocation</code>.
     *
     * @return a <code>String</code> representation of the location.
     */
    private StringBuffer formatFuzzyRange(StringBuffer sb, FuzzyLocation fl)
    {
        if (! fl.hasBoundedMin())
        {
            // <123
            sb.append('<');
            sb.append(fl.getMin());
        }
        else if (fl.getOuterMin() != fl.getInnerMin())
        {
            // (123.567)
            sb.append('(');
            sb.append(fl.getOuterMin());
            sb.append('.');
            sb.append(fl.getInnerMin());
            sb.append(')');
        }
        else
        {
            // 123
            sb.append(fl.getMin());
        }

        sb.append("..");

        if (! fl.hasBoundedMax())
        {
            // >567
            sb.append('>');
            sb.append(fl.getMax());
        }
        else if (fl.getInnerMax() != fl.getOuterMax())
        {
            // (567.789)
            sb.append('(');
            sb.append(fl.getInnerMax());
            sb.append('.');
            sb.append(fl.getOuterMax());
            sb.append(')');
        }
        else
        {
            // 567
            sb.append(fl.getMax());
        }

        return sb;
    }

    /**
     * <code>formatFuzzyPoint</code> creates an EMBL/Genbank style
     * String representation of a <code>FuzzyPointLocation</code>.
     *
     * @param sb a <code>StringBuffer</code> to format the location into.
     * @param fpl a <code>FuzzyPointLocation</code>.
     *
     * @return a <code>String</code> representation of the location.
     */
    private StringBuffer formatFuzzyPoint(StringBuffer sb, FuzzyPointLocation fpl)
    {
        if (! fpl.hasBoundedMin())
        {
            // <123
            sb.append('<');
            sb.append(fpl.getMax());
        }
        else if (! fpl.hasBoundedMax())
        {
            // >567
            sb.append('>');
            sb.append(fpl.getMin());
        }
        else
        {
            // (567.789)
            sb.append('(');
            sb.append(fpl.getMin());
            sb.append('.');
            sb.append(fpl.getMax());
            sb.append(')');
        }

        return sb;
    }

    /**
     * <code>formatRange</code> creates an EMBL/Genbank style String
     * representation of a <code>RangeLocation</code>.
     *
     * @param sb a <code>StringBuffer</code> to format the location into.
     * @param rl a <code>Location</code>.
     *
     * @return a <code>String</code> representation of the location.
     */
    private StringBuffer formatRange(StringBuffer sb, Location rl)
    {
        // 123..567
        sb.append(rl.getMin());
        sb.append("..");
        sb.append(rl.getMax());

        return sb;
    }

    /**
     * <code>formatPoint</code> creates an EMBL/Genbank style String
     * representation of a <code>PointLocation</code>.
     *
     * @param sb a <code>StringBuffer</code> to format the location into.
     * @param pl a <code>PointLocation</code>.
     *
     * @return a <code>String</code> representation of the location.
     */
    private StringBuffer formatPoint(StringBuffer sb, PointLocation pl)
    {
        sb.append(Integer.toString(pl.getMin()));
        return sb;
    }

    /**
     * Formats a between location x y into x^y.
     *
     * @param sb a <code>StringBuffer</code> to format the location into.
     * @param theLocation The between location object to be formatted
     *
     * @return A string representation of the location
     */
    private StringBuffer formatBetween(StringBuffer sb, BetweenLocation theLocation)
    {
        sb.append(theLocation.getMin());
        sb.append('^');
        sb.append(theLocation.getMax());
        return sb;
    }

    /**
     * <code>toComplement</code> accepts an EMBL/Genbank style String
     * representation of a <code>Location</code> and returns the
     * complementary strand version.
     *
     * @param value a <code>String</code>.
     *
     * @return a <code>String</code> representation of the
     * complementary strand location.
     */
    private String toComplement(String value)
    {
        return "complement(" + value + ")";
    }

    /**
     * <code>loadFeatureData</code> reads data describing EMBL/Genbank
     * features and qualifiers from an XML file and populates two data
     * Maps, one for features, one for qualifiers. The file describes
     * which qualifiers are optional or mandatory for a feature type
     * and which qualifiers are written within quotes (or
     * parantheses). The DTD used by the XML file is stored in the
     * file's internal subset.
     *
     * @param featureDataFile a <code>String</code> indicating the
     * name of the file.
     * @param featureData a <code>Map</code> to populate with
     * feature data.
     * @param qualifierData a <code>Map</code> to populate with
     * qualifier data.
     */
    static void loadFeatureData(String featureDataFile,
                                Map    featureData,
                                Map    qualifierData)
    {
        try
        {
            InputStream featureDataStream  =
                ClassTools.getClassLoader(EmblFileFormer.class).getResourceAsStream(featureDataFile);
            if (featureDataStream == null)
                throw new BioError("Failed to find resource: "
                                   + featureDataFile);

            InputSource   is = new InputSource(featureDataStream);
            DocumentBuilder parser = DocumentBuilderFactory.newInstance().newDocumentBuilder();

            // Get document and then the root element
            Document doc          = parser.parse(is);
            NodeList featureNodes = doc.getDocumentElement().getChildNodes();

            // For nodes in root element (features)
            int fNodeCount = featureNodes.getLength();

            for (int i = 0; i < fNodeCount; i++)
            {
                Node featureNode = featureNodes.item(i);
                if (! (featureNode instanceof Element))
                    continue;

                Element  feature = (Element) featureNode;
                String fNodeName = feature.getNodeName();

                if (fNodeName.equals("feature"))
                {
                    String featureKey = feature.getAttribute("key");

                    NodeList qualifierNodes = feature.getChildNodes();

                    // For nodes in each feature (qualifiers)
                    int qNodeCount = qualifierNodes.getLength();
                    for (int j = 0; j < qNodeCount; j++)
                    {
                        Node qualifierNode = qualifierNodes.item(j);
                        if (! (qualifierNode instanceof Element))
                            continue;

                        Element qualifier = (Element) qualifierNode;
                        String  qNodeName = qualifier.getNodeName();

                        if (qNodeName.equals("qualifier"))
                        {
                            Map qData = new HashMap();

                            qData.put("form", qualifier.getAttribute("form"));
                            qData.put("mandatory",
                                      new Boolean(qualifier.getAttribute("mandatory")));

                            qualifierData.put(qualifier.getAttribute("name"), qData);
                        }
                    }
                    featureData.put(featureKey, qualifierData.keySet());
                }

                featureDataStream.close();
            }
        }
        catch (IOException ioe)
        {
            ioe.printStackTrace();
        }
        catch (SAXException se)
        {
            se.printStackTrace();
        }
        catch (BioError be)
        {
            be.printStackTrace();
        }
        catch (ParserConfigurationException ex)
        {
            ex.printStackTrace();
        }
    }
}
