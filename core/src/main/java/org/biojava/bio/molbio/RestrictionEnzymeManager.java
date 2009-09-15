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

package org.biojava.bio.molbio;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.ResourceBundle;
import java.util.Set;
import java.util.regex.Pattern;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.program.tagvalue.ChangeTable;
import org.biojava.bio.program.tagvalue.LineSplitParser;
import org.biojava.bio.program.tagvalue.Parser;
import org.biojava.bio.program.tagvalue.RegexSplitter;
import org.biojava.bio.program.tagvalue.TagDropper;
import org.biojava.bio.program.tagvalue.TagValueContext;
import org.biojava.bio.program.tagvalue.TagValueListener;
import org.biojava.bio.program.tagvalue.TagValueParser;
import org.biojava.bio.program.tagvalue.ValueChanger;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ParserException;
import org.biojava.utils.SmallSet;

/**
 * <p><code>RestrictionEnzymeManager</code> manages collections of
 * static <code>RestrictionEnzyme</code> instances. A properties file
 * should be placed in the CLASSPATH containing a key
 * "rebase.data.file" and a corresponding value of a REBASE file
 * (standard REBASE format #31 conventionally named withrefm.### where
 * ### is the version number). This file will be loaded by the
 * <code>RestrictionEnzymeManager</code> <code>ClassLoader</code>. The
 * properties are loaded as a <code>ResourceBundle</code>, so the file
 * should be named "RestrictionEnzymeManager.properties".</p>
 * <p>Since 1.5, a format #31 REBASE file can be loaded at anytime
 * using the method <code>loadEnzymeFile</code> and optionally filtered
 * for commercially available enzymes.</p>
 *  
 * @author Keith James
 * @author George Waldon
 * @since 1.3

 */
public final class RestrictionEnzymeManager
{
     /**
     * <code>REBASE_DATA_KEY</code> the ResourceBundle key which
     * specifies the location of the REBASE flat file.
     */
    public static final String REBASE_DATA_KEY = "rebase.data.file";

    /**
     * <code>REBASE_TAG_NAME</code> the REBASE tag containing the
     * enzyme name.
     */
    public static final String REBASE_TAG_NAME = "<1>";

    /**
     * <code>REBASE_TAG_ISZR</code> the REBASE tag containing the
     * enzyme isoschizomers.
     */
    public static final String REBASE_TAG_ISZR = "<2>";

    /**
     * <code>REBASE_TAG_SITE</code> the REBASE tag containing the
     * enzyme site.
     */
    public static final String REBASE_TAG_SITE = "<3>";

    /**
     * <code>REBASE_TAG_METH</code> the REBASE tag containing the
     * methylation site.
     */
    public static final String REBASE_TAG_METH = "<4>";

    /**
     * <code>REBASE_TAG_ORGN</code> the REBASE tag containing the
     * organism.
     */
    public static final String REBASE_TAG_ORGN = "<5>";

    /**
     * <code>REBASE_TAG_SRCE</code> the REBASE tag containing the
     * source.
     */
    public static final String REBASE_TAG_SRCE = "<6>";

    /**
     * <code>REBASE_TAG_COMM</code> the REBASE tag containing the
     * commercial suppliers.
     */
    public static final String REBASE_TAG_COMM = "<7>";

    /**
     * <code>REBASE_TAG_REFS</code> the REBASE tag containing the
     * references.
     */
    public static final String REBASE_TAG_REFS = "<8>";
    
    
    private static boolean loadCommercialOnly = false;

    private static ResourceBundle bundle =
        ResourceBundle.getBundle(RestrictionEnzymeManager.class.getName());

    static
    {
        String rebaseDataFileName = bundle.getString(REBASE_DATA_KEY);
        InputStream is = RestrictionEnzymeManager.class.getResourceAsStream(rebaseDataFileName);
        loadData(is);
    }

    private static Map nameToSite;
    private static Map nameToEnzyme;
    private static Map nameToIsoschizomers;
    private static Map sizeToCutters;
    private static Map enzymeToPattern;
    private static Map enzymeToAnnotation;
    private static Map enzymeToSuppliers;

    /**
     * <code>RestrictionEnzymeManager</code> is a static utility
     * method class and no instances should be created.
     */
    private RestrictionEnzymeManager() { }
    
    /**
     * <code>loadEnzymeFile</code> loads a new REBASE file (or any file using
     * REBASE format #31).
     *
     * @param is an InputStream over the file to load.
     * @param commercialOnly indicates whether or not only commercially available 
     * enzymes are loaded.
     *
     * @since 1.5
     */
    public static synchronized void loadEnzymeFile(InputStream is, boolean commercialOnly) {
        loadCommercialOnly = commercialOnly;
        loadData(is);
    }

    /**
     * <code>getAllEnzymes</code> returns an unmodifable set of all
     * available enzymes.
     *
     * @return a <code>Set</code> of <code>RestrictionEnzyme</code>s.
     */
    public static Set getAllEnzymes()
    {
        return Collections.unmodifiableSet(enzymeToPattern.keySet());
    }

    /**
     * <code>getEnzyme</code> returns an enzyme by name.
     *
     * @param name a <code>String</code> such as EcoRI, case
     * sensitive.
     *
     * @return a <code>RestrictionEnzyme</code>.
     */
    public static RestrictionEnzyme getEnzyme(String name)
    {
        if (! nameToEnzyme.containsKey(name))
            throw new IllegalArgumentException("Unknown RestrictionEnzyme name '"
                                               + name
                                               + "'");

        return (RestrictionEnzyme) nameToEnzyme.get(name);
    }

    /**
     * <code>getIsoschizomers</code> returns an unmodifable set of the
     * isoschizomers of this enzyme.
     *
     * @param name a <code>String</code> such as EcoRI, case
     * sensitive.
     *
     * @return a <code>Set</code> of <code>RestrictionEnzyme</code>s.
     */
    public static Set getIsoschizomers(String name)
    {
        if (! nameToIsoschizomers.containsKey(name))
            throw new IllegalArgumentException("Unknown RestrictionEnzyme name '"
                                               + name
                                               + "'");
        Set result = (Set) nameToIsoschizomers.get(name);
        if(result.contains(null))
            return Collections.EMPTY_SET;
        return Collections.unmodifiableSet(result);
    }

    /**
     * <code>getRecognitionSequence</code> returns a string that describes
     * the recognition site of this enzyme. It corresponds to the field <3>
     * of the REBASE file.
     *
     * @param name a <code>String</code> such as EcoRI, case
     * sensitive.
     * @return a <code>String</code> describing the recognition sequence, 
     * e.g. "G^AATTC" for EcoRI.
     * @since 1.5
     */
    public static String getRecognitionSequence(String name)
    {
        if (! nameToSite.containsKey(name))
            throw new IllegalArgumentException("Unknown RestrictionEnzyme name '"
                                               + name
                                               + "'");
        return (String) nameToSite.get(name);
    }

    /**
     * <code>getNCutters</code> returns an unmodifable set of all
     * enzymes with a cut site of size n.
     *
     * @param n an <code>int</code> cut site size.
     *
     * @return a <code>Set</code> of <code>RestrictionEnzyme</code>s.
     */
    public static Set getNCutters(int n)
    {
        Integer size = new Integer(n);
        if (! sizeToCutters.containsKey(size))
            return Collections.EMPTY_SET;

        return Collections.unmodifiableSet((Set) sizeToCutters.get(size));
    }

    /**
     * <code>getPatterns</code> returns two <code>Pattern</code>
     * objects for an enzyme, one matches the forward strand and one
     * the reverse. This enables searching of both strands of a
     * sequence without reverse-complementing it. As
     * <code>Pattern</code> objects are thread-safe these may be used
     * for all searches.
     *
     * @param enzyme a <code>RestrictionEnzyme</code>.
     *
     * @return a <code>Pattern []</code> array with the forward strand
     * <code>Pattern</code> at index 0 and the reverse at index 1.
     */
    public static Pattern [] getPatterns(RestrictionEnzyme enzyme)
    {
        if (! enzymeToPattern.containsKey(enzyme))
            throw new IllegalArgumentException("RestrictionEnzyme '"
                                               + enzyme.getName()
                                               + "' is not registered. No precompiled Pattern is available");

        return (Pattern []) enzymeToPattern.get(enzyme);
    }

    /**
     * <code>getAnnotation</code> returns an immutable, static
     * annotation describing the enzyme. This is suitable for adding
     * to <code>Feature</code>s which represent restriction sites. The
     * annotation produced currently contains one key "dbxref" in line
     * with the GenBank/EMBL qualifier for the "misc_binding" feature
     * key. The key has a corresponding value "REBASE:&lt;enzyme
     * name&gt;".
     *
     * @param enzyme a <code>RestrictionEnzyme</code>.
     *
     * @return an <code>Annotation</code>.
     */
    public static Annotation getAnnotation(RestrictionEnzyme enzyme)
    {
        if (! enzymeToAnnotation.containsKey(enzyme))
            throw new IllegalArgumentException("RestrictionEnzyme '"
                                               + enzyme.getName()
                                               + "' is not registered. No Annotation is available");

        return (Annotation) enzymeToAnnotation.get(enzyme);
    }


    /**
     * <code>getSuppliers</code> returns a string describing the suppliers
     * of this enzyme according to REBASE encoding for commercial sources 
     * or an empty String if the enzyme is not commecially available.
     *
     * <P>REBASE #31 version 604 code: </P>
     * <P>A GE Healthcare (8/05) <BR>
     * B Invitrogen Corporation(8/05)<BR>
     * C Minotech Biotechnology (9/05)<BR>
     * E Stratagene (9/05)<BR>
     * F Fermentas International Inc. (2/06)<BR>
     * G Qbiogene (9/05)<BR>
     * H American Allied Biochemical, Inc. (9/05)<BR>
     * I SibEnzyme Ltd. (2/06)<BR>
     * J Nippon Gene Co., Ltd. (8/05)<BR>
     * K Takara Bio Inc. (9/05)<BR>
     * M Roche Applied Science (8/05)<BR>
     * N New England Biolabs (2/06)<BR>
     * O Toyobo Biochemicals (9/05)<BR>
     * Q Molecular Biology Resources (8/05)<BR>
     * R Promega Corporation (9/05)<BR>
     * S Sigma Chemical Corporation (9/05)<BR>
     * U Bangalore Genei (9/05)<BR>
     * V Vivantis Technologies (1/06)<BR>
     * X EURx Ltd. (9/05)<BR>
     * Y CinnaGen Inc. (9/05)
     * </P>
     *
     * @param enzyme a <code>RestrictionEnzyme</code>.
     *
     * @return a <code>String</code>.
     * @since 1.5
     */
    public static String getSuppliers(RestrictionEnzyme enzyme)
    {
        if (! enzymeToSuppliers.containsKey(enzyme))
            return "";
        return (String) enzymeToSuppliers.get(enzyme);
    }

    /**
     * <code>register</code> regisiters a new
     * <code>RestrictionEnzyme</code> with the manager. It does not
     * check that the isoschizomers are known to the manager. If there
     * are custom isoschizomers in the <code>Set</code>, they should
     * be also be registered.
     *
     * @param enzyme a <code>RestrictionEnzyme</code> to register.
     *
     * @param isoschizomers a <code>Set</code> of
     * <code>RestrictionEnzyme</code>s which are isoschizomers.
     */
    public synchronized static void register(RestrictionEnzyme enzyme,
                                             Set               isoschizomers)
    {
        for (Iterator ii = isoschizomers.iterator(); ii.hasNext();)
        {
            Object o = ii.next();

            if (! (o instanceof RestrictionEnzyme))
            {
                throw new IllegalArgumentException("Isoschizomers set may contain only RestrictionEnzymes. Found '"
                                                   + o
                                                   + "'");
            }
        }

        registerEnzyme(enzyme);

        String name = enzyme.getName();
        nameToIsoschizomers.put(name, isoschizomers);
    }

    /**
     * <code>registerEnzyme</code> registers an enzyme, but does not
     * populate its isoschizomers. This is because registering the
     * contents of a REBASE file and registering a custom enzyme
     * handle addition of isoschizomers differently, but both use this
     * method for all other registration functions.
     *
     * @param enzyme a <code>RestrictionEnzyme</code>.
     */
    private static void registerEnzyme(RestrictionEnzyme enzyme)
    {
        String name = enzyme.getName();
        nameToEnzyme.put(name, enzyme);

        Integer sizeKey = new Integer(enzyme.getRecognitionSite().length());
        if (sizeToCutters.containsKey(sizeKey))
        {
            Set s = (Set) sizeToCutters.get(sizeKey);
            s.add(enzyme);
        }
        else
        {
            Set s = new HashSet();
            s.add(enzyme);
            sizeToCutters.put(sizeKey, s);
        }

        Pattern forward = Pattern.compile(enzyme.getForwardRegex());
        Pattern reverse = Pattern.compile(enzyme.getReverseRegex());
        enzymeToPattern.put(enzyme, new Pattern [] { forward, reverse });

        Annotation annotation = new SmallAnnotation();
        try
        {
            annotation.setProperty("dbxref", "REBASE:" + name);
        }
        catch (ChangeVetoException cve)
        {
            throw new BioError("Assertion Failure: failed to modify Annotation", cve);
        }

        annotation.addChangeListener(ChangeListener.ALWAYS_VETO,ChangeType.UNKNOWN);
        enzymeToAnnotation.put(enzyme, annotation);
    }

    private static void loadData(InputStream is) {
        nameToSite          = new HashMap();
        nameToEnzyme        = new HashMap();
        nameToIsoschizomers = new HashMap();
        sizeToCutters       = new HashMap();
        enzymeToPattern     = new HashMap();
        enzymeToAnnotation  = new HashMap();
        enzymeToSuppliers   = new HashMap();
        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(is));

            // Basic linesplit parser
            LineSplitParser lsParser = new LineSplitParser();
            lsParser.setEndOfRecord(TagValueParser.EMPTY_LINE_EOR);
            lsParser.setSplitOffset(3);
            lsParser.setContinueOnEmptyTag(true);
            lsParser.setMergeSameTag(true);

            // The end of the chain
            RebaseEnzymeBuilder builder = new RebaseEnzymeBuilder();

            // Create isoschizomer value splitter
            RegexSplitter iso =
                new RegexSplitter(Pattern.compile("([^,]+)"), 1);
            // Create site value splitter
            RegexSplitter site =
                new RegexSplitter(Pattern.compile("(\\(-?\\d+/-?\\d+\\)|[A-Za-z^]+)"), 1);

            ChangeTable table = new ChangeTable();
            table.setSplitter(REBASE_TAG_ISZR, iso);
            table.setSplitter(REBASE_TAG_SITE, site);
            ValueChanger changer = new ValueChanger(builder, table);

            // Filter tags
            TagDropper rebaseTags = new TagDropper(changer);
            // Retain the enzyme name
            rebaseTags.addTag(REBASE_TAG_NAME);
            // Retain isoschizomers
            rebaseTags.addTag(REBASE_TAG_ISZR);
            // Retain recognition sequence
            rebaseTags.addTag(REBASE_TAG_SITE);
            // Retain commercial supplier
            rebaseTags.addTag(REBASE_TAG_COMM);



            Parser parser = new Parser();
            while (parser.read(br, lsParser, rebaseTags))
            {
                continue;
            }

            // Replace isoschizomer names with RestrictionEnzymes
            Map tempMap = new HashMap();
            Set tempSet = null;
            for (Iterator ni = nameToIsoschizomers.keySet().iterator(); ni.hasNext();)
            {
                Object name = ni.next();
                Set isoschizomers = (Set) nameToIsoschizomers.get(name);

                if (isoschizomers.size() == 0)
                    tempSet = Collections.EMPTY_SET;
                else
                    tempSet = (Set) isoschizomers.getClass().newInstance();

                tempMap.put(name, tempSet);

                for (Iterator ii = isoschizomers.iterator(); ii.hasNext();) {
                    String isoName = (String) ii.next();
                    Object re = nameToEnzyme.get(isoName);
                    //bug fix suggested by George Waldon
                    if(re!=null)
                        tempSet.add(re);
                }
            }

            nameToIsoschizomers = tempMap;
        }
        catch (Exception e)
        {
            throw new BioError("Failed to read REBASE data file",e);
        }
    }

    /**
     * <code>RebaseEnzymeBuilder</code> creates enzyme instances and
     * populates the maps.
     */
    private static class RebaseEnzymeBuilder implements TagValueListener
    {
        private String recseq;
        private String name;
        private Set isoschizomers;
        private List isoBuffer;
        private SymbolList site;
        private int [] usCutPositions;
        private int [] dsCutPositions;
        private boolean isCommerciallyAvailable;

        private String tagState;
        private String suppliers;
        private boolean unknownSite;

        RebaseEnzymeBuilder() { }

        public void startRecord() throws ParserException
        {
            isoBuffer = new ArrayList(30);
            recseq         = "";
            site           = null;
            dsCutPositions = null;
            usCutPositions = null;
            unknownSite = false;
            isCommerciallyAvailable = false;
        }

        public void endRecord() throws ParserException
        {
            if (! getRecordState())
                return;
            if (unknownSite || site == null)
                return;

            int isoCount = isoBuffer.size();
            if (isoCount < 30)
            {
                isoschizomers = new SmallSet(isoCount);
                for (int i = 0; i < isoCount; i++)
                    isoschizomers.add(isoBuffer.get(i));
            }
            else
            {
                isoschizomers = new HashSet(isoBuffer);
            }

            if(!loadCommercialOnly || isCommerciallyAvailable) {
                RestrictionEnzyme re = createEnzyme();
                registerEnzyme(re);
                nameToIsoschizomers.put(name, isoschizomers);
                enzymeToSuppliers.put(re,suppliers);
                nameToSite.put(name,recseq);
            }
        }

        public void startTag(Object tag) throws ParserException
        {
            tagState = (String) tag;
        }

        public void endTag() throws ParserException { }

        public void value(TagValueContext context, Object value)
            throws ParserException
        {
            if (tagState.equals(REBASE_TAG_NAME))
                name = (String) value;
            else if (tagState.equals(REBASE_TAG_ISZR))
                isoBuffer.add(value);
            else if (tagState.equals(REBASE_TAG_SITE)) {
                recseq += (String) value;
                processSite(value);
            } else if (tagState.equals(REBASE_TAG_COMM))
                processSuppliers(value);
            else
                throw new ParserException("Unable to handle value for tag '"
                                          + tagState
                                          + "'");
        }

        boolean getRecordState()
        {
            return tagState != null;
        }

        RestrictionEnzyme createEnzyme()
        {
            RestrictionEnzyme enzyme = null;

            try
            {
                if (usCutPositions != null)
                {
                    enzyme = new RestrictionEnzyme(name, site,
                                                   usCutPositions[0],
                                                   usCutPositions[1],
                                                   dsCutPositions[0],
                                                   dsCutPositions[1]);
                }
                else
                {
                    enzyme = new RestrictionEnzyme(name, site,
                                                   dsCutPositions[0],
                                                   dsCutPositions[1]);
                }
            }
            catch (IllegalAlphabetException iae)
            {
                throw new BioError("New DNA SymbolList no longer consists on DNA Alphabet",iae);
            }

            return enzyme;
        }

        private void processSuppliers(Object value) throws ParserException {
            suppliers = (String) value;
            if(suppliers.length()!=0)
                isCommerciallyAvailable = true;
        }

        private void processSite(Object value) throws ParserException
        {
            StringBuffer sb = new StringBuffer((String) value);
            int div, forIdx, revIdx;

            // REBASE marks enzymes whose site is not known with '?'
            if (sb.charAt(0) == '?')
            {
                unknownSite = true;
                return;
            }

            if (sb.charAt(0) == '(')
            {
                // Index separator
                div = sb.indexOf("/");

                try
                {
                    forIdx = Integer.parseInt(sb.substring(1, div));
                    revIdx = Integer.parseInt(sb.substring(div + 1,
                                                           sb.length() - 1));
                }
                catch (NumberFormatException nfe)
                {
                    throw new ParserException("Failed to parse cut site index",nfe);
                }

                // Indices before the site indicate a double cutter
                if (site == null)
                {
                    usCutPositions = new int [2];
                    usCutPositions[0] = forIdx;
                    usCutPositions[1] = revIdx;
                }
                else
                {
                    dsCutPositions = new int [2];
                    dsCutPositions[0] = forIdx + site.length();
                    dsCutPositions[1] = revIdx + site.length();
                }
            }
            else
            {
                // Explicit cut site marker
                int cut = sb.indexOf("^");
                dsCutPositions = new int [2];

                try
                {
                    if (cut == -1)
                    {
                        site = DNATools.createDNA(sb.substring(0));
                        dsCutPositions[0] = 1;
                        dsCutPositions[1] = 1;
                    }
                    else
                    {
                        sb.deleteCharAt(cut);
                        site = DNATools.createDNA(sb.substring(0));
                        dsCutPositions[0] = cut;
                        dsCutPositions[1] = site.length() - cut;
                    }
                }
                catch (IllegalSymbolException iae)
                {
                    throw new ParserException("Illegal DNA symbol in recognition site",iae);
                }
            }
        }
    }
}
