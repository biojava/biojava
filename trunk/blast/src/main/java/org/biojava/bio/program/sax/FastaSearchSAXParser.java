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

package org.biojava.bio.program.sax;

import java.io.BufferedReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.search.SearchContentHandler;
import org.biojava.utils.ParserException;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;

/**
 * <p><code>FastaSearchSAXParser</code> is a SAX2 compliant parser for
 * '-m 10' format output from the the Fasta search program (see the
 * Fasta documentation for details of this format).</p>
 *
 * <p>Versions of Fasta supported are as follows. Note that the compile
 * time option -DM10_CONS should be used to allow correct reporting of
 * the number of matches in HSPSummary elements</p>
 *
 * <ul>
 *   <li>33t07</li>
 *   <li>33t08 (current tests are against output from this version)</li>
 * </ul>
 *
 * <p>The SAX2 events produced are as if the input to the parser was
 * an XML file validating against the BioJava
 * BlastLikeDataSetCollection DTD. There is no requirement for an
 * intermediate conversion of native output to XML format.</p>
 *
 * @author Keith James
 * @since 1.2
 */
public class FastaSearchSAXParser extends AbstractNativeAppSAXParser
    implements SearchContentHandler
{
    private FastaSearchParser fastaParser;
    private Map               searchProperties;
    private Map               hitProperties;

    private String queryID;
    private String databaseID;

    private AttributesImpl  attributes;
    private QName                qName;

    private boolean firstHit = true;

    // Set/reset by callback from main parser
    private boolean moreSearchesAvailable = true;

    // For formatting rounded numbers
    private NumberFormat nFormat;

    // Platform independent linebreaks
    private String nl;

    // For creating character events
    private StringBuffer props;
    private StringBuffer seqTokens;
    private String stringOut;
    private char [] charOut;

    /**
     * Creates a new <code>FastaSearchSAXParser</code> instance.
     */
    public FastaSearchSAXParser()
    {
        this.setNamespacePrefix("biojava");
        this.addPrefixMapping("biojava", "http://www.biojava.org");

        fastaParser = new FastaSearchParser();
        attributes  = new AttributesImpl();
        qName       = new QName(this);
        nFormat     = new DecimalFormat("###.0");
        props       = new StringBuffer(1024);
        seqTokens   = new StringBuffer(2048);
        nl          = System.getProperty("line.separator");
    }

    public void parse(InputSource source)
        throws IOException, SAXException
    {
        BufferedReader content = getContentStream(source);

        if (oHandler == null)
            throw new SAXException("Running FastaSearchSAXParser with null ContentHandler");

        try
        {
            attributes.clear();
            // Namespace attribute
            qName.setQName("xmlns");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    "");
            // Namespace attribute
            qName.setQName("xmlns:biojava");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    "http://www.biojava.org");

            // Start the BlastLikeDataSetCollection
            startElement(new QName(this, this.prefix("BlastLikeDataSetCollection")),
                         (Attributes) attributes);

            while (moreSearchesAvailable)
            {
                // This method returns once a single result is
                // parsed. The parser also informs us of subsequent
                // results via the setMoreSearches() method.
                fastaParser.parseSearch(content, this);
            }

            // End the BlastLikeDataSetCollection
            endElement(new QName(this, this.prefix("BlastLikeDataSetCollection")));
        }
        catch (BioException be)
        {
            throw new SAXException(be);
        }
        catch (ParserException pe)
        {
            throw new SAXException(pe);
        }
    }

    public boolean getMoreSearches()
    {
        return moreSearchesAvailable;
    }

    public void setMoreSearches(boolean value)
    {
        moreSearchesAvailable = value;
    }

    /**
     * <code>setQuerySeq</code> identifies the query sequence by a
     * name, ID or URN.
     *
     * @param identifier a <code>String</code> which should be an
     * unique identifer for the sequence.
     *
     * @deprecated use <code>setQueryID</code> instead.
     */
    public void setQuerySeq(String identifier)
    {
        setQueryID(identifier);
    }

    public void setQueryID(String queryID)
    {
        this.queryID = queryID;
    }

    /**
     * <code>setSubjectDB</code> identifies the database searched by a
     * name, ID or URN.
     *
     * @param identifier a <code>String</code> which should be an unique
     * identifier for the database searched.
     *
     * @deprecated use <code>setDatabaseID</code> instead.
     */
    public void setSubjectDB(String identifier)
    {
        setDatabaseID(identifier);
    }

    public void setDatabaseID(String databaseID)
    {
        this.databaseID = databaseID;
    }

    public void startSearch()
    {
        searchProperties = new HashMap();
    }

    public void addSearchProperty(Object key, Object value)
    {
        searchProperties.put(key, value);
    }

    public void endSearch()
    {
        try
        {
	    // If we found any hits then we need to close a Detail
            // element too
            if (! firstHit)
            {
                endElement(new QName(this, this.prefix("Detail")));

                // Prime to get next firstHit
                firstHit = true;
            }

            endElement(new QName(this, this.prefix("BlastLikeDataSet")));
        }
        catch (SAXException se)
        {
            System.err.println("An error occurred while creating an endElement SAX event: ");
            se.printStackTrace();
        }
    }

    public void startHeader() { }

    public void endHeader()
    {
        try
        {
            attributes.clear();
            // Program name attribute
            qName.setQName("program");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    (String) searchProperties.get("pg_name"));

            // Program version attribute
            qName.setQName("version");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    (String) searchProperties.get("pg_ver"));

            // Start the BlastLikeDataSet
            startElement(new QName(this, this.prefix("BlastLikeDataSet")),
                         (Attributes) attributes);

            attributes.clear();
            // Start the Header
            startElement(new QName(this, this.prefix("Header")),
                         (Attributes) attributes);

            attributes.clear();
            // Query id attribute
            qName.setQName("id");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    queryID);

            // metaData attribute for QueryId
            qName.setQName("metaData");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    "none");

            // Start the QueryId
            startElement(new QName(this, this.prefix("QueryId")),
                         (Attributes) attributes);
            // End the QueryId
            endElement(new QName(this, this.prefix("QueryId")));

            attributes.clear();
            // id attribute for DatabaseId
            qName.setQName("id");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    databaseID);

            // metaData attribute for DatabaseId
            qName.setQName("metaData");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    "none");

            // Start the DatabaseId
            startElement(new QName(this, this.prefix("DatabaseId")),
                         (Attributes) attributes);
            // End the DatabaseId
            endElement(new QName(this, this.prefix("DatabaseId")));

            attributes.clear();
            // Whitespace attribute for raw data
            qName.setQName("xml:space");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    "preserve");

            // Start the RawOutput
            startElement(new QName(this, this.prefix("RawOutput")),
                         (Attributes) attributes);

            // Reconstitute the 'raw' header from the properties Map
            Set spKeys = searchProperties.keySet();

            String [] searchPropKeys =
                (String []) spKeys.toArray(new String [spKeys.size() - 1]);
            Arrays.sort(searchPropKeys);

            // Clear StringBuffer
            props.setLength(0);

            props.append(nl);
            for (int i = 0; i < searchPropKeys.length; i++)
            {
                props.append(searchPropKeys[i] + ": ");
                props.append((String) searchProperties.get(searchPropKeys[i]) + nl);
            }

            charOut = new char [props.length()];
            props.getChars(0, props.length(), charOut, 0);

            // Characters of raw header
            characters(charOut, 0, charOut.length);

            // End the RawOutput
            endElement(new QName(this, this.prefix("RawOutput")));

            // End the Header
            endElement(new QName(this, this.prefix("Header")));
        }
        catch (SAXException se)
        {
            System.err.println("An error occurred while creating SAX events from header data: ");
            se.printStackTrace();
        }
    }

    public void startHit()
    {
        // Hit elements must be wrapped in a Detail element so we
        // start one at the first hit
        if (firstHit)
        {
            firstHit = false;
            attributes.clear();

            try
            {
                startElement(new QName(this, this.prefix("Detail")),
                             (Attributes) attributes);
            }
            catch (SAXException se)
            {
                System.err.println("An error occurred while creating startElement SAX event from hit data: ");
                se.printStackTrace();
            }
        }

        hitProperties = new HashMap();
    }

    public void addHitProperty(Object key, Object value)
    {
        hitProperties.put(key, value);
    }

    public void endHit() { }

    public void startSubHit() { }

    public void addSubHitProperty(Object key, Object value)
    {
        hitProperties.put(key, value);
    }

    public void endSubHit()
    {
        attributes.clear();
        // Query sequence length attribute
        qName.setQName("sequenceLength");
        attributes.addAttribute(qName.getURI(),
                                qName.getLocalName(),
                                qName.getQName(),
                                "CDATA",
                                (String) hitProperties.get("subject_sq_len"));

        try
        {
            // Start the Hit
            startElement(new QName(this, this.prefix("Hit")),
                         (Attributes) attributes);

            attributes.clear();
            // Hit id attribute
            qName.setQName("id");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    (String) hitProperties.get("id"));
            // Metadata attribute
            qName.setQName("metaData");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    "none");
            // Start the HitId
            startElement(new QName(this, this.prefix("HitId")),
                         (Attributes) attributes);
            // End the HitId
            endElement(new QName(this, this.prefix("HitId")));

            attributes.clear();
            // Start the HitDescription
            startElement(new QName(this, this.prefix("HitDescription")),
                         (Attributes) attributes);

            stringOut = (String) hitProperties.get("desc");

            charOut = new char [stringOut.length()];
            stringOut.getChars(0, stringOut.length(), charOut, 0);

            // Characters of description
            characters(charOut, 0, charOut.length);

            // End the HitDescription
            endElement(new QName(this, this.prefix("HitDescription")));

            // Start the HSPCollection
            startElement(new QName(this, this.prefix("HSPCollection")),
                         (Attributes) attributes);

            // Start the HSP (for Fasta, we use one "HSP" to represent the hit
            startElement(new QName(this, this.prefix("HSP")),
                         (Attributes) attributes);

            String score;
            if (hitProperties.containsKey("fa_z-score"))
                score = (String) hitProperties.get("fa_z-score");
            else
                throw new SAXException("Failed to retrieve hit z-score from search data");

            // Score attribute
            qName.setQName("score");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    score);
            // expectValue attribute
            qName.setQName("expectValue");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    (String) hitProperties.get("fa_expect"));
            // numberOfIdentities attribute
            qName.setQName("numberOfIdentities");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    countTokens(':', (String) hitProperties.get("matchTokens")));

            String overlap;
            if (hitProperties.containsKey("fa_overlap"))
                overlap = hitProperties.get("fa_overlap").toString();
            else if (hitProperties.containsKey("sw_overlap"))
                overlap = hitProperties.get("sw_overlap").toString();
            else
                throw new SAXException("Failed to retrieve hit overlap from search data");

            // alignmentSize attribute
            qName.setQName("alignmentSize");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    overlap);

            float percentId;
            if (hitProperties.containsKey("fa_ident"))
                percentId = Float.parseFloat((String) hitProperties.get("fa_ident"));
            else
                percentId = Float.parseFloat((String) hitProperties.get("sw_ident"));

            // percentageIdentity attribute
            qName.setQName("percentageIdentity");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    nFormat.format(percentId * 100));

            // Maybe proper RNA check? Should be same for query and subject
            String seqType;
            if (hitProperties.get("query_sq_type").equals("dna"))
                seqType = "dna";
            else
                seqType = "protein";

            // querySequenceType attribute
            qName.setQName("querySequenceType");
            attributes.addAttribute(qName.getURI(),
				    qName.getLocalName(),
				    qName.getQName(),
				    "CDATA",
				    seqType);

            // Maybe proper RNA check? Maybe raise exception if not
            // same as query type?
            if (hitProperties.get("subject_sq_type").equals("dna"))
                seqType = "dna";
            else
                seqType = "protein";

            // hitSequenceType attribute
            qName.setQName("hitSequenceType");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    seqType);

            // Strand information only valid for DNA
            if (seqType.equals("dna"))
            {
                // queryStrand attribute (always plus for Fasta)
                qName.setQName("queryStrand");
                attributes.addAttribute(qName.getURI(),
                                        qName.getLocalName(),
                                        qName.getQName(),
                                        "CDATA",
                                        "plus");

                String strand;
                if (hitProperties.get("fa_frame").equals("f"))
                    strand = "plus";
                else
                    strand = "minus";

                // hitStrand attribute (may be minus for Fasta vs. nt sequence)
                qName.setQName("hitStrand");
                attributes.addAttribute(qName.getURI(),
                                        qName.getLocalName(),
                                        qName.getQName(),
                                        "CDATA",
                                        strand);
            }

            // Start the HSPSummary
            startElement(new QName(this, this.prefix("HSPSummary")),
                         (Attributes) attributes);

            attributes.clear();
            // Start the RawOutput
            startElement(new QName(this, this.prefix("RawOutput")),
                         (Attributes) attributes);

            // Reconstitute the 'raw' header from the properties Map
            Set hpKeys = hitProperties.keySet();

            String [] hitPropKeys =
                (String []) hpKeys.toArray(new String [hpKeys.size() - 1]);
            Arrays.sort(hitPropKeys);

            // Clear StringBuffer
            props.setLength(0);
            props.append(nl);
            for (int i = 0; i < hitPropKeys.length; i++)
            {
                // Skip the sequence and consensus tokens
                if (hitPropKeys[i].endsWith("Tokens"))
                    continue;
                props.append(hitPropKeys[i] + ": ");
                props.append((String) hitProperties.get(hitPropKeys[i]) + nl);
            }

            charOut = new char [props.length()];
            props.getChars(0, props.length(), charOut, 0);

            // Characters of raw header
            characters(charOut, 0, charOut.length);

            // End the RawOutput
            endElement(new QName(this, this.prefix("RawOutput")));

            // End the HSPSummary
            endElement(new QName(this, this.prefix("HSPSummary")));

            // Start the BlastLikeAlignment
            startElement(new QName(this, this.prefix("BlastLikeAlignment")),
                         (Attributes) attributes);

            String alStart     = (String) hitProperties.get("query_al_start");
            String alStop      = (String) hitProperties.get("query_al_stop");
            String alDispStart = (String) hitProperties.get("query_al_display_start");

            // Query sequence startPosition attribute
            qName.setQName("startPosition");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    alStart);

            // Query sequence stopPosition attribute
            qName.setQName("stopPosition");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    alStop);

            // Start the QuerySequence
            startElement(new QName(this, this.prefix("QuerySequence")),
                         (Attributes) attributes);

            seqTokens.setLength(0);
            seqTokens.append((String) hitProperties.get("querySeqTokens"));

            // Fasta includes context sequence which we need to trim
            stringOut = prepSeqTokens(seqTokens,
                                      Integer.parseInt(alStart),
                                      Integer.parseInt(alStop),
                                      Integer.parseInt(alDispStart));

            charOut = new char [stringOut.length()];
            stringOut.getChars(0, stringOut.length(), charOut, 0);

            // Characters of QuerySequence
            characters(charOut, 0, charOut.length);

            // End the QuerySequence
            endElement(new QName(this, this.prefix("QuerySequence")));

            attributes.clear();
            // Whitespace attribute for MatchConsensus
            qName.setQName("xml:space");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    "preserve");

            // Start the MatchConsensus
            startElement(new QName(this, this.prefix("MatchConsensus")),
                         (Attributes) attributes);

            stringOut = ((String) hitProperties.get("matchTokens")).trim();
            
            charOut = new char [stringOut.length()];
            stringOut.getChars(0, stringOut.length(), charOut, 0);

            // Characters of MatchConsensus
            characters(charOut, 0, charOut.length);

            // End the MatchConsensus
            endElement(new QName(this, this.prefix("MatchConsensus")));

            alStart     = (String) hitProperties.get("subject_al_start");
            alStop      = (String) hitProperties.get("subject_al_stop");
            alDispStart = (String) hitProperties.get("subject_al_display_start");

            attributes.clear();
            // Hit sequence startPosition attribute
            qName.setQName("startPosition");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    alStart);

            // Hit sequence stopPosition attribute
            qName.setQName("stopPosition");
            attributes.addAttribute(qName.getURI(),
                                    qName.getLocalName(),
                                    qName.getQName(),
                                    "CDATA",
                                    alStop);

            // Start the HitSequence
            startElement(new QName(this, this.prefix("HitSequence")),
                         (Attributes) attributes);

            seqTokens.setLength(0);
            seqTokens.append((String) hitProperties.get("subjectSeqTokens"));

            // Fasta includes context sequence which we need to trim
            stringOut = prepSeqTokens(seqTokens,
                                      Integer.parseInt(alStart),
                                      Integer.parseInt(alStop),
                                      Integer.parseInt(alDispStart));

            charOut = new char [stringOut.length()];
            stringOut.getChars(0, stringOut.length(), charOut, 0);

            // Characters of HitSequence
            characters(charOut, 0, charOut.length);

            // End the HitSequence
            endElement(new QName(this, this.prefix("HitSequence")));

            // End the BlastLikeAlignment
            endElement(new QName(this, this.prefix("BlastLikeAlignment")));

            // End the HSP
            endElement(new QName(this, this.prefix("HSP")));

            // End the HSPCollection
            endElement(new QName(this, this.prefix("HSPCollection")));

            // End the hit
            endElement(new QName(this, this.prefix("Hit")));
        }
        catch (SAXException se)
        {
            System.err.println("An error occurred while creating SAX events from hit data: ");
            se.printStackTrace();
        }
    }

    /**
     * <code>countTokens</code> counts up the occurrences of a char in
     * a <code>String</code>.
     *
     * @param token a <code>char</code> to count.
     * @param string a <code>String</code> to count within.
     *
     * @return a <code>String</code> representation of the total count.
     */
    private String countTokens(char token, String string)
    {
        int count = 0;
        for (int i = string.length(); --i >= 0;)
        {
            if (string.charAt(i) == token)
                count++;
        }
        return String.valueOf(count);
    }

    /**
     * The <code>prepSeqTokens</code> method prepares the sequence
     * data extracted from the Fasta output. Two things need to be
     * done; firstly, the leading gaps are removed from the sequence
     * (these are just format padding and not really part of the
     * alignment) and secondly, as Fasta supplies some flanking
     * sequence context for its alignments, this must be removed
     * too. See the Fasta documentation for an explanation of the
     * format.
     *
     * @param name a <code>StringBuffer</code> containing the
     * unprepared sequence tokens.
     * @param alStart an <code>int</code> indicating the start
     * position of the alignment in the original sequence.
     * @param alStop an <code>int</code> indicating the stop
     * position of the alignment in the original sequence.
     * @param alDispStart an <code>int</code> indicating the start
     * of a flanking context in the original sequence.
     *
     * @return a <code>String</code> value consisting of a subsequence
     * containing only the interesting alignment.
     */
    private String prepSeqTokens(StringBuffer seqTokens,
                                 int          alStart,
                                 int          alStop,
                                 int          alDispStart)
    {
        // Strip leading gap characters
        while (seqTokens.charAt(0) == '-')
            seqTokens.deleteCharAt(0);
        
        int gapCount = 0;
        // Count gaps to add to number of chars returned
        for (int i = seqTokens.length(); --i >= 0;)
        {
            if (seqTokens.charAt(i) == '-')
                gapCount++;
        }

        // Calculate the position at which the real alignment
        // starts/stops, allowing for the gaps, which are not counted
        // in the numbering system
        return seqTokens.substring(alStart - alDispStart,
                                   alStop  - alDispStart + gapCount + 1);
    }
}
