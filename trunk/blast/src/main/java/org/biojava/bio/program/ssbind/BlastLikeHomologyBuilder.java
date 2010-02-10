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

package org.biojava.bio.program.ssbind;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.alignment.SimpleAlignment;
import org.biojava.bio.search.SearchContentHandler;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.biojava.bio.seq.homol.Homology;
import org.biojava.bio.seq.homol.HomologyFeature;
import org.biojava.bio.seq.homol.SimpleHomology;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.utils.ChangeVetoException;

/**
 * <p><code>BlastLikeHomologyBuilder</code> populates a
 * <code>List</code> with <code>Homology</code> instances created from
 * SAX events supplied via a <code>SeqSimilarityAdapter</code>. The
 * SAX events should describe elements conforming to the BioJava
 * BlastLikeDataSetCollection DTD. Suitable sources are
 * <code>BlastLikeSAXParser</code> or
 * <code>FastaSAXParser</code>. Annotated <code>ViewSequence</code>s
 * wrapping both query and subject sequences are created and populated
 * with <code>HomologyFeature</code>s. See the documentation of
 * <code>Homology</code> and <code>HomologyFeature</code>.</p>
 *
 * <p>As <code>SimpleHomologyFeature</code>s are created on views of
 * the query and subject sequences, both query and subject should be
 * nucleotide sequences (<code>SimpleHomologyFeature</code> extends
 * <code>StrandedFeature</code>.). This limits the searches currently
 * handled to BLASTN, TBLASTX and Fasta DNA.</p>
 *
 * @author Keith James
 * @author Greg Cox
 * @since 1.2
 */
public class BlastLikeHomologyBuilder extends ViewSequenceFactory
    implements SearchContentHandler
{
    /**
     * <code>HOMOLOGY_FEATURE_TYPE</code> is the type String used by
     * <code>BlastLikeHomologyBuilder</code> when creating
     * <code>HomologyFeature</code>s. This is the String which is
     * returned when an <code>HomologyFeature</code>'s
     * <code>getType()</code> method is called.
     */
    public static final String HOMOLOGY_FEATURE_TYPE = "homology";

    // Identifiers for query and database
    private String queryID;

    // Data holders for search result properties
    private Map resultData;
    private Map hitData;
    private Map subHitData;

    private SymbolTokenization tokenParser;
    private StringBuffer       tokenBuffer;

    // List for holding homologies from current search. There may be
    // more than one search result in a stream
    private List homologies;
    // Flag indicating whether there are more results in the stream
    private boolean moreSearchesAvailable = false;
    // List to accept homologies from all results in the stream
    private List target;

    /**
     * Creates a new <code>BlastLikeHomologyBuilder</code> which will
     * instantiate <code>Homology</code> objects into the
     * <code>List</code> target.
     *
     * @param target a <code>List</code>.
     */
    public BlastLikeHomologyBuilder(List target)
    {
        this.target = target;

        resultData       = new HashMap();
        hitData          = new HashMap();
        subHitData       = new HashMap();
        queryViewCache   = new HashMap();
        subjectViewCache = new HashMap();
        tokenBuffer      = new StringBuffer(1024);
    }

    public void setQueryID(String queryID)
    {
        this.queryID = queryID;
    }

    public boolean getMoreSearches()
    {
        return moreSearchesAvailable;
    }

    public void setMoreSearches(boolean value)
    {
        moreSearchesAvailable = value;
    }

    public void startSearch()
    {
        subjectViewCache.clear();
        homologies = new ArrayList();
    }

    public void endSearch()
    {
        target.addAll(homologies);
    }

    public void startHeader()
    {
        resultData.clear();
    }

    public void endHeader() { }

    public void startHit()
    {
        hitData.clear();
        subHitData.clear();
    }

    public void endHit() { }

    public void startSubHit() { }

    public void endSubHit()
    {
        try
        {
            homologies.add(makeHomology());
        }
        catch (BioException be)
        {
            System.err.println("Failed to build Homology:");
            be.printStackTrace();
        }
    }

    public void addSearchProperty(Object key, Object value)
    {
        resultData.put(key, value);
    }

    public void addHitProperty(Object key, Object value)
    {
        hitData.put(key, value);
    }

    public void addSubHitProperty(Object key, Object value)
    {
        subHitData.put(key, value);
    }

    /**
     * <code>makeHomology</code> creates a new
     * <code>SimpleHomology</code> describing the similarity between
     * the query and subject sequences. The
     * <code>HomologyFeatures</code> created are added to
     * <code>ViewSequence</code>s wrapping the query and subject
     * sequences.
     *
     * @return an <code>Homology</code>.
     *
     * @exception BioException if an error occurs.
     */
    private Homology makeHomology() throws BioException
    {
        subHitData.putAll(resultData);
        subHitData.putAll(hitData);

        if (tokenParser == null)
        {
            String identifier;
            // Try explicit sequence type first
            if (subHitData.containsKey("hitSequenceType"))
                identifier = (String) subHitData.get("hitSequenceType");
            // Otherwise try to resolve from the program name (only
            // works for Blast)
            else if (subHitData.containsKey("program"))
                identifier = (String) subHitData.get("program");
            else
                throw new BioException("Failed to determine sequence type");

            FiniteAlphabet alpha = AlphabetResolver.resolveAlphabet(identifier);
            tokenParser = alpha.getTokenization("token");
        }

        Strand qStrand = StrandedFeature.POSITIVE;
        Strand sStrand = StrandedFeature.POSITIVE;

        // In cases where an explicit strand is given (FASTA DNA, BLASTN)
        if (subHitData.containsKey("queryStrand") &&
            subHitData.get("queryStrand").equals("minus"))
            qStrand = StrandedFeature.NEGATIVE;

        if (subHitData.containsKey("subjectStrand") &&
            subHitData.get("subjectStrand").equals("minus"))
            sStrand = StrandedFeature.NEGATIVE;

        // In cases where a frame is given as this contains strand
        // information (TBLASTN for hit, TBLASTX for both query and
        // hit)
        if (subHitData.containsKey("queryFrame") &&
            ((String) subHitData.get("queryFrame")).startsWith("minus"))
            qStrand = StrandedFeature.NEGATIVE;

        if (subHitData.containsKey("subjectFrame") &&
            ((String) subHitData.get("subjectFrame")).startsWith("minus"))
            sStrand = StrandedFeature.NEGATIVE;

        int qStart = Integer.parseInt((String) subHitData.get("querySequenceStart"));
        int   qEnd = Integer.parseInt((String) subHitData.get("querySequenceEnd"));
        int sStart = Integer.parseInt((String) subHitData.get("subjectSequenceStart"));
        int   sEnd = Integer.parseInt((String) subHitData.get("subjectSequenceEnd"));

        // The start/end coordinates from BioJava XML don't follow the
        // BioJava paradigm of start < end, with orientation given by
        // the strand property. Rather, they present start/end as
        // displayed in BLAST output, with the coordinates being
        // inverted on the reverse strand. We account for this here.
        if (qStrand == StrandedFeature.NEGATIVE)
        {
            int swap = qStart;
            qStart = qEnd;
            qEnd   = swap;
        }

        if (sStrand == StrandedFeature.NEGATIVE)
        {
            int swap = sStart;
            sStart = sEnd;
            sEnd   = swap;
        }

        String subjectID = (String) hitData.get("subjectId");

        Sequence   queryView = makeQueryViewSequence(queryID);
        Sequence subjectView = makeSubjectViewSequence(subjectID);

        SimpleHomology homology = new SimpleHomology();

        // Map of HomologyFeatures to Alignment sequences
        Map labelMap = new HashMap();

        try
        {
            String source = "unknown";
            if (subHitData.containsKey("program"))
                source = (String) subHitData.get("program");

            tokenBuffer.setLength(0);
            tokenBuffer.append((String) subHitData.get("querySequence"));

            // Query sequence feature
            HomologyFeature.Template qt = new HomologyFeature.Template();
            qt.type       = HOMOLOGY_FEATURE_TYPE;
            qt.source     = source;
            qt.location   = new RangeLocation(qStart, qEnd);
            qt.strand     = qStrand;
            qt.annotation = AnnotationFactory.makeAnnotation(subHitData);
            qt.homology   = homology;

            // Map the new feature to the alignment SymbolList
            labelMap.put(queryView.createFeature(qt),
                         new SimpleSymbolList(tokenParser, tokenBuffer.substring(0)));

            tokenBuffer.setLength(0);
            tokenBuffer.append((String) subHitData.get("subjectSequence"));

            // Subject sequence feature
            HomologyFeature.Template st = new HomologyFeature.Template();
            st.type       = HOMOLOGY_FEATURE_TYPE;
            st.source     = source;
            st.location   = new RangeLocation(sStart, sEnd);
            st.strand     = sStrand;
            st.annotation = AnnotationFactory.makeAnnotation(subHitData);
            st.homology   = homology;

            // Map the new feature to the alignment SymbolList
            labelMap.put(subjectView.createFeature(st),
                         new SimpleSymbolList(tokenParser, tokenBuffer.substring(0)));

            Alignment a = new SimpleAlignment(labelMap);
            homology.setAlignment(a);

            return homology;
        }
        catch (ChangeVetoException cve)
        {
            throw new BioException( "Failed to create HomologyFeature",cve);
        }
    }
}
