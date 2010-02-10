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

import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.alignment.SimpleAlignment;
import org.biojava.bio.search.SearchContentHandler;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.homol.SimilarityPairFeature;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * <p><code>SimilarityPairBuilder</code> annotates query and subject
 * <code>Sequence</code> with <code>SimilarityPairFeature</code>s
 * created from SAX events supplied via a
 * <code>SeqSimilarityAdapter</code>. The objective is to describe a
 * simple pairwise relationship between the two sequences. This
 * differs slightly from using <code>HomologyFeature</code>s which are
 * slightly heavier, have to contain a full alignment and don't have
 * an explicit distinction between query and subject sequences in the
 * alignment. The SAX events should describe elements conforming to
 * the BioJava BlastLikeDataSetCollection DTD. Suitable sources are
 * <code>BlastLikeSAXParser</code> or <code>FastaSAXParser</code>.</p>
 *
 * <p>Annotated <code>ViewSequence</code>s wrapping both query and
 * subject sequences are created.</p>
 *
 * <p><strong>The current implementation should be used with care on
 * streams containing more than one search output</strong>. This is
 * because the builder will not stop after each report has been
 * processed and as a result all the subject sequences get
 * instantiated and a large object network could be created during
 * processing.</p>
 *
 * @author Keith James
 * @author Greg Cox
 * @since 1.2
 */
public class SimilarityPairBuilder extends ViewSequenceFactory
    implements SearchContentHandler
{
    /**
     * Constant <code>SIMILARITY_PAIR_FEATURE_TYPE</code> the type
     * String used by <code>SimilarityPairBuilder</code> when creating
     * <code>SimilarityPairFeature</code>s. This is the String which
     * is returned when a <code>SimilarityPairFeature</code>'s
     * <code>getType()</code> method is called.
     */
    public static final String SIMILARITY_PAIR_FEATURE_TYPE = "similarity";

    // Identifiers for query and database
    private String queryID;

    // Data holders for search result properties
    private Map resultData;
    private Map hitData;
    private Map subHitData;

    private SymbolTokenization tokenParser;
    private StringBuffer       tokenBuffer;

    // Flag indicating whether there are more results in the stream
    private boolean moreSearchesAvailable = false;

    public SimilarityPairBuilder()
    {
        resultData       = new HashMap();
        hitData          = new HashMap();
        subHitData       = new HashMap();
        queryViewCache   = new HashMap();
        subjectViewCache = new HashMap();
        tokenBuffer      = new StringBuffer(1024);
    }

    public Sequence getAnnotatedQuerySeq(String queryID)
        throws IllegalIDException
    {
        if (! queryViewCache.containsKey(queryID))
            throw new IllegalIDException("Failed to retrieve annotated query sequence from cache using ID '"
                                         + queryID
                                         + "' (unknown ID");

        return (Sequence) queryViewCache.get(queryID);
    }

    public Sequence getAnnotatedSubjectSeq(String subjectID)
        throws IllegalIDException
    {
        if (! subjectViewCache.containsKey(subjectID))
            throw new IllegalIDException("Failed to retrieve annotated subject sequence from cache using ID '"
                                         + subjectID
                                         + "' (unknown ID");

        return (Sequence) subjectViewCache.get(subjectID);
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
    }

    public void endSearch() { }

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
            makeSimilarity();
        }
        catch (BioException be)
        {
            System.err.println("Failed to build Similarity:");
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

    private void makeSimilarity() throws BioException
    {
        subHitData.putAll(resultData);
        subHitData.putAll(hitData);

        // Try to get a valid TokenParser
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

        // Set strands of hit on query and subject
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

        // Get start/end
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

        Sequence   queryView = makeQueryViewSequence(queryID);

        // Map of Alignment sequences
        Map labelMap = new HashMap();

        try
        {
            // Set source to the program name
            String source = "unknown";
            if (subHitData.containsKey("program"))
                source = (String) subHitData.get("program");

            tokenBuffer.setLength(0);
            tokenBuffer.append((String) subHitData.get("querySequence"));
            labelMap.put(SimilarityPairFeature.QUERY_LABEL,
                         new SimpleSymbolList(tokenParser, tokenBuffer.substring(0)));

            tokenBuffer.setLength(0);
            tokenBuffer.append((String) subHitData.get("subjectSequence"));
            labelMap.put(SimilarityPairFeature.SUBJECT_LABEL,
                         new SimpleSymbolList(tokenParser, tokenBuffer.substring(0)));

            double score = 0.0;
            if (subHitData.containsKey("score"))
                score = Double.parseDouble((String) subHitData.get("score"));

            // Query sequence feature
            SimilarityPairFeature.Template qt =
                new SimilarityPairFeature.Template();
            qt.type       = SIMILARITY_PAIR_FEATURE_TYPE;
            qt.source     = source;
            qt.location   = new RangeLocation(qStart, qEnd);
            qt.strand     = qStrand;
            qt.score      = score;
            qt.annotation = AnnotationFactory.makeAnnotation(subHitData);

            // Subject sequence feature
            SimilarityPairFeature.Template st =
                new SimilarityPairFeature.Template();
            st.type       = SIMILARITY_PAIR_FEATURE_TYPE;
            st.source     = source;
            st.location   = new RangeLocation(sStart, sEnd);
            st.strand     = sStrand;
            st.score      = score;
            st.annotation = AnnotationFactory.makeAnnotation(subHitData);

            Alignment a = new SimpleAlignment(labelMap);
            qt.alignment = a;
            st.alignment = a;

            SimilarityPairFeature qf =
                (SimilarityPairFeature) queryView.createFeature(qt);

            SimilarityPairFeature sf =
                (SimilarityPairFeature) queryView.createFeature(qt);

            sf.setSibling(qf);
            qf.setSibling(sf);

            qf.addChangeListener(ChangeListener.ALWAYS_VETO,
                                 ChangeType.UNKNOWN);
            sf.addChangeListener(ChangeListener.ALWAYS_VETO,
                                 ChangeType.UNKNOWN);
        }
        catch (ChangeVetoException cve)
        {
            throw new BioError("Assertion failure creating "
                               + "SimilarityPairFeature. Template "
                               + "modification vetoed",cve);
        }
    }
}
