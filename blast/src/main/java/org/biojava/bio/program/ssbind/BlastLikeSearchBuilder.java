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
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.alignment.SimpleAlignment;
import org.biojava.bio.search.SearchBuilder;
import org.biojava.bio.search.SeqSimilaritySearchHit;
import org.biojava.bio.search.SeqSimilaritySearchResult;
import org.biojava.bio.search.SeqSimilaritySearchSubHit;
import org.biojava.bio.search.SimpleSeqSimilaritySearchHit;
import org.biojava.bio.search.SimpleSeqSimilaritySearchResult;
import org.biojava.bio.search.SimpleSeqSimilaritySearchSubHit;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.db.SequenceDBInstallation;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.utils.SmallMap;

/**
 * <p><code>BlastLikeSearchBuilder</code> will create
 * <code>SeqSimilaritySearchResult</code>s from SAX events via a
 * <code>SeqSimilarityAdapter</code>. The SAX events should describe
 * elements conforming to the BioJava BlastLikeDataSetCollection
 * DTD. Suitable sources are <code>BlastLikeSAXParser</code> or
 * <code>FastaSearchSAXParser</code>. The result objects are placed in
 * the <code>List</code> supplied to the constructor.</p>
 *
 * <p>The start/end/strand of <code>SeqSimilaritySearchHit</code>s are
 * calculated from their constituent
 * <code>SeqSimilaritySearchSubHit</code>s as follows:</p>
 *
 * <ul>
 * <li>The query start is the lowest query start coordinate of its
 *     sub-hits, regardless of strand</li>
 * <li>The query end is the highest query end coordinate of its sub-hits,
 *     regardless of strand</li>
 * <li>The hit start is the lowest hit start coordinate of its sub-hits,
 *     regardless of strand</li>
 * <li>The hit end is the highest hit end coordinate of its sub-hits,
 *     regardless of strand</li>
 * <li>The query strand is null for protein sequences. Otherwise it is
 *     equal to the query strand of its sub-hits if they are all on the
 *     same strand, or <code>StrandedFeature.UNKNOWN</code> if the sub-hits
 *     have mixed query strands</li>
 * <li>The hit strand is null for protein sequences. Otherwise it is
 *     equal to the hit strand of its sub-hits if they are all on the same
 *     strand, or <code>StrandedFeature.UNKNOWN</code> if the sub-hits have
 *     mixed hit strands</li>
 * </ul>
 *
 * <p>
 * This class has special meanings for particular keys: if you want to
 * adapt this class for another parser, you will need to be aware of
 * this. These originate from and are fully described in the
 * BlastLikeDataSetCollection DTD.
 * </p>
 * <table>
 * <tr>
 *   <th>Key</th>
 *   <th>Meaning</th>
 * </tr>
 * <tr>
 *   <td>program</td>
 *   <td>either this value or the subjectSequenceType value must be set. This can take values
 *       acceptable to AlphabetResolver. These are BLASTN, BLASTP, BLASTX, TBLASTN,
 *       TBLASTX, DNA and PROTEIN. </td>
 * </tr>
 *
 * <tr>
 *   <td>databaseId</td>
 *   <td>Identifier of database searched (in SequenceDBInstallation).</td>
 * </tr>
 * <tr>
 *   <td>subjectSequenceType</td>
 *   <td>type of sequence that hit is. Can be DNA or PROTEIN.</td>
 * </tr>
 * <tr>
 *   <td>subjectId</td>
 *   <td>id of sequence that is hit</td>
 * </tr>
 * <tr>
 *   <td>subjectDescription</td>
 *   <td>description of sequence that is hit</td>
 * </tr>
 * <tr>
 *   <td>queryStrand</td>
 *   <td>Strandedness of query in alignment. Takes values of "plus" and "minus"</td>
 * </tr>
 * <tr>
 *   <td>subjectStrand</td>
 *   <td>Strandedness of query in alignment. Takes values of "plus" and "minus"</td>
 * </tr>
 * <tr>
 *   <td>queryFrame</td>
 *   <td>self-evident</td>
 * </tr>
 * <tr>
 *   <td>subjectFrame</td>
 *   <td>self-evident</td>
 * </tr>
 * <tr>
 *   <td>querySequenceStart</td>
 *   <td>self-evident</td>
 * </tr>
 * <tr>
 *   <td>querySequenceEnd</td>
 *   <td>self-evident</td>
 * </tr>
 * <tr>
 *   <td>subjectSequenceStart</td>
 *   <td>self-evident</td>
 * </tr>
 * <tr>
 *   <td>subjectSequenceEnd</td>
 *   <td>self-evident</td>
 * </tr>
 * <tr>
 *   <td>score</td>
 *   <td>self-evident</td>
 * </tr>
 * <tr>
 *   <td>expectValue</td>
 *   <td>self-evident</td>
 * </tr>
 * <tr>
 *   <td>pValue</td>
 *   <td>self-evident</td>
 * </tr>
 * </table>
 *
 * @author Keith James
 * @author Greg Cox
 * @since 1.2
 */
public class BlastLikeSearchBuilder implements SearchBuilder
{
    // Supplier of instances of searched databases
    private SequenceDBInstallation subjectDBs;
    // Holder for all query sequences
    private SequenceDB querySeqHolder;

    // The ID of the database searched
    private String databaseID;
    // The ID of the query sequence
    private String queryID;

    // Hit and Result annotation
    private Annotation resultAnnotation;

    // Data holders for search result properties
    private Map resultPreAnnotation;
    private Map searchParameters;
    private Map hitData;
    private Map subHitData;

    private SymbolTokenization tokenParser;

    private List hits;
    private List subHits;

    private SeqSimilaritySearchSubHit [] subs;

    // Flag indicating whether there are more results in the stream
    private boolean moreSearchesAvailable = false;

    // List to accept all results in the stream
    private List target;

    /**
     * Creates a new <code>BlastLikeSearchBuilder</code> which will
     * instantiate results into the <code>List</code> target.
     *
     * @param target a <code>List</code>.
     */
    public BlastLikeSearchBuilder(List target)
    {
        this.target = target;

        resultPreAnnotation = new HashMap();
        searchParameters    = new HashMap();
        hitData             = new HashMap();
        subHitData          = new HashMap();
    }

    /**
     * Creates a new <code>BlastLikeSearchBuilder</code> which will
     * instantiate results into the <code>List</code> target.
     *
     * @param target a <code>List</code>.
     * @param querySeqHolder a <code>SequenceDB</code> of query
     * sequences.
     * @param subjectDBs a <code>SequenceDBInstallation</code> of
     * databases searched.
     */
    public BlastLikeSearchBuilder(List                   target,
                                  SequenceDB             querySeqHolder,
                                  SequenceDBInstallation subjectDBs)
    {
        this(target);
        this.querySeqHolder = querySeqHolder;
        this.subjectDBs = subjectDBs;
    }

    public SeqSimilaritySearchResult makeSearchResult()
        throws BioException
    {
        if (querySeqHolder == null)
            throw new BioException("Running BlastLikeSearchBuilder with null query SequenceDB");

        if (subjectDBs == null)
            throw new BioException("Running BlastLikeSearchBuilder with null subject SequenceDB installation");

        Sequence query = querySeqHolder.getSequence(queryID);
        if (query == null)
            throw new BioException("Failed to retrieve query sequence from queryDB using ID '"
                                   + queryID
                                   + "' (sequence was null)");

        SequenceDB subjectDB = (SequenceDB) subjectDBs.getSequenceDB(databaseID);
        if (subjectDB == null)
            throw new BioException("Failed to retrieve database from installation using ID '"
                                   + databaseID
                                   + "' (database was null)");

        return new SimpleSeqSimilaritySearchResult(query,
                                                   subjectDB,
                                                   searchParameters,
                                                   hits,
                                                   resultAnnotation);
    }

    /**
     * <code>setQuerySeqHolder</code> sets the query sequence holder
     * to a specific database.
     *
     * @param querySeqHolder a <code>SequenceDB</code> containing the
     * query sequence(s).
     */
    public void setQuerySeqHolder(SequenceDB querySeqHolder)
    {
        this.querySeqHolder = querySeqHolder;
    }

    /**
     * <code>setSubjectDBInstallation</code> sets the subject database
     * holder to a specific installation.
     *
     * @param subjectDBs a <code>SequenceDBInstallation</code>
     * containing the subject database(s)
     */
    public void setSubjectDBInstallation(SequenceDBInstallation subjectDBs)
    {
        this.subjectDBs = subjectDBs;
    }

    public void setQueryID(String queryID)
    {
        this.queryID = queryID;
        addSearchProperty("queryId", queryID);
    }

    public void setDatabaseID(String databaseID)
    {
        this.databaseID = databaseID;
        addSearchProperty("databaseId", databaseID);
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
        hits = new ArrayList();
    }

    public void endSearch()
    {
        try
        {
            resultAnnotation = AnnotationFactory.makeAnnotation(resultPreAnnotation);
            target.add(makeSearchResult());
        }
        catch (BioException be)
        {
            System.err.println("Failed to build SeqSimilaritySearchResult:");
            be.printStackTrace();
        }
    }

    public void startHeader()
    {
        resultPreAnnotation.clear();
        searchParameters.clear();
    }

    public void endHeader() { }

    public void startHit()
    {
        hitData.clear();
        subHits = new ArrayList();
    }

    public void endHit()
    {
        hits.add(makeHit());
    }

    public void startSubHit()
    {
        subHitData.clear();
    }

    public void endSubHit()
    {
        try
        {
            subHits.add(makeSubHit());
        }
        catch (BioException be)
        {
            be.printStackTrace();
        }
    }

    public void addSearchProperty(Object key, Object value)
    {
        resultPreAnnotation.put(key, value);
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
     * <code>makeHit</code> creates a new hit. The hit's strand data
     * is the same as that of the highest-scoring sub-hit. The hit's
     * start/end data are the same as the extent of the sub-hits on
     * that strand.
     *
     * @return a <code>SeqSimilaritySearchHit</code>.
     */
    private SeqSimilaritySearchHit makeHit()
    {
        double sc = Double.NaN;
        double ev = Double.NaN;
        double pv = Double.NaN;

        subs = (SeqSimilaritySearchSubHit []) subHits
            .toArray(new SeqSimilaritySearchSubHit [subHits.size() - 1]);

        // Sort to get highest score
        Arrays.sort(subs, SeqSimilaritySearchSubHit.byScore);
        sc = subs[subs.length - 1].getScore();
        ev = subs[subs.length - 1].getEValue();
        pv = subs[subs.length - 1].getPValue();

        // Check for any mixed or null strands
        boolean    mixQueryStrand = false;
        boolean  mixSubjectStrand = false;
        boolean   nullQueryStrand = false;
        boolean nullSubjectStrand = false;

        // Start with index 0 value (arbitrarily)
        Strand qStrand = subs[0].getQueryStrand();
        Strand sStrand = subs[0].getSubjectStrand();

        int qStart = subs[0].getQueryStart();
        int qEnd   = subs[0].getQueryEnd();
        int sStart = subs[0].getSubjectStart();
        int sEnd   = subs[0].getSubjectEnd();

        if (qStrand == null)
            nullQueryStrand = true;
        if (sStrand == null)
            nullSubjectStrand = true;

        // Compare all other values
        for (int i = subs.length; --i > 0;)
        {
            Strand qS = subs[i].getQueryStrand();
            Strand sS = subs[i].getSubjectStrand();

            if (qS == null)
                nullQueryStrand = true;
            if (sS == null)
                nullSubjectStrand = true;

            if (qS != qStrand)
                mixQueryStrand = true;
            if (sS != sStrand)
                mixSubjectStrand = true;

            qStart = Math.min(qStart, subs[i].getQueryStart());
            qEnd   = Math.max(qEnd,   subs[i].getQueryEnd());

            sStart = Math.min(sStart, subs[i].getSubjectStart());
            sEnd   = Math.max(sEnd,   subs[i].getSubjectEnd());
        }

        // Note any mixed strand hits as unknown strand
        if (mixQueryStrand)
            qStrand = StrandedFeature.UNKNOWN;
        if (mixSubjectStrand)
            sStrand = StrandedFeature.UNKNOWN;

        // Any null strands from protein sequences
        if (nullQueryStrand)
            qStrand = null;
        if (nullSubjectStrand)
            sStrand = null;

        String subjectID = (String) hitData.get("subjectId");

        return new SimpleSeqSimilaritySearchHit(sc, ev, pv,
                                                qStart, qEnd, qStrand,
                                                sStart, sEnd, sStrand,
                                                subjectID,
                                                AnnotationFactory.makeAnnotation(hitData),
                                                subHits);
    }

    /**
     * <code>makeSubHit</code> creates a new sub-hit.
     *
     * @return a <code>SeqSimilaritySearchSubHit</code>.
     *
     * @exception BioException if an error occurs.
     */
    private SeqSimilaritySearchSubHit makeSubHit() throws BioException
    {
        // Try to get a valid TokenParser
        if (tokenParser == null)
        {
            String identifier;

            // Try explicit sequence type first
            if (subHitData.containsKey("subjectSequenceType"))
                identifier = (String) subHitData.get("subjectSequenceType");
            // Otherwise try to resolve from the program name (only
            // works for Blast)
            else if (resultPreAnnotation.containsKey("program"))
                identifier = (String) resultPreAnnotation.get("program");
            else
                throw new BioException("Failed to determine sequence type");

            FiniteAlphabet alpha = AlphabetResolver.resolveAlphabet(identifier);
            tokenParser = alpha.getTokenization("token");
        }

        // BLASTP output has the strands set null (protein sequences)
        Strand qStrand = null;
        Strand sStrand = null;

        // Override where an explicit strand is given (FASTA DNA,
        // BLASTN)
        if (subHitData.containsKey("queryStrand"))
            if (subHitData.get("queryStrand").equals("plus"))
                qStrand = StrandedFeature.POSITIVE;
            else
                qStrand = StrandedFeature.NEGATIVE;

        if (subHitData.containsKey("subjectStrand"))
            if (subHitData.get("subjectStrand").equals("plus"))
                sStrand = StrandedFeature.POSITIVE;
            else
                sStrand = StrandedFeature.NEGATIVE;

        // Override where a frame is given as this contains strand
        // information (BLASTX for query, TBLASTN for hit, TBLASTX for
        // both)
        if (subHitData.containsKey("queryFrame"))
            if (((String) subHitData.get("queryFrame")).startsWith("plus"))
                qStrand = StrandedFeature.POSITIVE;
            else
                qStrand = StrandedFeature.NEGATIVE;

        if (subHitData.containsKey("subjectFrame"))
            if (((String) subHitData.get("subjectFrame")).startsWith("plus"))
                sStrand = StrandedFeature.POSITIVE;
            else
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

        // Get scores
        double sc = Double.NaN;
        double ev = Double.NaN;
        double pv = Double.NaN;

        if (subHitData.containsKey("score"))
            sc = Double.parseDouble((String) subHitData.get("score"));

        if (subHitData.containsKey("expectValue"))
        {
            String val = (String) subHitData.get("expectValue");
            // Blast sometimes uses invalid formatting such as 'e-156'
            // rather than '1e-156'
            if (val.startsWith("e"))
                ev = Double.parseDouble("1" + val);
            else
                ev = Double.parseDouble(val);
        }

        if (subHitData.containsKey("pValue"))
            pv = Double.parseDouble((String) subHitData.get("pValue"));

        Map labelMap = new SmallMap();

        // Note that the following is removing the raw sequences
        StringBuffer tokenBuffer = new StringBuffer(1024);
        tokenBuffer.append((String) subHitData.remove("querySequence"));
        labelMap.put(SeqSimilaritySearchSubHit.QUERY_LABEL,
                     new SimpleSymbolList(tokenParser, tokenBuffer.substring(0)));

        tokenBuffer = new StringBuffer(1024);
        tokenBuffer.append((String) subHitData.remove("subjectSequence"));
        labelMap.put(hitData.get("subjectId"),
                     new SimpleSymbolList(tokenParser, tokenBuffer.substring(0)));

        return new SimpleSeqSimilaritySearchSubHit(sc, ev, pv,
                                                   qStart, qEnd, qStrand,
                                                   sStart, sEnd, sStrand,
                                                   new SimpleAlignment(labelMap),
                                                   AnnotationFactory.makeAnnotation(subHitData));
    }
}
