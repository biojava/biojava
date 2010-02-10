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

import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.search.SeqSimilaritySearchHit;
import org.biojava.bio.search.SeqSimilaritySearchResult;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.biojava.bio.seq.db.DummySequenceDB;
import org.biojava.bio.seq.db.DummySequenceDBInstallation;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.db.SequenceDBInstallation;

/**
 * <code>SSBindCase</code> is a base class for tests of object
 * bindings for Blast-like SAX events.
 *
 * @author Keith James
 * @since 1.2
 */
public class SSBindCase extends TestCase
{
    protected SequenceDB             queryDB;
    protected SequenceDBInstallation dbInstallation;

    protected SeqSimilarityAdapter adapter;
    protected InputStream          searchStream;
    protected List                 searchResults;

    protected double   topHitScore;
    protected String   topHitSeqID;
    protected int     topHitQStart;
    protected int       topHitQEnd;
    protected Strand topHitQStrand;
    protected int     topHitSStart;
    protected int       topHitSEnd;
    protected Strand topHitSStrand;

    protected double   botHitScore;
    protected String   botHitSeqID;
    protected int     botHitQStart;
    protected int       botHitQEnd;
    protected Strand botHitQStrand;
    protected int     botHitSStart;
    protected int       botHitSEnd;
    protected Strand botHitSStrand;

    public SSBindCase(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {
    	if (this.getClass().getName().equals("org.biojava.bio.program.ssbind.SSBindCase")) return;
        queryDB        = new DummySequenceDB("query");
        dbInstallation = new DummySequenceDBInstallation();
        searchResults  = new ArrayList();

        // Set builder to build into a List
        BlastLikeSearchBuilder builder =
            new BlastLikeSearchBuilder(searchResults);

        // Set the holder for query sequences and databases
        builder.setQuerySeqHolder(queryDB);
        builder.setSubjectDBInstallation(dbInstallation);

        // Adapter from SAX -> search result construction interface
        adapter = new SeqSimilarityAdapter();

        // Set the handler which will instantiate result objects
        adapter.setSearchContentHandler(builder);
    }

    protected void tearDown() throws Exception
    {
    	if (this.getClass().getName().equals("org.biojava.bio.program.ssbind.SSBindCase")) return;
        searchStream.close();
        
        searchStream = null;
        queryDB = null;
        dbInstallation = null;
        adapter = null;
        searchResults = null;
        topHitSeqID = null;
        botHitSeqID = null;
    }
    
    public void testResultCount()
    {
    	if (this.getClass().getName().equals("org.biojava.bio.program.ssbind.SSBindCase")) return;
        assertEquals(1, searchResults.size());
    }

    public void testResultGetQuerySequence() throws Exception
    {
    	if (this.getClass().getName().equals("org.biojava.bio.program.ssbind.SSBindCase")) return;
        SeqSimilaritySearchResult result =
            (SeqSimilaritySearchResult) searchResults.get(0);

        assertEquals(queryDB.getSequence(""), result.getQuerySequence());
    }

    public void testResultGetSequenceDB()
    {
    	if (this.getClass().getName().equals("org.biojava.bio.program.ssbind.SSBindCase")) return;
        SeqSimilaritySearchResult result =
            (SeqSimilaritySearchResult) searchResults.get(0);

        assertEquals(dbInstallation.getSequenceDB(""), result.getSequenceDB());
    }

    public void testResultGetAnnotation()
    {
    	if (this.getClass().getName().equals("org.biojava.bio.program.ssbind.SSBindCase")) return;
        SeqSimilaritySearchResult result =
            (SeqSimilaritySearchResult) searchResults.get(0);
        int nrAnnotations = result.getAnnotation().keys().size();
        // there are either 4 or 5 annotations now with the new support for query length
        assertTrue( nrAnnotations > 3);
        assertTrue( nrAnnotations < 6);
    }

    public void testTopHit()
    {
    	if (this.getClass().getName().equals("org.biojava.bio.program.ssbind.SSBindCase")) return;
        SeqSimilaritySearchResult result =
            (SeqSimilaritySearchResult) searchResults.get(0);

        SeqSimilaritySearchHit hit =
            (SeqSimilaritySearchHit) result.getHits().get(0);

        assertEquals(topHitScore, hit.getScore(), 0.0);
        assertEquals(topHitSeqID, hit.getSubjectID());

        assertEquals(topHitQStart, hit.getQueryStart());
        assertEquals(topHitQEnd,   hit.getQueryEnd());
        assertSame(topHitQStrand,  hit.getQueryStrand());

        assertEquals(topHitSStart, hit.getSubjectStart());
        assertEquals(topHitSEnd,   hit.getSubjectEnd());
        assertSame(topHitSStrand,  hit.getSubjectStrand());
    }

    public void testBotHit()
    {
    	if (this.getClass().getName().equals("org.biojava.bio.program.ssbind.SSBindCase")) return;
        SeqSimilaritySearchResult result =
            (SeqSimilaritySearchResult) searchResults.get(0);

        List hits = result.getHits();

        SeqSimilaritySearchHit hit =
            (SeqSimilaritySearchHit) hits.get(hits.size() - 1);

        assertEquals(botHitScore, hit.getScore(), 0.0);
        assertEquals(botHitSeqID, hit.getSubjectID());

        assertEquals(botHitQStart, hit.getQueryStart());
        assertEquals(botHitQEnd,   hit.getQueryEnd());
        assertSame(botHitQStrand,  hit.getQueryStrand());

        assertEquals(botHitSStart, hit.getSubjectStart());
        assertEquals(botHitSEnd,   hit.getSubjectEnd());
        assertSame(botHitSStrand,  hit.getSubjectStrand());
    }

    protected void setTopHitValues(double score, String id,
                                   int qStart, int qEnd, Strand qStrand,
                                   int sStart, int sEnd, Strand sStrand)
    {
    	if (this.getClass().getName().equals("org.biojava.bio.program.ssbind.SSBindCase")) return;
        topHitScore   = score;
        topHitSeqID   = id;
        topHitQStart  = qStart;
        topHitQEnd    = qEnd;
        topHitQStrand = qStrand;
        topHitSStart  = sStart;
        topHitSEnd    = sEnd;
        topHitSStrand = sStrand;
    }

    protected void setBotHitValues(double score, String id,
                                   int qStart, int qEnd, Strand qStrand,
                                   int sStart, int sEnd, Strand sStrand)
    {
    	if (this.getClass().getName().equals("org.biojava.bio.program.ssbind.SSBindCase")) return;
        botHitScore   = score;
        botHitSeqID   = id;
        botHitQStart  = qStart;
        botHitQEnd    = qEnd;
        botHitQStrand = qStrand;
        botHitSStart  = sStart;
        botHitSEnd    = sEnd;
        botHitSStrand = sStrand;
    }
}
