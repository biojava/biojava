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

package org.biojava.bio.program.sax.blastxml;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.program.ssbind.BlastLikeSearchBuilder;
import org.biojava.bio.program.ssbind.SeqSimilarityAdapter;
import org.biojava.bio.search.SearchContentHandler;
import org.biojava.bio.search.SeqSimilaritySearchResult;
import org.biojava.bio.seq.db.DummySequenceDB;
import org.biojava.bio.seq.db.DummySequenceDBInstallation;
import org.xml.sax.InputSource;

public class BlastXMLTest
    extends TestCase
{
    private List results;

    protected void setUp()
        throws Exception
    {
        // get test input file
        String resName = "org/biojava/bio/program/sax/blastxml/input.xml";
        InputStream resStream = getClass().getClassLoader().getResourceAsStream(resName);
        assert resStream != null
                : "Resource " + resName + " could not be located";
        InputSource is = new InputSource(resStream);

        //make a BlastLikeSAXParser
        BlastXMLParserFacade parser = new BlastXMLParserFacade();

        //make the SAX event adapter that will pass events to a Handler.
        SeqSimilarityAdapter adapter = new SeqSimilarityAdapter();
  
        //set the parsers SAX event adapter
        parser.setContentHandler(adapter);
  
        //The list to hold the SeqSimilaritySearchResults
        results = new ArrayList();
  
        //create the SearchContentHandler that will build SeqSimilaritySearchResults
        //in the results List
        SearchContentHandler builder = new BlastLikeSearchBuilder(results,
            new DummySequenceDB("queries"), new DummySequenceDBInstallation());
  
        //register builder with adapter
        adapter.setSearchContentHandler(builder);
  
        //parse the file, after this the result List will be populated with
        //SeqSimilaritySearchResults
        parser.parse(is);
    }

    private static int EXPECTED_QUERY_COUNT = 1;
    private static int EXPECTED_HIT_COUNT = 5;

    public void testHitCount()
    {
        // check that there are right numbers of queries
        assertEquals(EXPECTED_QUERY_COUNT, results.size());

        // check that number of hits is correct
        SeqSimilaritySearchResult result =
            (SeqSimilaritySearchResult) results.get(0);
        assertEquals(EXPECTED_HIT_COUNT, result.getHits().size());
    }
}

