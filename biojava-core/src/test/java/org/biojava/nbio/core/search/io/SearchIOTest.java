package org.biojava.nbio.core.search.io;

import java.io.File;
import java.net.URL;
import org.biojava.nbio.core.search.io.blast.BlastXMLQuery;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Designed by Paolo Pavan.
 * You may want to find my contacts on Github and LinkedIn for code info 
 * or discuss major changes.
 * https://github.com/paolopavan
 * 
 * @author Paolo Pavan
 */

public class SearchIOTest {
    
    public SearchIOTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    /**
     * Constructor test with GuessFactory
     */
    @Test
    public void testConstructorWithFactoryGuess() {
        String resource = "/org/biojava/nbio/core/search/io/blast/test.two-query.blasttxt";
        URL resourceURL = getClass().getResource(resource);
        File file = new File(resourceURL.getFile());
        
        final SearchIO instance;
        try {
            instance = new SearchIO(file);
        } catch (Exception e) {
            fail("test failed:\n"+e.getMessage());
        }
    }
    /**
     * Constructor test specifying Factory
     */
    @Test
    public void testConstructorWithoutFactoryGuess() {
        String resource = "/org/biojava/nbio/core/search/io/blast/testBlastReport.xml";
        URL resourceURL = getClass().getResource(resource);
        File file = new File(resourceURL.getFile());
        
        ResultFactory blastResultFactory = new BlastXMLQuery();
        final SearchIO instance;
        try {
            instance = new SearchIO(file, blastResultFactory);
        } catch (Exception e) {
            fail("test failed:\n"+e.getMessage());
        }
    }
    /**
     * Constructor test specifying Factory and using a evalue threshold filter
     */
    @Test
    public void testConstructorWithEvalueHspFilter() {
        //
        String resource = "/org/biojava/nbio/core/search/io/blast/testBlastReport.xml";
        URL resourceURL = getClass().getResource(resource);
        File file = new File(resourceURL.getFile());
        
        ResultFactory blastResultFactory = new BlastXMLQuery();
        final SearchIO instance;
        try {
            instance = new SearchIO(file, blastResultFactory, 10e-10);
        } catch (Exception e) {
            fail("test failed:\n"+e.getMessage());
        }
    }
}
