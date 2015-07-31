/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.nbio.core.search.io;

import java.io.File;
import java.net.URL;
import java.util.Iterator;
import org.biojava.nbio.core.search.io.blast.BlastXMLQuery;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author pavanpa
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
    
    @Test
    public void testConstructorWithFactoryGuess() {
        System.out.println("Constructor test with GuessFactory");
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
    
    @Test
    public void testConstructorWithoutFactoryGuess() {
        System.out.println("Constructor test specifying Factory");
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
    
    @Test
    public void testConstructorWithEvalueHspFilter() {
        System.out.println("Constructor test specifying Factory");
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
