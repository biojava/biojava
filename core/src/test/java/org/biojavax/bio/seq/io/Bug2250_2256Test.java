/*
 * Bug2250Test.java
 *
 * Created on March 26, 2007, 4:44 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package org.biojavax.bio.seq.io;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;

import junit.framework.TestCase;

import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojavax.Namespace;
import org.biojavax.Note;
import org.biojavax.RichObjectFactory;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.seq.RichLocation;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.io.GenbankFormat.Terms;

/**
 * Test for locus line parsing bug #2250 and #2256
 * @author Mark Schreiber
 * @author Richard Holland
 */
public class Bug2250_2256Test extends TestCase{
    private GenbankFormat gbFormat;

    /**
     * @see junit.framework.TestCase#setUp()
     */
    protected void setUp() {
        this.gbFormat = new GenbankFormat();
    }

    protected void tearDown() throws Exception {
        gbFormat = null;
    }    

    
    /** Creates a new instance of Bug2250_2256Test */
    public Bug2250_2256Test() {
    }
    
    /**
     * Mainly tests the parsing of the locus line.
     */
    public void testBug2250_2256(){
        String filename = "/AL121964.gb";
        RichSequence seq = readFile(filename);
        
        assertEquals("AL121964",seq.getAccession());
        assertEquals("HSJ154G14", seq.getName());
        assertFalse(seq.getCircular());
        assertEquals("PRI", seq.getDivision());
        
        String udat = null;
        String molType = seq.getAlphabet().getName();
        for (Iterator i = seq.getNoteSet().iterator(); i.hasNext(); ) {
            Note n = (Note)i.next();
            if (n.getTerm().equals(Terms.getDateUpdatedTerm())) udat=n.getValue();
            else if (n.getTerm().equals(Terms.getMolTypeTerm())) molType=n.getValue();
        }
        assertNotNull(udat);
        assertEquals("18-MAY-2005", udat);
        assertNotNull(molType);
        assertEquals("DNA", molType);
        assertEquals(1,
        		((RichLocation)((RichFeature)seq.getFeatureSet().iterator().next()).getLocation()).getRank());
    }
    
    /**
     * Read a genbank file, return a RichSequence
     * @param filename name of file to read
     * @return a RichSequence instance
     */
    private RichSequence readFile(String filename) {
	InputStream inStream = this.getClass().getResourceAsStream(filename);
        BufferedReader br = new BufferedReader(new InputStreamReader(inStream));
        SymbolTokenization tokenization = RichSequence.IOTools.getDNAParser();
        Namespace namespace = RichObjectFactory.getDefaultNamespace();
        SimpleRichSequenceBuilder builder = new SimpleRichSequenceBuilder();
        RichSequence sequence = null;
        try {
            this.gbFormat.readRichSequence(br, tokenization, builder, namespace);
            sequence = builder.makeRichSequence();
        } catch (Exception e) {
            e.printStackTrace();
            fail("Unexpected exception: "+e);
        }
		return sequence;
	}    
}
