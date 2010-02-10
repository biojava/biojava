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
import java.util.ArrayList;
import java.util.Iterator;

import junit.framework.TestCase;

import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojavax.Namespace;
import org.biojavax.Note;
import org.biojavax.RichObjectFactory;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.io.GenbankFormat.Terms;

/**
 * Test for EMBL KW parsing bug #2255
 * @author Mark Schreiber
 * @author Richard Holland
 */
public class Bug2255Test extends TestCase{
    private EMBLFormat gbFormat;

    /**
     * @see junit.framework.TestCase#setUp()
     */
    protected void setUp() {
        this.gbFormat = new EMBLFormat();
    }

    protected void tearDown() throws Exception {
        gbFormat = null;
    }    

    
    /** Creates a new instance of Bug2250Test */
    public Bug2255Test() {
    }
    
    /**
     * Mainly tests the parsing of the locus line.
     */
    public void testBug2255(){
        String filename = "/X56734.embl";
        RichSequence seq = readFile(filename);
        
        ArrayList kws = new ArrayList();
        for (Iterator i = seq.getNoteSet().iterator(); i.hasNext(); ) {
            Note n = (Note)i.next();
            if (n.getTerm().equals(Terms.getKeywordTerm())) kws.add(n.getValue());
        }
        assertEquals(kws.get(0),"Hello");
        assertEquals(kws.get(1),"World");
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
