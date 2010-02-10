/*
 * Bug2249_2248Test.java
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
import java.io.OutputStream;
import java.util.Iterator;

import junit.framework.TestCase;

import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojavax.Namespace;
import org.biojavax.Note;
import org.biojavax.RichObjectFactory;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.io.GenbankFormat;
import org.biojavax.bio.seq.io.SimpleRichSequenceBuilder;
import org.biojavax.bio.seq.io.GenbankFormat.Terms;

/**
 * Test for Genpept->Uniprot docref parsing bugs #2249 and #2248.
 * @author Mark Schreiber
 * @author Richard Holland
 */
public class Bug2249_2248Test extends TestCase{
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

    
    /** Creates a new instance of Bug2249_2248Test */
    public Bug2249_2248Test() {
    }
    
    /**
     * Mainly tests the parsing of the locus line.
     */
    public void testBug2249_2248(){
        String filename = "/AAX56332.gp";
        RichSequence seq = readFile(filename);
        
        assertEquals("AAX56332",seq.getAccession());
        assertFalse(seq.getCircular());
        assertEquals("PRI", seq.getDivision());
        
        String lengthType = null;
        String molType = null;
        for (Iterator i = seq.getNoteSet().iterator(); i.hasNext(); ) {
            Note n = (Note)i.next();
            if (n.getTerm().equals(Terms.getMolTypeTerm())) molType=n.getValue();
        }
        assertNull(molType);

        // Make sure we can convert.
        OutputStream nos = new OutputStream() {
              public void write( int b ) {}
              public void write( byte b[] ) {}
              public void write( byte b[], int off, int len ) {}
              public void flush() {}
              public void close() {}
        };
        try {
            RichSequence.IOTools.writeUniProt(nos, seq, seq.getNamespace());
        } catch (Exception e) {
            e.printStackTrace();
            fail("Unexpected exception: "+e);
        }
    }
    
    /**
     * Read a genpept file, return a RichSequence
     * @param filename name of file to read
     * @return a RichSequence instance
     */
    private RichSequence readFile(String filename) {
	InputStream inStream = this.getClass().getResourceAsStream(filename);
        BufferedReader br = new BufferedReader(new InputStreamReader(inStream));
        SymbolTokenization tokenization = RichSequence.IOTools.getProteinParser();
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
