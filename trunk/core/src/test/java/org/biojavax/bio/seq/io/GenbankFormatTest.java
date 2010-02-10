package org.biojavax.bio.seq.io;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.StringReader;
import java.util.Iterator;

import junit.framework.TestCase;

import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojavax.Namespace;
import org.biojavax.Note;
import org.biojavax.RichObjectFactory;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.io.GenbankFormat.Terms;
import org.biojavax.bio.taxa.NCBITaxon;

/**
 * Tests for GenbankFormat. Ain't parsing fun?
 * @author Bubba Puryear
 * @author George Waldon
 */
public class GenbankFormatTest extends TestCase {
    private GenbankFormat gbFormat;

    /**
     * @see junit.framework.TestCase#setUp()
     */
    protected void setUp() {
        this.gbFormat = new GenbankFormat();
    }
		
		public void testGenbankParsingWithBondFeatures() {
				readProteinFile("/BondFeature.gb");
		}

    public void testGenbankParsing_oldStyleFile() {
		RichSequence sequence = readDNAFile("/NoAccession.gb");
        assertEquals("NoAccess", sequence.getName());
        assertTrue(sequence.getCircular());
        assertEquals(null, sequence.getDescription());
        assertEquals(null, sequence.getDivision());
        assertEquals(null, sequence.getTaxon());
        assertEquals("NoAccess", sequence.getURN());
        assertEquals(0, sequence.getVersion());
        String stranded = null;
        String udat = null;
        String molType = sequence.getAlphabet().getName();
        for (Iterator i = sequence.getNoteSet().iterator(); i.hasNext(); ) {
            Note n = (Note)i.next();
            if (n.getTerm().equals(Terms.getStrandedTerm())) stranded=n.getValue();
            else if (n.getTerm().equals(Terms.getDateUpdatedTerm())) udat=n.getValue();
            else if (n.getTerm().equals(Terms.getMolTypeTerm())) molType=n.getValue();
        }
        assertNull(stranded);
        assertNotNull(udat);
        assertEquals("30-JUN-2006", udat);
        assertNotNull(molType);
        assertEquals("DNA", molType);
    }


    public void testGenbankParsing_contemporaryApp() {
		RichSequence sequence = readDNAFile("/AY069118.gb");
        assertEquals("AY069118", sequence.getName());
        assertFalse(sequence.getCircular());
        assertEquals("Drosophila melanogaster GH13089 full length cDNA.", sequence.getDescription());
        assertEquals("INV", sequence.getDivision());
        NCBITaxon taxon = sequence.getTaxon();
        assertNotNull(taxon);
        assertEquals("Drosophila melanogaster", taxon.getDisplayName());
        assertEquals("AY069118", sequence.getURN());
        assertEquals(1, sequence.getVersion());
        String stranded = null;
        String udat = null;
        String molType = sequence.getAlphabet().getName();
        for (Iterator i = sequence.getNoteSet().iterator(); i.hasNext(); ) {
            Note n = (Note)i.next();
            if (n.getTerm().equals(Terms.getStrandedTerm())) stranded=n.getValue();
            else if (n.getTerm().equals(Terms.getDateUpdatedTerm())) udat=n.getValue();
            else if (n.getTerm().equals(Terms.getMolTypeTerm())) molType=n.getValue();
        }
        assertNull(stranded);
        assertNotNull(udat);
        assertEquals("17-DEC-2001", udat);
        assertNotNull(molType);
        assertEquals("mRNA", molType);
    }


    public void testGenbankWithNoAccession() {
    	RichSequence sequence = readDNAFile("/NoAccession.gb");
        assertNotNull(sequence);
        assertEquals("NoAccess", sequence.getAccession());
    }

    public void testCanReadWhatIsWritten() {
    	// Read a genbank file
    	RichSequence sequence = readDNAFile("/AY069118.gb");
        assertNotNull(sequence);

        // Write the file to an in-memory buffer
        OutputStream output = new ByteArrayOutputStream();
		RichSequenceFormat genbank = new GenbankFormat();
		RichStreamWriter seqsOut = new RichStreamWriter(output, genbank);
		SequenceIterator seqIterator = new RichSequence.IOTools.SingleRichSeqIterator(sequence);
		try {
			seqsOut.writeStream(seqIterator, null);
		} catch (IOException e) {
        	fail("Unexpected exception: "+e);
		}

		// Re-read the generated output
		String newContent = output.toString();
        SymbolTokenization dna = RichSequence.IOTools.getDNAParser();
        Namespace defaultNs = RichObjectFactory.getDefaultNamespace();
		BufferedReader input = new BufferedReader(new StringReader(newContent));
		RichSequence rereadSeq = null;
        try {
            RichStreamReader reader = new RichStreamReader(input, new GenbankFormat(), dna, RichSequenceBuilderFactory.FACTORY, defaultNs);
            rereadSeq = reader.nextRichSequence();
        } catch (Exception e) {
        	e.printStackTrace();
        	fail("Unexpected exception: "+e);
        }
        assertNotNull(rereadSeq);
        assertEquals(sequence.getAccession(), rereadSeq.getAccession());
        assertEquals(sequence.getName(), rereadSeq.getName());
        assertEquals(sequence.seqString(), rereadSeq.seqString());
    }

    /** Test whether the parser reads minimal sequences. The sequence prototype
     * was generated by writing a sequence read in fasta format 
     * (">testempty no sequence") under the tested format.
     */
    public void testReadEmptySequence() {
        RichSequence sequence = readDNAFile("/empty_genbank.gb");
        assertNotNull(sequence);
        assertEquals(sequence.getName(), "testempty");
        assertEquals(sequence.getAccession(), "");
        assertEquals(sequence.getVersion(), 0);
        assertEquals(sequence.getDescription(), "no sequence");
        assertEquals(sequence.getInternalSymbolList().length(), 0);
    }

    /**
     * Read a genbank file, return a RichSequence
     * @param filename name of file to read
     * @return a RichSequence instance
     */
    private RichSequence readDNAFile(String filename) {
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
	
	private RichSequence readProteinFile(String filename) {
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
