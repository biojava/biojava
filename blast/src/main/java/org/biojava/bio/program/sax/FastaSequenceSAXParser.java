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
package org.biojava.bio.program.sax;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.StringTokenizer;

import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;

/**
 * A SAX2 parser for dealing with multiple sequences in
 * FASTA format.
 *
 * For example:
 * <pre>
 * >Seq1
 * GATCGATCGTAGCTAGATGCTAGCATGCTAGCTGACTGATCGATCGTAGCTAGCTAGCTGACTG
 * >Seq2
 * GATCGATCGTAGCTAGATGCTAGCATGCTAGCTGACTGATCGATCGTAGCTAGCTAGCTGACTG
 * </pre>
 * <p>
 *
 * Copyright &copy; 2000,2001 Cambridge Antibody Technology.
 
 * <p>
 * Primary author -<ul>
 * <li>Simon Brocklehurst (CAT)
 * </ul>
 * Other authors  -<ul>
 * <li>Neil Benn          (CAT)
 * <li>Lawrence Bower     (CAT)
 * <li>Derek Crockford    (CAT)
 * <li>Tim Dilks          (CAT)
 * <li>Colin Hardman      (CAT)
 * <li>Stuart Johnston    (CAT)
 *</ul>
 *
 * @author Cambridge Antibody Technology (CAT)
 * @author Greg Cox
 * @version 1.0
 *
 */
public class FastaSequenceSAXParser extends AbstractNativeAppSAXParser {



    private AttributesImpl          oAtts     = new AttributesImpl();
    private QName                   oAttQName = new QName(this);
    private char[]                  aoChars;

    private StringBuffer            oSeqName  = new StringBuffer();
    private StringBuffer            oSeq      = new StringBuffer();
    private boolean                 tOnFirst  = true;

    private static final int        STARTUP            = 0;
    private static final int        IN_STREAM          = 1;


    /**
     * Initialises internal state
     * Sets namespace prefix to "biojava"
     */
    public FastaSequenceSAXParser() {
	iState = STARTUP;
	this.setNamespacePrefix("biojava");
    }

    /**
     * Describe 'parse' method here.
     *
     * @param poSource	 -
     */
    public void parse(InputSource poSource )
	throws IOException,SAXException {

	BufferedReader            oContents;
	String                    oLine = null;

	//Use method form superclass
	oContents = this.getContentStream(poSource);

	// loop over file
	try {
	    // loop over file
	    oLine = oContents.readLine();
	    while (oLine != null) {
		//System.out.println(oLine);
		this.interpret(oContents,oLine);
		oLine = oContents.readLine();
	    } // end while
	} catch (java.io.IOException x) {
	    System.out.println(x.getMessage());
	    System.out.println("Stream read interupted");
	} // end try/catch

	//at end of stream...
	//do final sequence
	this.emitSequence();

	this.endElement(new QName(this,
				  this.prefix("SequenceCollection")));
	oContents.close();

    }

    /**
     * Describe <code>interpret</code> method here.
     *
     * @param poContents a <code>BufferedReader</code> value
     * @param poLine a <code>String</code> value
     * @exception SAXException if an error occurs
     */
    private void interpret(BufferedReader poContents, String poLine)
	throws SAXException {


	if (iState == STARTUP) {
	    oAtts.clear();
	    this.startElement(
	      new QName(this,
			this.prefix("SequenceCollection")),
				  (Attributes)oAtts);
	    this.changeState(IN_STREAM);
	}

	if (iState == IN_STREAM) {
	    //look for the start of first record i.e.a header
	    if ( poLine.startsWith(">") ) {
		if (!tOnFirst) {
		    this.emitSequence();
		}
		this.parseHeaderLine(poLine);
		oSeq.setLength(0);
		return;
	    } else {
		this.appendSequence(poLine);
	    }

	}
    }
    /**
     * Parse the header part of a record i.e. >myName, and
     * emit messages.
     *
     * @param poLine a <code>String</code> value
     */
    private void parseHeaderLine(String poLine) {
	oSeqName.setLength(0);
	oSeqName.append(poLine.substring(1));

	//flip flag to begin emitting sequence elements
	tOnFirst = false;
	//System.out.println(oSeqName);
    }
    /**
     * Builds up sequence data - NB white space is
     * removed.
     *
     * @param poLine a <code>String</code> value
     */
    private void appendSequence(String poLine) {
	StringTokenizer oSt = new StringTokenizer(poLine,"\n\t\r ");
        while (oSt.hasMoreTokens()) {
          oSeq.append(oSt.nextToken());
        }
    }
    /**
     * Describe <code>emitSequence</code> method here.
     *
     */
    private void emitSequence() throws SAXException {
	    oAtts.clear();

	    oAttQName.setQName("sequenceName");
	    oAtts.addAttribute(oAttQName.getURI(),
			   oAttQName.getLocalName(),
			   oAttQName.getQName(),
			   "CDATA",oSeqName.substring(0));

	    this.startElement(
	      new QName(this,
			this.prefix("Sequence")),
				  (Attributes)oAtts);

	    aoChars = oSeq.substring(0).toCharArray();
	    this.characters(aoChars,0,aoChars.length);
	    this.endElement(new QName(this,this.prefix("Sequence")));

    }
}
