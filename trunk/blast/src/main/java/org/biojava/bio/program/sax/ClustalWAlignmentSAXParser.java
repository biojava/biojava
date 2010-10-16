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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;

/**
 * A SAX2 parser for dealing with a multiple sequence
 * alignment as produced by ClustalW outputing .aln format.
 * For example,
 * <pre>
  K1C0_XENLA/125-441      DKVHALETANTELERKIKEWYEKQRPGSSSGDGAKDYSKYYT
  K1C4_XENLA/81-396       EKVRALEAANADLELKIREWYEKQK-GSGIGAGSKDFSKYFE
  K1C5_XENLA/73-384       DRVRSLEQANHELELKIREYLDKK-----AAVGSLDYSGYYN
  keratin15               DKVRALEEANADLEVKIHDWYQKQTP----ASPECDYSQYFK

  K1C0_XENLA/125-441      -----AKFLLQNDNARLAADDFKMKFEN--------------
  K1C4_XENLA/81-396       -----SRVVLQIDNAKLAADDFRLKFEN--------------
  K1C5_XENLA/73-384       -----TRLVLSIDNAKLAADDFKIKYES--------------
  keratin15               -----SRVILEIDNARLAADDFRLKYEN--------------
 * </pre>
 * <p>
 * Please note, this parser reads the whole alignment in to
 * core memory and thus does not scale to work with very large
 * alignments on low-end hardware.
 * <p>
 * Please also note that this class has not been tested with
 * many version of CLUSTAL W.
 *
 * Copyright &copy; 2000,2001 Cambridge Antibody Technology.
 * 
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
public class ClustalWAlignmentSAXParser extends AbstractNativeAppSAXParser {



    private AttributesImpl          oAtts      = new AttributesImpl();
    private QName                   oAttQName  = new QName(this);
    private char[]                  aoChars;

    private String                  oSeqName;
    private String                  oTmpSeq;
    private StringBuffer            oSeq         = new StringBuffer();
    private HashMap                 oAlignment   = new HashMap();
    private ArrayList               oSeqNameList = new ArrayList();

    private static final int        STARTUP            = 0;
    private static final int        IN_STREAM          = 1;


    /**
     * Initialises internal state
     * Sets namespace prefix to "biojava"
     */
    public ClustalWAlignmentSAXParser() {
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
	    System.out.println("Stream read interrupted");
	} // end try/catch

	//at end of stream...

	//at this point, alignment is parsed, now cycle through
	//and emit elements
	for (int i = 0; i < oSeqNameList.size(); i++) {
	    oSeqName = (String) oSeqNameList.get(i);
	    this.emitSequence(oSeqName,
			      (String) oAlignment.get(oSeqName));


	}

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

	    if (this.lineIsRelevant(poLine)) {
		//build aligment in memory
		this.appendToAlignment(poLine);
	    }

	}
    }
    /**
     * Parse a relevant line, and add to alignment
     *
     * @param poLine a <code>String</code> value
     */
    private void appendToAlignment(String poLine) {
	//System.out.println(poLine);
	StringTokenizer oSt = new StringTokenizer(poLine,"\n\t\r ");

	//First token is sequence name
	oSeqName = oSt.nextToken();
	//System.out.println(oSeqName);

	oSeq.setLength(0);
        while (oSt.hasMoreTokens()) {
          oSeq.append(oSt.nextToken());
        }
	//System.out.println(oSeq);

	//At this point, have name of sequence, and a segment of the sequence

	//Update object...

	if (oAlignment.get(oSeqName) == null) {
	    //Here if on first occurence of this sequence
	    //Add to alignment
	    oAlignment.put(oSeqName,oSeq.substring(0));
	    //maintain ordered list of sequence names
	    oSeqNameList.add(oSeqName);
	} else {
	    //Here if building up an existing sequence
	    oTmpSeq = (String) oAlignment.get(oSeqName);
	    oAlignment.put(oSeqName,oTmpSeq.concat(oSeq.substring(0)));
	}
    }

    /**
     * Only interested in lines that are part of the alignment.
     * Returns true if line is in alignment, false if not.
     *
     * @param poLine a <code>String</code> value
     * @return a <code>boolean</code> value
     */
    private boolean lineIsRelevant(String poLine) {

	//blank lines not relevant
	//lines that starts with a space  not relevant (consensus line)
	//lines that start with "CLUSTAL W (" not relevant (title)

	if ( (poLine.trim().equals("")) ||
	     (poLine.startsWith(" ")) ||
	     (poLine.startsWith("CLUSTAL W (")) ) {

	    //System.out.println("Irrelevant|"+poLine+"|");
	    return false;
	}

	//if here,line is part of alignment, so return true
	return true;
    }
    /**
     * Emit a sequence element
     *
     * @param poSequenceName a <code>String</code> value
     * @param poSequence a <code>String</code> value
     * @exception SAXException if an error occurs
     */
    private void emitSequence(String poSequenceName, String poSequence) throws SAXException {
	    oAtts.clear();

	    oAttQName.setQName("sequenceName");
	    oAtts.addAttribute(oAttQName.getURI(),
			   oAttQName.getLocalName(),
			   oAttQName.getQName(),
			   "CDATA",poSequenceName);

	    this.startElement(
	      new QName(this,
			this.prefix("Sequence")),
				  (Attributes)oAtts);

	    aoChars = poSequence.toCharArray();
	    this.characters(aoChars,0,aoChars.length);
	    this.endElement(new QName(this,this.prefix("Sequence")));

    }
}
