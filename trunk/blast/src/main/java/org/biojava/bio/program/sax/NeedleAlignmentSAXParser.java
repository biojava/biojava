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
 * A SAX2 parser for dealing with a pairwise sequence
 * alignment as produced by Needle (EMBOSS) outputing .needle format.
 * For example,
 * <pre>
TGFBG4_frame2

TGFB1           136      NTSELREAVPEPVLLSRAELRLLRLKLKVEQHVELYQKYSNNSWR 180

TGFBG4_frame2

TGFB1           181      YLSNRLLAPSDSPEWLSFDVTGVVRQWLSRG.............. 211
                                                       |
TGFBG4_frame2   1                                     LGELHSQTGFPLATPT 16

TGFB1           212      GEIEGFRLSAHCSCDSRDNTLQVDIN................... 237
                         ||||||||||||||||||||||||||
TGFBG4_frame2   17       GEIEGFRLSAHCSCDSRDNTLQVDINGEACFPGHAQLRVCVCVFP 61

TGFB1           238      ........................GFTTGRRGDLATIHGMNRPFL 258
                                                 |||||||||||||||||||||
TGFBG4_frame2   62       SAPRPTYLSLECVCMSPIPLPHKAGFTTGRRGDLATIHGMNRPFL 106

TGFB1           259      LLMATPLERAQHLQSSRHRRALDTNYCFSSTEKNCCVRQLYIDFR 303
                         |||||||||||||||||||||||||||| :
TGFBG4_frame2   107      LLMATPLERAQHLQSSRHRRALDTNYCFRALP............. 138

TGFB1           304      KDLGWKWIHEPKGYHANFCLGPCPYIWSLDTQYSKVLALYNQHNP 348
                            ||:    |    |
TGFBG4_frame2   139      ...GWR....PSRLGAL                             148

TGFB1           349      GASAAPCCVPQALEPLPIVYYVGRKPKVEQLSNMIVRSCKCS    390

TGFBG4_frame2
 * </pre>
 * <p>
 * Please note, this parser reads the whole alignment in to
 * core memory and thus does not scale to work with very large
 * alignments on low-end hardware.
 * <p>
 * Please also note that this class has not been tested with
 * many version of Needle.
 *
 * WARNING: The parser currently assumes a line length of 45
 * This restriction is easy to remove, just haven't done it.
 *
 * Copyright &copy; 2001,2002 Cambridge Antibody Technology.
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
class NeedleAlignmentSAXParser extends AbstractNativeAppSAXParser {



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
    public NeedleAlignmentSAXParser() {
	iState = STARTUP;
	this.setNamespacePrefix("biojava");
    }

    /**
     * Describe 'parse' method here.
     *
     * @param nil	 -
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

	//at this point, alignment is parsed - nibble back extraneous end
	//gaps

	this.nibbleEndGaps();


	// now cycle through and emit elements
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
     * Parse a relevant line i.e. alignment line, and add to alignment
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

	//Take from particular positions - assumes default output
	//line length etc.

	//Postions 25 to 69


	oTmpSeq = poLine.substring(25,70).replace('.','-').replace(' ','-');

	oSeq.append(oTmpSeq);


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
	//lines that start with "Needle(" not relevant (title)

	if ( (poLine.trim().equals("")) ||
	     (poLine.startsWith(" ")) ||
	     (poLine.startsWith("Global: ")) ||
	     (poLine.startsWith("Score: ")) ||
	     (poLine.startsWith("%id = ")) ||
	     (poLine.startsWith("Overall %id = ")) ) {

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
    /**
     * After the alignment is parsed, the ends of the sequences are
     * padded out with blanks to the end of the last line:
     *
     * e.g.
     * GCACACA---------------------------------
     * CAC-------------------------------------
     *
     * This method removes extraneous gaps to produce
     *
     * GCACACA
     * CAC----
     *
     */
    private void nibbleEndGaps() {


	//first get length of alignment
        String oSeqName1 = (String) oSeqNameList.get(0);
        String oSeqName2 = (String) oSeqNameList.get(1);
	String oSeq1 = (String) oAlignment.get(oSeqName1);
	String oSeq2 = (String) oAlignment.get(oSeqName2);

	int iLength = oSeq1.length();
	int iNibbles = 0;
	boolean tNibble = true;

	//index position of last character in String
	int iPos = iLength-1;
	//System.out.println("Seq1 length = " + oSeq1.length());
	//System.out.println("Seq2 length = " + oSeq2.length());

	while (tNibble) {
	    //System.out.println("-->iPos = " + iPos);
	    if ( (oSeq1.charAt(iPos) == '-') &&
                 (oSeq2.charAt(iPos) == '-') ) {
		iNibbles++;
		iPos--;
	    } else {
		tNibble = false;
	    }

	}

	//System.out.println("----------->NIBBLE " + iNibbles);
	if (iNibbles > 0) {

	    oAlignment.put(oSeqName1,oSeq1.substring(0,iLength-iNibbles));
	    oAlignment.put(oSeqName2,oSeq2.substring(0,iLength-iNibbles));

	}

    }

}
