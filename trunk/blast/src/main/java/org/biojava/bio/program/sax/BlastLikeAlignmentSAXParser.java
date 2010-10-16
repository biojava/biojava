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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.StringTokenizer;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;

/**
 * A reusable class for SAX parsing Blast-alignments...
 *
 * Primary author -
 *                 Simon Brocklehurst (CAT)
 * Other authors  -
 *                 Tim Dilks          (CAT)
 *                 Colin Hardman      (CAT)
 *                 Stuart Johnston    (CAT)
 *                 Mathieu Wiepert    (Mayo Foundation)
 *
 * Copyright 2000 Cambridge Antibody Technology Group plc.
 * 
 *
 * This code released to the biojava project, May 2000
 *
 * @author Cambridge Antibody Technology Group plc
 * @author Greg Cox
 * @version 0.1
 *
 */
final class BlastLikeAlignmentSAXParser extends AbstractNativeAppSAXParser {

    private AttributesImpl       oAtts            = new AttributesImpl();
    private QName                oAttQName        = new QName(this);
    private char[]               aoChars;
    private ArrayList            oAlignment;
    private String               oLine;
    private String               oSeq;
    private StringBuffer         oQuery           = new StringBuffer();
    private StringBuffer         oHit             = new StringBuffer();
    private StringBuffer         oMatchConsensus  = new StringBuffer();
    private StringBuffer         oStartId         = new StringBuffer();
    private StringBuffer         oStopId          = new StringBuffer();
    private StringBuffer         oHitStartId      = new StringBuffer();
    private StringBuffer         oHitStopId       = new StringBuffer();
    private StringBuffer         oQueryStartId    = new StringBuffer();
    private StringBuffer         oQueryStopId     = new StringBuffer();
    private StringTokenizer      oSt;
    private String               oParsedSeq;
    private int                  iOffset;
    private int                  iEnd;
    private boolean              tJustDoneConsensus;

    private static final int STARTUP             = 0;
    private static final int ON_FIRST_SEGMENT    = 1;
    private static final int DONE_FIRST_SEGMENT  = 2;



    public BlastLikeAlignmentSAXParser(String poNamespacePrefix) {
    this.changeState(STARTUP);
    this.setNamespacePrefix(poNamespacePrefix);
    this.addPrefixMapping("biojava", "http://www.biojava.org");

    }

    public void parse(ArrayList poAlignment)
    throws SAXException {

    oAlignment = poAlignment;

    oAtts.clear();
    this.startElement(new QName(this,this.prefix("BlastLikeAlignment")),
              (Attributes)oAtts);

    this.changeState(ON_FIRST_SEGMENT);

    //for a new alignment initialise

    oQuery.setLength(0);
    oQueryStartId.setLength(0);
    oQueryStopId.setLength(0);
    oHit.setLength(0);
    oHitStartId.setLength(0);
    oHitStopId.setLength(0);
    oMatchConsensus.setLength(0);

    tJustDoneConsensus = false;
    //Loop over all alignment lines
    int iAlSize = oAlignment.size();
    for (int i = 0; i < iAlSize;i++) {
        //System.out.println(oAlignment.get(i));
        oLine = (String)oAlignment.get(i);
        this.parseLine(oLine);

    }

    //at this point alignment is parsed

//  System.out.println("QueryStart:"+oQueryStartId);
//  System.out.println("QueryStop:"+oQueryStopId);
//  System.out.println("HitStart:"+oHitStartId);
//  System.out.println("HitStop:"+oHitStopId);
//  System.out.println("Query:"+oQuery);
//  System.out.println("Match:"+oMatchConsensus);
//  System.out.println("Hit  :"+oHit);


    //output elements
    //QuerySequence

    oAtts.clear();
    oAttQName.setQName("startPosition");
    oAtts.addAttribute(oAttQName.getURI(),
               oAttQName.getLocalName(),
               oAttQName.getQName(),
               "CDATA",oQueryStartId.substring(0));

    oAttQName.setQName("stopPosition");
    oAtts.addAttribute(oAttQName.getURI(),
               oAttQName.getLocalName(),
               oAttQName.getQName(),
               "CDATA",oQueryStopId.substring(0));

    this.startElement(new QName(this,this.prefix("QuerySequence")),
              (Attributes)oAtts);
    aoChars = oQuery.substring(0).toCharArray();
    this.characters(aoChars,0,aoChars.length);
    this.endElement(new QName(this,this.prefix("QuerySequence")));


    //Match consensus
    oAtts.clear();
    oAttQName.setQName("xml:space");
    oAtts.addAttribute(oAttQName.getURI(),
               oAttQName.getLocalName(),
               oAttQName.getQName(),
               "NMTOKEN","preserve");
    this.startElement(new QName(this,this.prefix("MatchConsensus")),
              (Attributes)oAtts);
    aoChars = oMatchConsensus.substring(0).toCharArray();
    this.characters(aoChars,0,aoChars.length);

    this.endElement(new QName(this,this.prefix("MatchConsensus")));

    //HitSequence

    oAtts.clear();
    oAttQName.setQName("startPosition");
    oAtts.addAttribute(oAttQName.getURI(),
               oAttQName.getLocalName(),
               oAttQName.getQName(),
               "CDATA",oHitStartId.substring(0));

    oAttQName.setQName("stopPosition");
    oAtts.addAttribute(oAttQName.getURI(),
               oAttQName.getLocalName(),
               oAttQName.getQName(),
               "CDATA",oHitStopId.substring(0));

    this.startElement(new QName(this,this.prefix("HitSequence")),
              (Attributes)oAtts);
    aoChars = oHit.substring(0).toCharArray();
    this.characters(aoChars,0,aoChars.length);
    this.endElement(new QName(this,this.prefix("HitSequence")));


    //end Alignment
    this.endElement(new QName(this,
                  this.prefix(this.prefix("BlastLikeAlignment"))));
    }
    /**
     * Describe 'parseLine' method here.
     *
     * @param poLine     -
     */
    private void parseLine(String poLine) throws SAXException{

    poLine = poLine.toUpperCase();

    if ( (poLine.startsWith("QUERY:")) ||
             (poLine.startsWith("SBJCT:")) ) {
        oSt = new StringTokenizer(poLine,":");

        //there should be two tokens at this point
        if (oSt.countTokens() != 2) {
        throw (new SAXException(
        "Failed to parse a line in a BlastLikeAlignment" +
        " due it having an unexpected format." +
        " The line is shown below.\n" +
        poLine));
        }

        //get here if Query line OK

        //skip first token (i.e. "Query")
        oSt.nextToken();

        //next token is the alignment - make it uppercase.
        oSeq = oSt.nextToken().trim();

       //System.out.println(oSeq);
       //To get numbers robustly, tokenize on letters, gaps, and unknowns

       oSt = new StringTokenizer(oSeq," ABCDEFGHIJKLMNOPQRSTUVWXYZ-*");

       //System.out.println("Token Count----->" + oSt.countTokens());

       //throw exception if there number of tokens is not two
       //(these correspond to start and stop ids
       if (oSt.countTokens() != 2) {
           throw (new SAXException(
           "Failed to parse a line of an alignment due to it having" +
        " an unexpected character."));
       }
       //here if tokens for start and stop OK

       oStartId.setLength(0);
       oStartId.append(oSt.nextToken().trim());

       oStopId.setLength(0);
       oStopId.append(oSt.nextToken().trim());

           //System.out.println("StartId="+oStartId+" : "+"StopId="+oStopId);


       //To get sequence robustly, tokenize on numbers only

       oSt = new StringTokenizer(oSeq,"0123456789");

       //System.out.println("Token Count----->" + oSt.countTokens());

       if (oSt.countTokens() != 1) {
           throw (new SAXException(
           "Failed to parse a line of an alignment due to it having" +
        " an unexpected character."));
       }

       oParsedSeq = oSt.nextToken().trim();

       //System.out.println(oParsedSeq);

       //get info for consensus only on Query: lines

       iOffset = poLine.indexOf(oParsedSeq);
       iEnd = iOffset + oParsedSeq.length();


       //System.out.println("Offset="+iOffset+" : End="+iEnd);
       //end if this is a query or sbjct line
    } else {
        //here if on a consensus sequence line
        //only time should get here is if
        //a Query: line has just been parsed.

        //deal with software that doesn't output spaces to end of
        //consensus line

        if (iEnd <= poLine.length()) {
	    oParsedSeq = poLine.substring(iOffset,iEnd);
        } else {
	    int iLen = iEnd - poLine.length();
	    char[] oPadding = new char[ iLen ];
	    Arrays.fill( oPadding,
			 0,
			 iLen,
			 ' ' );
	    oParsedSeq = poLine.substring( iOffset).concat
		( new String( oPadding ) );
	}
    }

    //get startIds for query and subject
    if (iState == ON_FIRST_SEGMENT) {
        //here if on first block of an alignment
        if (poLine.startsWith("QUERY:")) {
	    oQueryStartId.append(oStartId);
	    oQueryStopId.append(oStopId);
	    oQuery.append(oParsedSeq);
	    tJustDoneConsensus = false;
	    return;
        }
        if (poLine.startsWith("SBJCT:")) {
	    oHitStartId.append(oStartId);
	    oHitStopId.append(oStopId);
	    oHit.append(oParsedSeq);

	    if (!tJustDoneConsensus) {

		//handle rare case of a totally blank
		//consensus line
		char[] oPadding = new char[ iEnd-iOffset ];
		Arrays.fill( oPadding,
			     0,
			     iEnd-iOffset,
			     ' ' );
		oMatchConsensus.append( new String( oPadding) );
	    }

	    tJustDoneConsensus = false;

	    //here finished with block
	    this.changeState(DONE_FIRST_SEGMENT);
	    return;
        }
        oMatchConsensus.append(oParsedSeq);

        tJustDoneConsensus = true;

    } //end if onFirstSegment

    //if inside the alignment, set the stopids each time
    //so that they are correct for multi-block alignments
    if (iState == DONE_FIRST_SEGMENT) {

        if (poLine.startsWith("QUERY:")) {
	    oQueryStopId.setLength(0);
	    oQueryStopId.append(oStopId);
	    oQuery.append(oParsedSeq);
	    tJustDoneConsensus = false;
	    return;
        }
        if (poLine.startsWith("SBJCT:")) {
	    oHitStopId.setLength(0);
	    oHitStopId.append(oStopId);
	    oHit.append(oParsedSeq);
	    //here finished with block
	    if (!tJustDoneConsensus) {

		//handle rare case of a totally blank
		//consensus line
		char[] oPadding = new char[ iEnd-iOffset ];
		Arrays.fill( oPadding,
			     0,
			     iEnd-iOffset,
			     ' ' );
		oMatchConsensus.append( new String( oPadding) );

	    }
	    tJustDoneConsensus = false;
	    return;
        }
        //get here if on a match consensus
        oMatchConsensus.append(oParsedSeq);
        tJustDoneConsensus = true;
        return;
    } //end if inAlignment
    }

}
