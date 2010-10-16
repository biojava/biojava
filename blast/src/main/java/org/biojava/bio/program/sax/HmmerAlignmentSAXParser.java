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
import java.util.StringTokenizer;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;

/**
 * A reusable class for parsing Hmmr Aligment Section
 *
 * Primary author -
 *                 Colin Hardman      (CAT)
 * Other authors  -
 *                 Tim Dilks          (CAT)
 *                 Simon Brocklehurst (CAT)
 *                 Stuart Johnston    (CAT)
 *                 Lawerence Bower    (CAT)
 *                 Derek Crockford    (CAT)
 *                 Neil Benn          (CAT)
 *
 * Copyright 2001 Cambridge Antibody Technology Group plc.
 * 
 *
 * This code released to the biojava project, May 2001
 * under the LGPL license.
 *
 * @author Cambridge Antibody Technology Group plc
 * @author Greg Cox
 * @version 0.1
 *
 */
class HmmerAlignmentSAXParser extends AbstractNativeAppSAXParser {


    private BufferedReader       oContents;
    private AttributesImpl       oAtts              = new AttributesImpl();
    private QName                oAttQName          = new QName(this);
    private String               oLine;

    /**
     * Creates a new domain section parser.
     *
     * @param poVersion <code>BlastLikeVersionSupport</code>
     * @param poNamespacePrefix - the namespace prefix
     */
    HmmerAlignmentSAXParser( BlastLikeVersionSupport poVersion,
			    String poNamespacePrefix ) {
	this.setNamespacePrefix(poNamespacePrefix);
	//For XSLT Parser Compliance
	this.addPrefixMapping("biojava","http://www.biojava.org");

    }

    /**
     * Parse the buffer from the current position, until aligments
     * are parsed.
     *
     * @param poContents <code>BufferedReader</code> to parse.
     * @param poLine - the value of the current line.
     * @exception SAXException if an error occurs
     */
    public String parse( BufferedReader poContents, String poLine )
	throws SAXException {

	oContents = poContents;

	try {
	    oLine = poLine;

	    /*
	     * skip till linestarts with
	     * 'Alignments of top-scoring domains'
	     * while not line starts with 'Histogram' or '//'
	     *    case state = 0 if line not start with " "
	     *                grab to ":"
	     *                if key then set state 1, grab domain number
	     *                     else error
	     *
	     *               = 1 search for "*->" = index1
	     *                   search for "<" = index2, state 3
	     *                      else index2 = last nonwhitespace
	     *                           state =2
	     *                   grab index1+3 to index2 -> hmmmatch
	     *                   skip next line
	     *                   grab index1+3 to index2 -> seqmatch
	     *
	     *               = 2 search for "<-*" = index2, state 3
	     *                      else index2 = last nonwhitespace
	     *                   grab index1+3 to index2 -> hmmmatch
	     *                   skip next line
	     *                   grab index1 to index2 -> seqmatch
	     *
	     *               = 3 find the indexes of all '-' in seqmatch
	     *                   and store in domain's gaparray
	     *                   likewise for gaps ('.') in hmmmatch
	     *                   skip next and state = 0
	     *
	     */
	    oLine = oContents.readLine(); // skip 'Alignments...'
	    int state = 0;
	    StringBuffer oSequenceMatch = new StringBuffer();
	    StringBuffer oHmmMatch = new StringBuffer();
	    StringBuffer oMarkup = new StringBuffer();
	    int iAlignLen = 0;
	    String oScore = "";
	    String oEvalue = "";
	    String oRawSummary = "";
	    String oIdString = "";
	    String oFrom = "";
	    String oTo = "";
	    int index1 = 0;
	    int index2 = 0;

	    if ( oLine.trim().equals( "[no hits above thresholds]" ) ) {
		oLine = oContents.readLine();
		return oLine;
	    }

	    while ( (!oLine.trim().equals("//")) // hmmpfam
		    &&
		    (!oLine.trim().startsWith("Histogram of all scores:"))
		    // hmmsearch
		    ) {

		switch (state) {
		case 0:
		    if(!(oLine.trim().equals(""))){
			if(!(oLine.startsWith(" "))){

			    oMarkup.setLength( 0 );
			    oHmmMatch.setLength( 0 );
			    oSequenceMatch.setLength( 0 );

			    oRawSummary = oLine;
			    oIdString =
				oLine.substring(0, oLine.indexOf(":"));

			    StringTokenizer st =
				new StringTokenizer
				    (oLine.substring( oLine.indexOf(":") +1 ),
				     ",:" );
			    st.nextToken(); // metadata

			    String lenString = st.nextToken(); // from x to y
			    String scoreString = st.nextToken();
			    String eString = st.nextToken();

			    st = new StringTokenizer( lenString );
			    st.nextToken();  // from

			    oFrom = st.nextToken();
			    int iFrom = Integer.parseInt( oFrom );
			    st.nextToken();  // to
			    oTo = st.nextToken();
			    int iTo   = Integer.parseInt( oTo );

			    // score
			    st = new StringTokenizer( scoreString );
			    st.nextToken(); // score
			    oScore = st.nextToken();

			    st = new StringTokenizer( eString, "=" );
			    st.nextToken(); // e
			    oEvalue = st.nextToken();

			    oAtts.clear();
			    iAlignLen = (iTo-iFrom+1);

// 			    oAttQName.setQName("sequenceLength");
// 			    oAtts.addAttribute(oAttQName.getURI(),
// 					       oAttQName.getLocalName(),
// 					       oAttQName.getQName(),
// 					       "CDATA", iAlignLen + "");

			    this.startElement
				(new QName(this,this.prefix("Hit")),
				 (Attributes)oAtts);

			    oAtts.clear();
			    oAttQName.setQName("id");
			    oAtts.addAttribute(oAttQName.getURI(),
					       oAttQName.getLocalName(),
					       oAttQName.getQName(),
					       "CDATA", oIdString );

			    oAttQName.setQName("metaData");
			    oAtts.addAttribute(oAttQName.getURI(),
					       oAttQName.getLocalName(),
					       oAttQName.getQName(),
					       "CDATA", "none" );

			    this.startElement
				(new QName(this,this.prefix("HitId")),
				 (Attributes)oAtts);
			    this.endElement(new QName(this,this.prefix
						      ("HitId")));
			    // no hit description

			    oAtts.clear();
			    this.startElement
				(new QName(this,this.prefix("HSPCollection")),
				 (Attributes)oAtts);
			    this.startElement
				(new QName(this,this.prefix("HSP")),
				 (Attributes)oAtts);

			    // have to parse the aligment before
			    // can get number of positivies and identities.
			    // all goes into raw output.
			    state = 1;
			}
		    }
		    break;
		case 1:
		    index1 = oLine.indexOf("*->");
		    //	  index2 = oLine.indexOf("<-*");
		    index2 = oLine.indexOf("<");
		    if (index2 == -1) {
			index2 = oLine.trim().length() + index1;
			state =2;
		    }else{
			state =3;
		    }
		    oHmmMatch.append( oLine.substring(index1+3, index2) );
		    oLine = oContents.readLine();
		    oMarkup.append( oLine.substring(index1+3, index2) );
		    oLine = oContents.readLine();
		    oSequenceMatch.append( oLine.substring(index1+3, index2) );
		    break;
		case 2:
		    oLine = oContents.readLine(); 	  // skip blank
		    //	  index2 = oLine.indexOf("<-*");
		    index2 = oLine.indexOf("<");
		    if (index2 == -1) {
			index2 = oLine.trim().length() +index1;
		    }else{
			state =3;
		    }
		    oHmmMatch.append( oLine.substring( index1, index2 ) );
		    oLine = oContents.readLine();
		    oMarkup.append( oLine.substring(index1, index2) );
		    oLine = oContents.readLine();
		    oSequenceMatch.append( oLine.substring(index1, index2) );
		    //	  System.out.println(oSequenceMatch);
		    break;
		case 3:
		    //*****************************************************/

		    String oMarkupString = oMarkup.substring(0);
		    int iNumberOfPlus = this.countChar( oMarkupString, '+' );
		    int iNumberOfSpaces = this.countChar( oMarkupString, ' ' );

		    String oSequenceString = oSequenceMatch.substring(0);
		    String oHmmString      = oHmmMatch.substring(0);

		    int iNumberOfGaps = 0;
    //		    iNumberOfGaps +=  this.countChar( oHmmString, '.' );
		    iNumberOfGaps +=  this.countChar( oSequenceString, '-' );

		    int iAlignSize = (iAlignLen+iNumberOfGaps);

		    int iNumberOfPositives  = iAlignSize -  iNumberOfSpaces;
		    int iNumberOfIdentities = iNumberOfPositives - iNumberOfPlus;

// 		    System.err.println( "iAlignLen =\t" + iAlignLen );
// 		    System.err.println( "iNoGaps   =\t" + iNumberOfGaps );
// 		    System.err.println( "iAlignSize=\t" + iAlignSize );
// 		    System.err.println( "iNoSpaces =\t" + iNumberOfSpaces );
// 		    System.err.println( "iNoOfPlus =\t" + iNumberOfPlus );
// 		    System.err.println( "iNoOfPosi =\t" + iNumberOfPositives );

		    oAtts.clear();

		    oAttQName.setQName("score");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "CDATA", oScore );

		    oAttQName.setQName("expectValue");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "CDATA", oEvalue );

		    oAttQName.setQName("numberOfIdentities");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "CDATA", Integer.toString
				       ( iNumberOfIdentities ) );

		    oAttQName.setQName("alignmentSize");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "CDATA",  Integer.toString
				       ( iAlignSize) );

		    oAttQName.setQName("percentageIdentity");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "CDATA",
				       Integer.toString
				       (
					((int)(((double)iNumberOfIdentities/
						(double)iAlignSize)*100))
					)
				       );

		    oAttQName.setQName("numberOfPositives");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "CDATA", Integer.toString
				       ( iNumberOfPositives ) );

		    oAttQName.setQName("percentagePositives");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "CDATA",
				       Integer.toString
				       (
					((int)(((double)iNumberOfPositives/
						(double)iAlignSize)*100))
					)
				       );

		    oAttQName.setQName("numberOfGaps");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "CDATA", Integer.toString
				       ( iNumberOfGaps ) );

		    this.startElement(new QName(this,this.prefix("HSPSummary")),
				      (Attributes)oAtts);

		    //Raw HSPSummary Data
		    oAtts.clear();
		    oAttQName.setQName("xml:space");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "NMTOKEN","preserve");
		    this.startElement(new QName(this,this.prefix("RawOutput")),
				      (Attributes)oAtts);

		    char[] aoChars = oRawSummary.toCharArray();
		    this.characters( aoChars, 0, aoChars.length );

		    this.endElement(new QName(this,this.prefix("RawOutput")));
		    this.endElement(new QName(this,this.prefix("HSPSummary")));

		    // ALIGNMENT

		    oAtts.clear();
		    this.startElement(new QName(this,this.prefix
						("BlastLikeAlignment")),
				      (Attributes)oAtts);

		    oAtts.clear();
		    oAttQName.setQName("startPosition");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "CDATA", oFrom );

		    oAttQName.setQName("stopPosition");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "CDATA", oTo );

		    this.startElement(new QName(this,this.prefix
						("QuerySequence")),
				      (Attributes)oAtts);
		    aoChars = oSequenceString.toCharArray();
		    this.characters(aoChars,0,aoChars.length);
		    this.endElement(new QName(this,this.prefix
					      ("QuerySequence")));

		    //Match consensus
		    oAtts.clear();
		    oAttQName.setQName("xml:space");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "NMTOKEN","preserve");
		    this.startElement(new QName(this,this.prefix
						("MatchConsensus")),
				      (Attributes)oAtts);
		    aoChars = oMarkupString.toCharArray();
		    this.characters(aoChars,0,aoChars.length);

		    this.endElement(new QName(this,this.prefix
					      ("MatchConsensus")));

		    //HitSequence

		    oAtts.clear();
		    oAttQName.setQName("startPosition");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "CDATA", "100" );

		    oAttQName.setQName("stopPosition");
		    oAtts.addAttribute(oAttQName.getURI(),
				       oAttQName.getLocalName(),
				       oAttQName.getQName(),
				       "CDATA", "200" );

		    this.startElement(new QName(this,this.prefix("HitSequence")),
				      (Attributes)oAtts);
		    aoChars = oHmmString.toCharArray();
		    this.characters(aoChars,0,aoChars.length);
		    this.endElement(new QName(this,this.prefix("HitSequence")));

		    //end Alignment
		    this.endElement(new QName(this, this.prefix
					      (this.prefix
					       ("BlastLikeAlignment"))));
		    this.endElement(new QName(this,this.prefix("HSP")));
		    this.endElement(new QName(this,this.prefix
					      ("HSPCollection")));
		    this.endElement(new QName(this,this.prefix("Hit")));


		    state =0;
		    ; break;
		default: System.out.println("Can't reach here");
		    break;
		}
		oLine = oContents.readLine();
	    } // end while

	} catch (java.io.IOException x) {
	    System.out.println(x.getMessage());
	    System.out.println("File read interupted");
	} // end try/catch

	return oLine;
    }

    int countChar( String poString, char pcChar ) {

	int index = -1;
        int iCount = 0;
        while ( ( index = poString.indexOf( pcChar, index+1 ) ) != -1 ) {
            iCount++;
        }
        return iCount;
    }

}
