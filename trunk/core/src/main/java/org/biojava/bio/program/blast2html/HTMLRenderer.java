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

package org.biojava.bio.program.blast2html;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;

/**
 * Renders HTML version of blast-like output.<p>
 *
 * Makes an assumption that the gap character is '-',
 * should parameterize.
 *
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
 * This code released to the biojava project, May 2001
 * under the LGPL license.
 *
 * @author Cambridge Antibody Technology Group plc
 * @author Greg Cox
 * @version 1.0
 */
public class HTMLRenderer  {

    private StringBuffer align = new StringBuffer();
    private StringBuffer alignMarkUp1 = new StringBuffer();
    private StringBuffer alignMarkUp2 = new StringBuffer();

    private PrintWriter out = null;
    private String      oStyleDef;
    private StringBuffer pcDataBuffer = new StringBuffer( 80 );

    private boolean tAlternateSummary = false;
    private int     iSummaryCount =0;

    private boolean tNeedComma = false;

    /**
     * The width, in characters, of the sequence alignments
     */
    private int iAlignLen = 50;

    private List             oURLGeneratorList = null;
    private DatabaseURLGenerator oFirstURLGenerator = null;

    private AlignmentMarker oAlignmentMarker;

    private static String oJavaScriptDef
	= "<script language=\"JavaScript\">\nfunction setStatus(msgStr) { \n status=msgStr;\n  document.retVal = true;\n}\n</script>";

    private char[] padding = new char[100];

    {
	Arrays.fill( padding,
		     0,
		     100,
		     ' ' );
    }

    private Properties oOptions = null;

    /**
     * Flag to detect if HTML output was empty.
     */
    private boolean wasEmpty = true;

    /**
     * Creates an HTMLRenderer, that outputs the HTML to the specified
     * PrintWriter. <p>
     *
     * The style definition is expected to defines the following styles
     * <UL>
     *  <LI> footer
     *  <LI> alignment
     *  <LI> dbRetrieve
     *  <LI> titleLevel1
     *  <LI> titleLevel1Sub
     *  <LI> titleLevel2
     *  <LI> titleLevel3
     *  <LI> titleLevel3Sub
     *  <LI> titleLevel4
     *  <LI> summaryBodyLineOdd
     *  <LI> summaryBodyLineEven
     * </UL
     *
     * The Javascript function setStatus overrides status text for some links.
     *
     * The URLGeneratorFactory provides a list of objects that can convert a
     * database ID into a URL for fetching the entry. Obviously this requires
     * the ID to be parsed correctly.<B>
     * These URL's are used in two places, in the summary and detail sections.
     * The first URLGenerator is used for the summary, and all for the detail.
     *
     *
     * @param poPrintWriter - the output stream to write the HTML to.
     * @param poStyleDef - definition of styles, if null uses hard coded.
     * @param piAlignmentWidth width in characters of the alignment regions
     *                         of output
     * @param poFactory <code>URLGeneratorFactory</code> provides an array of
     *                  <code>DatabaseURLGenerator</code>s, null for no links
     * @param poAlignmentMarker - for configurable markup of alignments,
     *                            for no markup use null.
     * @param poOptions properties of
     */
    public HTMLRenderer( PrintWriter poPrintWriter,
			 String      poStyleDef,
			 int         piAlignmentWidth,
			 URLGeneratorFactory poFactory,
			 AlignmentMarker poAlignmentMarker,
			 Properties  poOptions ) {

	oOptions = poOptions;

	if ( poFactory != null ) {

	    oURLGeneratorList = poFactory.getDatabaseURLGenerators();
	    if ( !oURLGeneratorList.isEmpty() ) {
		oFirstURLGenerator
		    = (DatabaseURLGenerator) oURLGeneratorList.get( 0 );
	    }
	}

	if ( piAlignmentWidth <= 0 ) {
	    throw new IllegalArgumentException
		( "Alignment length must be > 0, not " + piAlignmentWidth );
	}

	iAlignLen = piAlignmentWidth;
	oStyleDef = poStyleDef ;

	if ( poPrintWriter == null ) {
	    throw new IllegalArgumentException
		( "PrintWriter cannot be null" );
	}

	out = poPrintWriter;
	oAlignmentMarker = poAlignmentMarker;
    }

    /**
     * Set the <CODE>PrintWriter</CODE> to output the HTML
     * to.
     *
     * @param poPrintWriter a <code>PrintWriter</code>
     */
    public void setPrintWriter( PrintWriter poPrintWriter ) {

	out = poPrintWriter;

	// initialize other stuff
	tAlternateSummary = false;
	iSummaryCount =0;
	tNeedComma = false;
	wasEmpty = true;

	align.setLength( 0 );
	alignMarkUp1.setLength( 0 );
	alignMarkUp2.setLength( 0 );
	pcDataBuffer.setLength( 0 );
    }

    /**
     * Returns true if  no HitSummary and no Hit elements were
     * encountered
     *
     * @return <code>boolean</code> - true if no summary of Hit elements
     *          where encountered.
     */
    public boolean wasEmpty() {
	return wasEmpty;
    }

    /**
     * Writes out the header section of the report.
     *
     * @param oProgram a <code>String</code> - eg 'ncbi-blastp'
     * @param oVersion a <code>String</code> - eg '2.2.1'
     * @param oQuery a <code>String</code>   - eg 'my query id'
     * @param oDatabase a <code>String</code>- eg 'embl'
     */
    void writeTitleAndHeader( String oProgram, String oVersion,
			      String oQuery,   String oDatabase ) {

	out.print( "<TITLE>" );
	out.print( oProgram );
	out.print( " version " );
	out.print( oVersion );
	out.print( " Search Results Query ID : " );
	out.print( oQuery );
	out.println( "</TITLE>" );

	out.println( "<table width=\"100%\" border=\"0\" cellpadding=\"3\" cellspacing=\"0\">" );
	out.println( "  <tr align=\"left\" valign=\"top\">" );
	out.println( "    <td class=\"titleLevel1\" rowspan=\"3\"><a name=\"#AnchorTop\"></a>Sequence Similarity Report</td>" );
	out.println( "    <td align=\"right\" class=\"titleLevel1Sub\"> Query Id:</td>" );
	out.print( "    <td class=\"titleLevel1Sub\">" );
	out.print( oQuery  );
	out.println( "&nbsp;</td>" );
	out.println( "  </tr><tr align=\"left\" valign=\"top\" class=\"level1titleContainer\"> " );
	out.println( "    <td class=\"titleLevel1Sub\" align=\"right\">Search&nbsp;Program:</td>" );
	out.print( "    <td class=\"titleLevel1Sub\">" );
        out.print( oProgram );
	out.print( " v.&nbsp;" );
	out.print( oVersion );
	out.println( "</td>" );
	out.println( "  </tr><tr align=\"left\" valign=\"top\" class=\"level1titleContainer\">" );
	out.println( "    <td class=\"titleLevel1Sub\" align=\"right\">Database:</td>" );
	out.print( "    <td class=\"titleLevel1Sub\">" );
	out.print( oDatabase );
	out.println( "</td>" );
	out.println( "  </tr>" );
	out.println( "</table>" );
    }


    /**
     * Returns the appropriate style and javascript definitions for this
     * renderer.
     *
     * @return a <code>String</code> value
     */
    public String getHeaderDefinitions() {

	StringBuffer sb = new StringBuffer();

	sb.append("<STYLE TYPE=\"text/css\">\n");
	sb.append("<!--\n");
	sb.append( oStyleDef );
 	if ( oAlignmentMarker != null ) {
 	    sb.append( oAlignmentMarker.getAlignmentStyles() );
 	}
	sb.append( "-->\n</STYLE>\n" );
	sb.append( oJavaScriptDef );
	return sb.substring(0);
    }


    /**
     * Called when first summary item is reached.
     *
     * @param oHitSummary a <code>HitSummary</code> - the first Summary item.
     */
    void startSummaryTable( HitSummary oHitSummary ) {

	wasEmpty = false;

	out.println( "<br>" );
	out.println( "<table width=\"100%\" border=\"0\" cellspacing=\"0\" cellpadding=\"0\">" );
	out.println( "<tr valign=\"middle\"> " );
     	out.println( "<td colspan=\"4\" class=\"titleLevel2\" align=\"left\">Overview of Results</td>" );
     	out.println( "</tr>" );
	out.println( "</table>" );

	out.println( "<table width=\"100%\" border=\"0\" cellspacing=\"0\" cellpadding=\"0\">" );

 	out.println( "<td class=\"titleLevel3\" >" );
        out.println( "<table  width=\"100%\" border=\"0\" cellspacing=\"0\" cellpadding=\"0\">" );
 	out.println( "<tr> " );

	out.println( "<td class=\"titleLevel3\" width=\"40%\" align=\"left\" >Hit Id</td>" );
 	out.println( "<td class=\"titleLevel3\" width=\"60%\" align=\"left\">Hit Description</td>" );
 	out.println( "</tr>" );

 	out.println( "</table>" );
 	out.println( "</td> " );

	if ( oHitSummary.readingFrame != null ) {
	    out.println( "<td width=\"7%\" class=\"titleLevel3\" align=\"center\">Frame</td>" );
	}

 	out.println( "<td width=\"7%\" class=\"titleLevel3\" align=\"center\">Score</td>" );

	if ( oHitSummary.expectValue != null ) {
	    out.println( "<td width=\"7%\" class=\"titleLevel3\" align=\"center\">E value</td>" );
	}
	if ( oHitSummary.smallestSumProbability != null ) {
	    out.println( "<td width=\"7%\" class=\"titleLevel3\" align=\"center\">P(N)</td>" );
	}
	if ( oHitSummary.numberOfHSPs != null ) {
	    out.println( "<td width=\"7%\" class=\"titleLevel3\" align=\"left\">HSPs</td>" );
	}
	if ( oHitSummary.numberOfContributingHSPs != null ) {
	    out.println( "<td width=\"7%\" class=\"titleLevel3\" align=\"center\">HSPs</td>" );
	}

 	out.println( "</tr>" );
	// start a new summary table
	iSummaryCount = 10;
    }


    /**
     * Called when a hit summary is reached ( except first ).
     *
     * @param oHitSummary a <code>HitSummary</code> value
     */
    void writeCurrentSummary( HitSummary oHitSummary ) {

	if ( iSummaryCount == 10 ) { // lets break up that big table
	    iSummaryCount = 0;

	    out.println( "</table>" );
	    out.println
	("<table width=\"100%\" border=\"0\" cellspacing=\"0\" cellpadding=\"0\">");
	}

	out.println( "<tr valign=\"middle\" align=\"left\" > " );

	out.print( "<td class=\"" );
	if ( tAlternateSummary ) {
	    out.print( "summaryBodyLineEven" );
	} else {
	    out.print( "summaryBodyLineOdd" );
	}
	out.print( "\">" );

 	out.println
 	    ("<table width=\"100%\" border=\"0\" cellspacing=\"0\" cellpadding=\"0\">");
 	out.println( "<tr> " );

	out.print( "<td width=\"40%\" class=\"" );

	if ( tAlternateSummary ) {
	    out.print( "summaryBodyLineEven" );
	} else {
	    out.print( "summaryBodyLineOdd" );
	}
	out.print( "\" >" );

	if ( oFirstURLGenerator!= null ) {

	    out.print( "<a href=\"" );
	    out.print( oFirstURLGenerator.toURL( oHitSummary.oHitId.id,
						 oOptions ));
	    out.print( "\" onMouseOver=\"setStatus('Retrieve hit from database');return document.retVal\" onMouseOut=\"setStatus(' ');return document.retVal\" class=\"dbRetrieve\">" );
	    out.print( oHitSummary.oHitId.id );
	    out.print( "</a></td>" );
	} else {
	    out.print( oHitSummary.oHitId.id );
	    out.println( "</td>" );
	}

	out.print( "<td width=\"60%\" class=\"" );
	if ( tAlternateSummary ) {
	    out.print( "summaryBodyLineEven" );
	} else {
	    out.print( "summaryBodyLineOdd" );
	}
	out.print( "\"><a href=\"#Anchor" );
	out.print( oHitSummary.oHitId.id );
	out.print( "\" onMouseOver=\"setStatus('Jump to detailed results');return document.retVal\" onMouseOut=\"setStatus(' ');return document.retVal\">" );

	this.writeContentChars( oHitSummary.oDesc.hitDescription );

	out.println( "</a></td>" );
	out.println( "</tr>" );

	out.println( "</table>" );
 	out.println( "</td> " );

	// reading frame
	this.writeSummaryColumn( oHitSummary.readingFrame, "center" );

	// score
 	out.print( "<td width=\"7%\" align=\"center\" class=\"" );
 	if ( tAlternateSummary ) {
 	    out.print( "summaryBodyLineEven" );
 	} else {
 	    out.print( "summaryBodyLineOdd" );
 	}
 	out.print( "\"><a href=\"#Anchor" );
 	out.print( oHitSummary.oHitId.id );
 	out.print( "\" onMouseOver=\"setStatus('Jump to detailed results');return document.retVal\" onMouseOut=\"setStatus(' ');return document.retVal\">" );
 	out.print(  oHitSummary.score );
 	out.println( "</a></td>" );

	// Expect Value
	this.writeSummaryColumn( oHitSummary.expectValue, "center" );
	// Smallest sum prob
	this.writeSummaryColumn( oHitSummary.smallestSumProbability, "center" );
	// numberOfHSPs
	this.writeSummaryColumn( oHitSummary.numberOfHSPs, "center");
	// numberOfContributingHSPs
	this.writeSummaryColumn( oHitSummary.numberOfContributingHSPs, "center" );
	out.println( "</tr>" );

	tAlternateSummary = !tAlternateSummary;
	iSummaryCount ++;
    }

    /**
     * Outputs a summary column item.
     *
     */
    void writeSummaryColumn( String poValue, String poAlign ) {

	if ( poValue != null ) {

	    out.print( "<td width=\"7%\" align=\"" );
	    out.print( poAlign );
	    out.print( "\" class=\"" );
	    if ( tAlternateSummary ) {
		out.print( "summaryBodyLineEven" );
	    } else {
		out.print( "summaryBodyLineOdd" );
	    }
	    out.print( "\">" );
	    out.print( poValue );
	    out.println( "</td>" );
	}
    }

    /**
     * Called when summary table is complete.
     *
     */
    void endSummaryTable() {

	out.println( "</table>" );
	out.println( "<br>" );

    }


    /**
     * Start the Detail table.
     *
     */
    void startDetailTable() {

	wasEmpty = false;

	out.println( "<br>" );
	out.println( "<table width=\"100%\" border=\"0\" cellspacing=\"0\" cellpadding=\"0\">" );
	out.println( "<tr valign=\"middle\"> " );
     	out.println( "<td class=\"titleLevel2\" align=\"left\">Detailed Analysis of Results</td>" );
     	out.println( "</tr>" );
	out.println( "</table>" );
    }

    /**
     * Write the current detail.
     *
     */
    void writeCurrentDetail( DetailHit oDetailHit ) {

	out.println
	    ( "<table width=\"100%\" border=\"0\" cellspacing=\"0\" cellpadding=\"5\">" );
	out.println( "<tr align=\"left\" valign=\"top\"> " );
	out.print( "<td class=\"titleLevel3\"> <a name=\"Anchor" );
	out.print( oDetailHit.oHitId.id );
	out.print( "\"></a>Hit Id : " );
	out.print( oDetailHit.oHitId.id );
	out.println( "<br><span class=\"titleLevel3Sub\">" );
	this.writeContentChars( oDetailHit.oDesc.hitDescription );
	out.println( "</span></td>" );
	out.println( "<td class=\"titleLevel3\" align=\"right\" color=\"#000000\" ><a href=\"#AnchorTop\" onMouseOver=\"setStatus('Jump to top of  document');return document.retVal\" onMouseOut=\"setStatus(' ');return document.retVal\">Top</a></td>" );
	out.println( "</tr>" );
	out.println( "<tr align=\"left\" valign=\"top\"> " );
	out.println( "<td colspan=\"2\"><span class=\"alignment\">" );
	out.print( "<p>Sequence length of hit = " );
	out.print( oDetailHit.sequenceLength);
	out.println( "</span></td>" );
	out.println( "<br>" );

	if ( oURLGeneratorList != null && !oURLGeneratorList.isEmpty() ) {
	    out.println( "<tr align=\"left\" valign=\"top\"> " );

	    out.print( "<td colspan=\"" );
	    out.print( oURLGeneratorList.size() + "" );
	    out.println( "\"> " );

	    for ( int i =0 ; i < oURLGeneratorList.size() ; i++ ) {

		DatabaseURLGenerator oURLGenerator = ( DatabaseURLGenerator )
		    oURLGeneratorList.get( i );
		out.print( oURLGenerator.toLink( oDetailHit.oHitId.id,
						 oOptions ));
	    }
	    out.println( "</p></td>" );
	    out.println( "</tr>" );
	}

	out.print( "</table>" );
	//	out.println( "<br>" );
     }



    /**
     * Utility for writing out each HSP info item, such as
     * score ot number of identities.
     *
     */
    void writeHSPInfo( String poName, String poValue ) {

	if ( poValue != null ) {
	    if ( tNeedComma ) {
	    out.print( ", " );
	    }
	    out.print( poName );
	    out.print( "&nbsp;=&nbsp;" );
	    out.print( poValue );
	    tNeedComma = true;
	}

    }

    /**
     * Writes out the current HSP.
     *
     */
    void writeCurrentHSP( HSPSummary oHSPSummary, BlastLikeAlignment oAlignment ) {
	 tNeedComma = false;

	 out.println
	     ( "<table width=\"100%\" border=\"0\" cellspacing=\"0\" cellpadding=\"5\">" );
	 out.println( "<tr align=\"left\" valign=\"top\"> " );
	 out.println
	     ( "<td class=\"titleLevel4\">High-scoring segment pair (HSP) group</td>" );
	 out.println( "</tr>" );
	 out.println( "</table>" );
	 out.println( "<div align=\"right\"><br>" );
	 out.println
	     ( "<table width=\"95%\" border=\"0\" cellspacing=\"0\" cellpadding=\"5\">" );
	 out.println( "<tr align=\"left\" valign=\"top\" class=\"level4title\"> " );
	 out.println( "<td class=\"titleLevel4\"> " );
	 out.println( "<p>HSP Information<br>" );

	 this.writeHSPInfo( "Score", oHSPSummary.score );
	 this.writeHSPInfo( "E", oHSPSummary.expectValue );
 	 this.writeHSPInfo( "P", oHSPSummary.pValue );

	 if ( oHSPSummary.numberOfIdentities != null ||
	      oHSPSummary.percentageIdentity != null ) {

	     out.print( ", Identities&nbsp;=&nbsp;" );
	     if ( oHSPSummary.numberOfIdentities != null ) {
		 out.print( oHSPSummary.numberOfIdentities );

		 if ( oHSPSummary.alignmentSize != null ) {
		     out.print( "/" );
		     out.print( oHSPSummary.alignmentSize );
		 }

		 if ( oHSPSummary.percentageIdentity != null ) {
		     out.print( "&nbsp;(" );
		     out.print( oHSPSummary.percentageIdentity );
		     out.print( "%)" );
		 }
	     }
	 }

	 if ( oHSPSummary.numberOfPositives != null ||
	      oHSPSummary.percentagePositives != null ) {

	     out.print( ", Positives&nbsp;=&nbsp;" );
	     if ( oHSPSummary.numberOfPositives != null ) {
		 out.print( oHSPSummary.numberOfPositives );

		 if ( oHSPSummary.alignmentSize != null ) {
		     out.print( "/" );
		     out.print( oHSPSummary.alignmentSize );
		 }

		 if ( oHSPSummary.percentagePositives != null ) {
		     out.print( "&nbsp;(" );
		     out.print( oHSPSummary.percentagePositives );
		     out.print( "%)" );
		 }
	     }
	 }

	 this.writeHSPInfo( "Length", oHSPSummary.alignmentSize );
	 this.writeHSPInfo( "Query Frame",
			    this.toSign( oHSPSummary.queryFrame ) );
	 this.writeHSPInfo( "Hit Frame",
			    this.toSign( oHSPSummary.hitFrame ) );

	 if ( oHSPSummary.queryFrame == null ) {
	     this.writeHSPInfo( "Query Strand",
				this.toSign( oHSPSummary.queryStrand ) );
	 }
	 if ( oHSPSummary.hitFrame == null ) {
	     this.writeHSPInfo( "Hit Strand",
				this.toSign( oHSPSummary.hitStrand ));
	 }

 	 this.writeHSPInfo( "P(N)", oHSPSummary.sumPValues );
 	 this.writeHSPInfo( "No.&nbsp;of&nbsp;gaps", oHSPSummary.numberOfGaps );

         out.println( "</p>" );
         out.println( "</td>" );
	 out.println( "</tr>" );
	 out.println( "</table>" );
	 out.println( "</div>" );

	 out.println( "<div align=\"right\">" );
	 out.println( "<table width=\"95%\" border=\"0\" cellspacing=\"0\" cellpadding=\"5\">" );
	 out.println( "<tr> " );
	 out.println( "<td> </td>" );
	 out.println( "</tr>" );

	 this.drawCurrentAlignment( oAlignment );

	 out.println( "</table>" );
	 out.println( "</div>" );
     }

    /**
     * Draws one block of the alignment.
     */
    String drawSubAlignment( String piQueryStart,
			     String piQueryStop,
			     String piHitStart,
			     String piHitStop,
			     String poQuery,
			     String poConsensus,
			     String poHit ) {

	align.setLength( 0 );

	int iMax = 4;

 	if ( piQueryStart.length() > iMax ) {
 	    iMax = piQueryStart.length();
 	}
 	if ( piHitStart.length() > iMax ) {
 	    iMax = piHitStart.length();
 	}
	iMax = iMax+2;

	String[] oFormattedSeq = new String[]{ poQuery, poHit, poConsensus };
	if ( oAlignmentMarker != null ) {
	    oAlignmentMarker.alignment2HTML( oFormattedSeq );
	}

	String oConsensusPad = "       ";
	oConsensusPad = oConsensusPad.concat( this.padTo( "", iMax) );

	align.append( "<tr align=\"left\" valign=\"top\"> " );
	align.append( "\n");
	align.append( "<td>" );
	align.append( "\n");
	align.append( "<pre><span class=\"alignment\">Query: " );
	align.append( this.padTo( piQueryStart, iMax) );

	align.append( oFormattedSeq[0] );
	align.append( "  ");
	align.append( this.padTo( piQueryStop, iMax) );

	align.append( "\n");

	align.append( oConsensusPad );

	align.append( oFormattedSeq[2] );
	align.append( "\n");
	align.append( "Hit  : " );
	align.append( this.padTo( piHitStart, iMax ) );

	align.append( oFormattedSeq[1] );
	align.append( "  ");
	align.append( this.padTo( piHitStop, iMax) );
	align.append( "</span></PRE> " );
	align.append( "\n");
	align.append( "" );
	align.append( "\n");
	align.append( "</td>" );
	align.append( "\n");
	align.append( "</tr>" );
	align.append( "\n");

	return align.substring(0);
    }


    /**
     * Draws a full alignment block.
     */
    void drawCurrentAlignment( BlastLikeAlignment oAlignment ) {

	int i = 0;

	int iCurrentQueryStart = Integer.parseInt
	    ( oAlignment.oQuerySeq.startPosition );
	int iCurrentHitStart = Integer.parseInt
	    ( oAlignment.oHitSeq.startPosition );

	int iQueryLen = oAlignment.oQuerySeq.seq.length();
	int iHitLen = oAlignment.oHitSeq.seq.length();

	int iNumberOfQueryGaps = 0;
	int iNumberOfHitGaps = 0;
	int index = -1;
	while ( ( index = oAlignment.oQuerySeq.seq.indexOf( '-', index+1 ) ) != -1 ) {
	    iNumberOfQueryGaps++;
	}
	index = -1;
	while ( ( index = oAlignment.oHitSeq.seq.indexOf( '-', index+1 ) ) != -1 ) {
	    iNumberOfHitGaps++;
	}

	int iQStop = Integer.parseInt( oAlignment.oQuerySeq.stopPosition );
	int iHStop = Integer.parseInt( oAlignment.oHitSeq.stopPosition );


        int queryDirection = 1;
	int hitDirection   = 1;

	if ( iQStop < iCurrentQueryStart ) {
	    queryDirection = -1;
	}
	if ( iHStop < iCurrentHitStart ) {
	    hitDirection = -1;
	}

  	int iQueryMultiplier = (( iQStop  - iCurrentQueryStart ) + queryDirection)/
 	    ( iQueryLen - iNumberOfQueryGaps );

  	int iHitMultiplier   = (( iHStop - iCurrentHitStart ) + hitDirection )/
 	    ( iHitLen - iNumberOfHitGaps );

	int iCurrentQueryEnd = 0;
	int iCurrentHitEnd   = 0;


	//
	//
	// Substring  ( i*iAlignLen, (i+1)*iAlignLen )
	//
	// Increment the end number by ( (iAlign-numberofgaps)* multiplier )
	//
	// The end check should be the current end number
	//

	while( ((i+1)*iAlignLen) < iQueryLen ) {

	    String oCurrentQueryString =  oAlignment.oQuerySeq.seq.substring
		( i*iAlignLen, (i+1)*iAlignLen );
	    String oCurrentHitString   =  oAlignment.oHitSeq.seq.substring
		( i*iAlignLen, (i+1)*iAlignLen );

	    iNumberOfQueryGaps = this.countNumberOfGaps( oCurrentQueryString );
	    iNumberOfHitGaps   = this.countNumberOfGaps( oCurrentHitString );

	    iCurrentQueryEnd = iCurrentQueryStart +
		( ( iAlignLen - iNumberOfQueryGaps ) * iQueryMultiplier );
	    iCurrentHitEnd =   iCurrentHitStart   +
		(( iAlignLen - iNumberOfHitGaps   ) * iHitMultiplier );

	    out.println( drawSubAlignment
			 ( iCurrentQueryStart + "",
			   (iCurrentQueryEnd - queryDirection) + "",
			   iCurrentHitStart + "",
			   (iCurrentHitEnd - hitDirection) + "",
			   oCurrentQueryString,
			   oAlignment.oConsensus.substring
			   ( i*iAlignLen, (i+1)*iAlignLen ) ,
			   oCurrentHitString )
			 );
	    i++;
	    iCurrentQueryStart += ( ( iAlignLen - iNumberOfQueryGaps )
				    * iQueryMultiplier );
	    iCurrentHitStart   += ( ( iAlignLen - iNumberOfHitGaps   )
				    * iHitMultiplier );

	} // end while

	if ( iQStop != iCurrentQueryEnd ) {

	    iCurrentQueryEnd = iQStop;
	    iCurrentHitEnd = iHStop;

	    out.println( drawSubAlignment( iCurrentQueryStart + "",
					   iCurrentQueryEnd + "",
					   iCurrentHitStart + "",
					   iCurrentHitEnd + "",
					   oAlignment.oQuerySeq.seq.substring
					   ( i*iAlignLen ) ,
					   oAlignment.oConsensus.substring
					   ( i*iAlignLen ) ,
					   oAlignment.oHitSeq.seq.substring
					   ( i*iAlignLen ) )
			 );
	}

    }

    /**
     * Renderers end of detail table
     *
     */
    void endDetailTable() {
	out.println
	    ( "<table width=\"100%\" border=\"0\" cellpadding=\"0\" cellspacing=\"0\">" );
	out.println( "<tr> " );
	out.println( "<td class=\"footer\" height=\"21\"> " );
	out.println
    ("<div align=\"center\"> <a href=\"http://www.biojava.org/\" class=\"footer\"> </a> " );
        out.println
     ( "Produced using <a href=\"http://www.biojava.org/\" class=\"footer\">biojava</a> " );
        out.println( "- www.biojava.org</div>" );
	out.println( "</td>" );
	out.println( "</tr>" );
	out.println( "</table>" );
    }


    /**
     * Convert from 'plus' and 'minus' to '+' and '-'
     *
     * @param poString <code>String</code> to convert.
     * @return new <code>String</code>
     */
    private String toSign( String poString ) {

	String oSign = null;

	if ( poString != null ) {
	    if (  poString.startsWith( "plus" ) ) {
		oSign = "+";
		if ( poString.length() > 4 ) {
		    oSign = oSign.concat( poString.substring( 4 ));
		}

	    } else if ( poString.startsWith( "minus" )) {
		oSign = "-";
		if ( poString.length() > 5 ) {
		    oSign = oSign.concat( poString.substring( 5 ));
		}
	    }
	}
	return oSign;
    }

    // ************************************************************ //
    // *                 Escape HTML chars                        * //
    // *                                                          * //
    // *   Probably replace this with something in a apache or    * //
    // *   java lib.                                              * //
    // *                                                          * //
    // *                                                          * //
    // ************************************************************ //

    private void writeContentChars( String poLine ) {

	poLine = this.replace( poLine, '&', "&amp;");
	poLine = this.replaceGtLtAndQuote( poLine );

	out.print( poLine );
    }

    private String replaceGtLtAndQuote( String poInputString ) {

	pcDataBuffer.setLength( 0 );
	pcDataBuffer.append( poInputString );
	int iLength = ( pcDataBuffer.length() );

	for (int i = iLength; --i >= 0; ) {
	    if (  (pcDataBuffer.charAt(i) == '<') ) {
		pcDataBuffer.deleteCharAt(i);
		// then insert escape char to the LHS
		pcDataBuffer.insert(i, "&lt;" );
	    } else if (  (pcDataBuffer.charAt(i) == '>') ) {
		pcDataBuffer.deleteCharAt(i);
		// then insert escape char to the LHS
		pcDataBuffer.insert(i, "&gt;" );
	    } else if (  (pcDataBuffer.charAt(i) == '\\') ) {
		pcDataBuffer.deleteCharAt(i);
		// then insert escape char to the LHS
		pcDataBuffer.insert(i, "&quot;" );
	    } // end if
	} // end for

	return pcDataBuffer.substring(0);
    }

    private String replace( String poInputString,
			    char pcCharToReplace,
			    String poReplacementString ) {

	pcDataBuffer.setLength( 0 );
	pcDataBuffer.append( poInputString );
	int iLength = ( pcDataBuffer.length() );

	for (int i = iLength; --i >= 0; ) {
	    if (  pcDataBuffer.charAt(i) == pcCharToReplace ) {
		pcDataBuffer.deleteCharAt(i);
		// then insert escape char to the LHS
		pcDataBuffer.insert(i, poReplacementString );
	    } // end if
	} // end for

	return pcDataBuffer.substring(0);
    }


    /**
     * Makes assumption about the gap character.
     *
     */
    int countNumberOfGaps( String poString ) {

	int index = -1;
	int iNumberOfGaps = 0;
	while ( ( index = poString.indexOf( '-', index+1 ) ) != -1 ) {
	    iNumberOfGaps++;
	}
	return iNumberOfGaps;
    }

    /**
     * Ensures the given string is the correct length.
     *
     */
    String padTo( String poString, int iNumberOfChars ) {

	int iLen = iNumberOfChars - poString.length();

	if ( iLen > padding.length ) {
	    padding = new char[ iLen ];
	    Arrays.fill( padding,
			 0,
			 iLen,
			 ' ' );
	}
	return poString.concat( String.copyValueOf( padding, 0, iLen ) );
    }


} // end class

