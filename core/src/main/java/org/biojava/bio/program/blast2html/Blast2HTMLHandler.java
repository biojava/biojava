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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

/**
 * Takes a SAX event stream and a HTMLRenderer to produce
 * a HTML Blast like program report.
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
public class Blast2HTMLHandler extends DefaultHandler  {

    /**
     * Variables to hold data while parsing.
     *
     */
    private StringBuffer sb = new StringBuffer();

    private String oProgram;
    private String oVersion;
    private String oQuery;
    private String oDatabase;

    private HitSummary oHitSummary = new HitSummary();
    private HitId      oHitId      = new HitId();
    private HitDescription oDesc   = new HitDescription();
    {
	oHitSummary.oHitId = oHitId;
	oHitSummary.oDesc  = oDesc;
    }

    private DetailHit oDetailHit = new DetailHit();
    private HSP       oHSP       = new HSP();
    private HSPSummary oHSPSummary = new HSPSummary();
    private BlastLikeAlignment oAlignment = new BlastLikeAlignment();
    private Sequence oQuerySeq = new Sequence();
    private Sequence oHitSeq = new Sequence();

    {
	oDetailHit.oHitId = oHitId;
	oDetailHit.oDesc  = oDesc;
	oHSP.oHSPSummary = oHSPSummary;
	oHSP.oAlignment   = oAlignment;
	oAlignment.oQuerySeq = oQuerySeq;
	oAlignment.oHitSeq = oHitSeq;

    }

    private String oRawOutput = null;

    // Flow control flags
    private boolean inCollection = false;
    private boolean firstSummary = true;
    private boolean firstDetail = true;

    /**
     * The Class to render the HTML
     */
    private HTMLRenderer oRenderer = null;



    /**
     * A content handler for rendering blast like outputs into
     * HTML.
     *
     * @param poRenderer <code>HTMLRenderer</code> - a configured
     *                   HTMLRenderer.
     */
    public Blast2HTMLHandler( HTMLRenderer poRenderer ) {

	if ( poRenderer == null ) {
	    throw new IllegalArgumentException
		( "HTMLRenderer cannot be null" );
	}
	oRenderer = poRenderer;
    }



    // ************************************************************ //
    // ****            ContentHandler overrides                **** //
    // ************************************************************ //

    /**
     * This is called when an element is entered. That is,
     * the parser has met the first tag of the tag pair.
     *
     * @param poNameSpace <code>String</code> - the name space.
     * @param poElementName <code>String</code> - the local name of the tag.
     * @param poQName  <code>String</code> - the fully qualified name
     *                                       with prefix
     * @param poAtts an <code>Attributes</code> - the tag attributes.
     * @exception SAXException if an error occurs
     */
    public void startElement ( String poNameSpace, String poElementName,
			       String poQName, Attributes poAtts)
	throws SAXException {

	if ( poElementName.equals( "BlastLikeDataSetCollection" ) ) {
	    inCollection = true; // only checking to do - assume once
	                         // inside collection it follows DTD
	}

	if ( inCollection ) {

	    if ( poElementName.equals( "BlastLikeDataSet" ) ) {

		oProgram = poAtts.getValue( "program" );
		oVersion = poAtts.getValue( "version" );
		sb.setLength( 0 );
	    } else if ( poElementName.equals( "HitSummary" ) ) {

		oHitSummary.score       = poAtts.getValue( "score" );
		oHitSummary.expectValue = poAtts.getValue( "expectValue" );

		oHitSummary.numberOfHSPs = poAtts.getValue( "numberOfHSPs" );
		oHitSummary.numberOfContributingHSPs = poAtts.getValue
		    ( "numberOfContributingHSPs" );
		oHitSummary.smallestSumProbability = poAtts.getValue
		    ( "smallestSumProbability" );
		oHitSummary.readingFrame = poAtts.getValue( "readingFrame" );
		oHitSummary.numberOfDomains = poAtts.getValue
		    ( "numberOfDomains" );

		if ( firstSummary ) {
		    oRenderer.startSummaryTable( oHitSummary );
		    firstSummary = false;
		}

		//
		// WHAT happens if there is more than one HSPCollection
		// per hit - probably won't work.
		//
		//
	    } else if ( poElementName.equals( "HSPCollection" ) ) {
		oRenderer.writeCurrentDetail( oDetailHit );
	    } else if ( poElementName.equals( "Hit" ) ) {

		if ( firstDetail ) {
		    oRenderer.startDetailTable();
		    firstDetail = false;
		}

		oDetailHit.sequenceLength = poAtts.getValue
		    ( "sequenceLength" );

	    } else if ( poElementName.equals( "HSPSummary" ) ) {

		oHSPSummary.percentageIdentity  = poAtts.getValue
		    ( "percentageIdentity" );
		oHSPSummary.score       = poAtts.getValue( "score" );
		oHSPSummary.expectValue = poAtts.getValue( "expectValue" );
		oHSPSummary.alignmentSize  = poAtts.getValue
		    ( "alignmentSize" );
		oHSPSummary.numberOfIdentities = poAtts.getValue
		    ( "numberOfIdentities" );

		oHSPSummary.hitStrand       = poAtts.getValue( "hitStrand" );
		oHSPSummary.queryStrand     = poAtts.getValue( "queryStrand" );
		oHSPSummary.queryFrame      = poAtts.getValue( "queryFrame" );
		oHSPSummary.hitFrame        = poAtts.getValue( "hitFrame" );

		oHSPSummary.numberOfPositives  = poAtts.getValue
		    ( "numberOfPositives" );
		oHSPSummary.percentagePositives= poAtts.getValue
		    ( "percentagePositives" );
		oHSPSummary.pValue          = poAtts.getValue( "pValue" );
		oHSPSummary.sumPValues      = poAtts.getValue( "sumPValues" );
		oHSPSummary.numberOfGaps    = poAtts.getValue( "numberOfGaps");


	    } else if ( poElementName.equals( "HitId" ) ) {
		oHitId.id = poAtts.getValue( "id" );
		oHitId.metaData = poAtts.getValue( "metaData" );

	    } else if ( poElementName.equals( "HitDescription" ) ) {
		sb.setLength( 0 );
	    } else if ( poElementName.equals( "QuerySequence" ) ) {
		oAlignment.oQuerySeq.startPosition =
		    poAtts.getValue( "startPosition" );
		oAlignment.oQuerySeq.stopPosition  =
		    poAtts.getValue( "stopPosition" );
		sb.setLength( 0 );
	    } else if ( poElementName.equals( "MatchConsensus" ) ) {
		sb.setLength( 0 );
	    } else if ( poElementName.equals( "HitSequence" ) ) {
		oAlignment.oHitSeq.startPosition
		    = poAtts.getValue( "startPosition" );
		oAlignment.oHitSeq.stopPosition
		    = poAtts.getValue( "stopPosition" );
		sb.setLength( 0 );
	    }  else if ( poElementName.equals( "RawOutput" ) ) {
		sb.setLength( 0 );
	    }

	} // end inCollection

    }

    /**
     * Called when the end of an element is reached.
     *
     * @param poNameSpace a <code>String</code> - the name space.
     * @param poElementName a <code>String</code> - the local element name.
     * @param poQName a <code>String</code> value - the qualified element name.
     */
    public void endElement ( String poNameSpace,
			     String poElementName,
			     String poQName ) {


	if ( poElementName.equals( "Header" ) ) {
	    this.getQueryIdAndDatabase();
	    oRenderer.writeTitleAndHeader( oProgram, oVersion,
				      oQuery, oDatabase );

	} else 	if ( poElementName.equals( "HitDescription" ) ) {
	    oDesc.hitDescription = sb.substring(0);

	} else 	if ( poElementName.equals( "Summary" ) ) {
	    oRenderer.endSummaryTable();
	} else 	if ( poElementName.equals( "HitSummary" ) ) {
	    oRenderer.writeCurrentSummary( oHitSummary );
	} else 	if ( poElementName.equals( "Detail" ) ) {
	    oRenderer.endDetailTable();
	} else 	if ( poElementName.equals( "RawOutput" ) ) {
	    oRawOutput  = sb.substring(0);
	} else 	if ( poElementName.equals( "HSPSummary" ) ) {
	    oHSPSummary.rawOutput = oRawOutput;
	} else if ( poElementName.equals( "QuerySequence" ) ) {
	    oAlignment.oQuerySeq.seq = sb.substring(0);
	} else if ( poElementName.equals( "MatchConsensus" ) ) {
	    oAlignment.oConsensus = sb.substring(0);
	} else if ( poElementName.equals( "HitSequence" ) ) {
	    oAlignment.oHitSeq.seq = sb.substring(0);
	} else if ( poElementName.equals( "HSP" ) ) {
	    oRenderer.writeCurrentHSP( oHSPSummary, oAlignment );
	} else 	if ( poElementName.equals( "BlastLikeDataSetCollection" ) ) {
	    this.reInit();
	}

    } // end

    /**
     * Describe <code>characters</code> method here.
     *
     * @param charBuffer - character array containing data.
     * @param start - the start position of relavent chars in passes array
     * @param length - the stop position of relavent chars in passes array
     */
    public void characters( char[] charBuffer, int start, int length) {

	// note this may be called more than once for a particular tag
	// ie loaded in chunks.
	// So the correct way to handle this stuff is to buffer contents
	// and deal with buffer when endElement is called.

	sb.append( charBuffer, start, length );
    }

    // ************************************************************ //
    // ****                 Utility functions                  **** //
    // ************************************************************ //

    /**
     * Re-initializes state, called between parsings.
     */
    void reInit() {
	inCollection = false;
	firstSummary = true;
	firstDetail = true;
	oRawOutput = null;
	sb.setLength( 0 );
    }


    /**
     * Parses out query and database id's from, the rawoutput.
     *
     * Changes in the Sax event generator may have made this redundant.
     *
     */
    void getQueryIdAndDatabase() {

	BufferedReader reader = new BufferedReader
            ( new StringReader( sb.substring(0) ) );

        String oLine;

        try {
            while (( oLine = reader.readLine() ) != null) {

		if ( oLine.startsWith( "Query=" ) ) {

		    int index = oLine.indexOf( "=" );
		    oQuery = oLine.substring( index+1 ).trim();
		    continue;
		}
		if ( oLine.startsWith( "Database:" ) ) {

		    int index = oLine.indexOf( ":" );
		    oDatabase = oLine.substring( index+1 ).trim();
		    break;
		}
            }
        }
        catch ( IOException e ) {
	    printError( e );
        }
    }


    /**
     * Print an error to System.err
     *
     * @param e an <code>Exception</code>
     */
    void printError( Exception e ) {

	System.out.println( e.getMessage() );
	e.printStackTrace();
    }



} // end class

// ************************************************************ //
// ****        Simple holder classes to hold temporary     **** //
// ****                 values  during parsing             **** //
// ************************************************************ //


class DetailHit {

    public String sequenceLength;
    public HitDescription oDesc;
    public HitId  oHitId;

}

class HitDescription {

    String hitDescription;
}

class HSP {

    public HSPSummary oHSPSummary;
    public BlastLikeAlignment oAlignment;
}

class HSPSummary {

    public String score;
    public String expectValue;
    public String numberOfIdentities;
    public String alignmentSize;
    public String percentageIdentity;

    public String hitStrand;
    public String queryStrand;

    public String queryFrame;
    public String hitFrame;

    public String numberOfPositives;
    public String percentagePositives;
    public String pValue;
    public String sumPValues;
    public String numberOfGaps;

    public String rawOutput;
}

class BlastLikeAlignment {

    public Sequence oQuerySeq;
    public Sequence oHitSeq;
    public String   oConsensus;
}

class Sequence {

    public String seq;
    public String startPosition;
    public String stopPosition;
}


class HitSummary {

    public String score;
    public String expectValue;
    public HitId  oHitId;
    public HitDescription oDesc;

    public String numberOfHSPs;
    public String numberOfContributingHSPs;
    public String smallestSumProbability;
    public String readingFrame;
    public String numberOfDomains;

}

class HitId {

    public String id;
    public String metaData;
}



