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

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

/**
 * <p>
 * In addition to marking up the two sequences in the same way as
 * <CODE>SimpleAlignmentStyler</CODE> this styles the markup between
 * the two sequences based on a Blast Matrix score.
 * </p>
 *
 * <p>
 * The colouring algorithm has a linear range of hues between blue and red,
 * from positive to negative scores with saturation affected by the
 * square of the difference from the mean score.
 * </p>
 *
 * <p>
 * So high scoring matches will be dark blue and large negative scoring
 * pairs will result in dark red. Scores in the middle will have a hue
 * between red and blue and be very pale.
 * </p>
 *
 * <p>
 * <CODE>BlastMatrixAlignmentStyler</CODE> reads the colourmap specified by
 * the property -DcolourMap for colouring the query and subject sequences.
 * </p>
 *
 * <p>
 * In addition will colour the consensus using a scale based on a
 * specified BLAST format scoring matrix.
 * </p>
 *
 * <p>
 * The location of the blast matrix is specified by,
 * -DblastMatrix=<filename>
 * </p>
 *
 * <p><pre>
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
 * </pre></p>
 *
 * <p>
 * This code released to the biojava project, May 2001
 * under the LGPL license.
 * </p>
 *
 * @author Cambridge Antibody Technology Group plc
 * @author Greg Cox
 * @version 1.0
 */
class BlastMatrixAlignmentStyler extends SimpleAlignmentStyler {


    private String[] aoChars = null;
    private int iNumberOfChars = 0;
    private int[][] aoScore = null;


    public BlastMatrixAlignmentStyler() {

	super( SimpleAlignmentStyler.SHOW_ALL );

	super.readColourMap();

	String oPropFilename = System.getProperty("blastMatrix");

	if( oPropFilename == null ) {
            System.err.println
		("No blast matrix specified " +
		 "with -DblastMatrix=<filename>" );

	} else {

	    try {
		this.readBlastMatrix( oPropFilename );
		this.createAlignmentColours();
	    } catch ( IOException e ) {
		System.err.println( e );
		e.printStackTrace();
	    }
	}
    }


    /**
     * Read the blast matrix into 2 data structures:
     * <ul>
     *   <li>array of chatacters in matrix ( aoChars );</li>
     *   <li>2D array of scores ( aoScore );</li>
     * </ul>
     *
     * @param poFileName - the filename of the matrix
     * @exception IOException if an error occurs
     */
    void readBlastMatrix( String poFileName )
 	throws IOException {

 	BufferedReader in = null;

 	try {
 	    File oFile = new File( poFileName );

 	    FileReader oFileReader = new FileReader( oFile );
 	    in = new BufferedReader( oFileReader );

 	    String oLine = in.readLine();

 	    boolean isFirst = true;

 	    int iCurrent = 0;

	    // 	     assuming a symmetrical layout to matrix

 	    while( oLine != null ) {

 		if ( oLine.startsWith( "#" ) || oLine.trim().equals("") )

		    { // skip comments
			oLine = in.readLine();
			continue;
 		}

 		if ( isFirst ) {
		    //     read header
 		    StringTokenizer st = new StringTokenizer( oLine );

 		    iNumberOfChars = st.countTokens();
 		    aoChars = new String[ iNumberOfChars ];
 		    aoScore = new int[iNumberOfChars][iNumberOfChars];
 		    int i = 0;
 		    while ( st.hasMoreTokens() ) {
 			aoChars[i] = st.nextToken();
 			i++;
 		    }
 		    isFirst = false;
 		} else { //read rest
		    //System.out.println( oLine );

 		    StringTokenizer st = new StringTokenizer( oLine );
 		    int i = 0;
 		    st.nextToken();
 		    while ( st.hasMoreTokens() ) {
 			aoScore[iCurrent][i] = Integer.parseInt (st.nextToken() );
 			i++;
 		    }
 		    iCurrent++;
 		}
 		oLine = in.readLine();
 	    } // end while
 	    in.close();

 	} catch ( IOException e ) {
 	    throw e;
 	} finally {
 	    try {
 		in.close();
 	    } catch ( Exception x ) {
		// 		 do nothing
 	    }
 	}

     } // end readBlastMatrix


    /**
     * Converts the blast scoring matrix into a set of colours.
     *
     */
    private void createAlignmentColours() {


	/*
	 * Find the Range of possible scores
	 *
	 */
 	int iMatchMin = 0;
 	int iMatchMax = 0;
 	for ( int i = 0; i < iNumberOfChars; i++ ) {
 	    for ( int j = 0; j < iNumberOfChars; j++ ) {

 		int iCurrent = aoScore[i][j];
 		if ( iCurrent < iMatchMin ) {
 		    iMatchMin = iCurrent;
 		} else if ( iCurrent > iMatchMax ) {
 		    iMatchMax = iCurrent;
 		}
 	    }
 	}  //end for

 	int iRange = iMatchMax - iMatchMin;

	/*
	 * For each pair calculate a normalised score, scale so between
	 * 0.0-0.33.<p>
	 * The hue is then 1.0 - this value.<BR>
	 *
	 * Higher saturation is weighted towards the extremes by taking the
	 * square of the deviation from the mean and multipling it by 20.<p>
	 *
	 * Consensus colours are put into the colour map as '!AA' '!AC' etc
	 * ie prefixed with an '!'
	 *
	 */
 	for ( int i = 0; i < iNumberOfChars; i++ ) {
 	    for ( int j = 0; j < iNumberOfChars; j++ ) {

 		int iCurrent = aoScore[i][j];
 		double dNorm = ((double)(iCurrent-iMatchMin))/(double)iRange;

		dNorm = dNorm * 0.33;
		// create saturation weighting
		float var = (float) (( dNorm- (0.33/2) ) *( dNorm- (0.33/2)));
		//	 create colour
  		Color oColor = Color.getHSBColor( (float)1.00 - (float)dNorm,
						  (float)0.0 + (var*20),
  						  (float)1.00    );
 		 // convert to hex
 		String oColourString = this.toHex( oColor );

 		String oColourClass  = getColourClass
 		    ( oColourString );

 		oColourMap.put( "!".concat( aoChars[i].concat( aoChars[j] ) ),
				oColourClass );
 	    }
 	}// end for

    }


    private String toHex( Color poColor ) {

	StringBuffer sb = new StringBuffer( 7 );
	sb.append( padHex(Integer.toHexString( poColor.getRed())) );
	sb.append( padHex(Integer.toHexString( poColor.getGreen())) );
	sb.append( padHex(Integer.toHexString( poColor.getBlue()) ));

	return sb.substring(0);
    }

    private String padHex( String poHex ) {
        if ( poHex.length() == 1 ) {
            return "0".concat( poHex );
        } else {
            return poHex;
        }
    }

    /**
     * <p>
     * Return the styles for the two aligned characters.
     * ( in the form of predefined font classes ).
     * </p>
     *
     * <p>
     * Null is acceptable value for no style.
     * </p>
     *
     * @param poFirst - the first char in the alignment
     * @param poSecond - the second char in the alignment
     * @param poStyleHolder - an array to hold the styles, [0] = first etc
     */
    public void getStyle( String poFirst, String poSecond,
			  String[] poStyleHolder ) {


	poStyleHolder[0] = (String)oColourMap.get( poFirst );
	poStyleHolder[1] = (String)oColourMap.get( poSecond );

	poStyleHolder[2] = (String)oColourMap.get( "!".concat(poFirst).concat
						   (poSecond) );
    }

} // end class
