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


/**
 * <p>
 * Class to do simple HTML colouring of sequence alignments.
 * </p>
 *
 * <p>
 * For any particular alignment position, the colouring depends only on the two
 * characters at that position. The decision is made in two stages:-
 * <ol>
 *   <li>Whether to colour or not in the <CODE>ColourCommand.isColoured()</CODE>
 *    method</li>
 *   <li>what the colours are in the <CODE>AlignmentStyler.getStyle()</CODE>
 *    method</li>
 * </ol>
 * </p>
 *
 * <p>
 * This allows simple choices to highlight the mismatches, the identities or
 * simply colour up everything.
 * </p>
 *
 * <p>
 * The current implemention of step 2. is a simple colour lookup.
 * </p>
 *
 * <p>
 * Limitations:
 * </p>
 *
 * <p>
 * As the FONT styles need to be defined before being used in the HTML,
 * it means the all colours to be used have to calculated up front.
 * </p>
 *
 * <p>
 * The position in the alignment is not passed in so position dependent
 * colouring ( say against a HMM profile ) would be involve either
 * changing the interfaces or implementing them such that they knew the
 * required information via another route.
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
public class AlignmentMarker {

    /**
     * Resuable StringBuffers for the markup of the alignments
     */
    private StringBuffer[] markedUp  = new StringBuffer[ 3 ];
    {
	markedUp[0] = new StringBuffer( 150 );
	markedUp[1] = new StringBuffer( 150 );
	markedUp[2] = new StringBuffer( 150 );
    }


    /**
     * Holds the current style for query, hit and consensus
     */
    private String[] oCurrentStyle = new String[3];
    {
	oCurrentStyle[0] = null;
	oCurrentStyle[1] = null;
	oCurrentStyle[2] = null;
    }


    /**
     * Holds the NewStyles for query, hit and consensus
     */
    private String[]  oNewStyle = new String[ 3 ];

    /**
     * Single method interface for deciding whether to colour
     * a particular alignment pair.
     */
    private ColourCommand oColourCommand;

    /**
     * Class that determines the style for each char
     */
    private AlignmentStyler oStyler;

    /**
     * Creates a new <code>AlignmentMarker</code> instance.
     *
     * @param poColourCommand - controls whether a particular alignment pair
     *                          should be coloured
     * @param poStyler - specifies all possible styles and where each style
     *                   is used
     */
    public AlignmentMarker( ColourCommand poColourCommand,
			    AlignmentStyler  poStyler ) {

	oColourCommand = poColourCommand;
	oStyler  = poStyler;
    }


    /**
     * <p>
     * Delegate to the AlignmentStyler
     * </p>
     *
     * <p>
     * Returns a fragment of HTML that defines the FONT
     * styles to be used in the alignment markup.
     * </p>
     *
     * <p>
     * For example:
     * <PRE>
     * FONT.C2-S{background-color:#FFFC50;color:#000000}
     * FONT.C4-S{background-color:#FC50FF;color:#000000}
     * FONT.C3-S{background-color:#FF7272;color:#000000}
     * FONT.C0-S{background-color:#50FF78;color:#000000}
     * FONT.C1-S{background-color:#FFCA50;color:#000000}
     * FONT.C5-S{background-color:#A5A5FF;color:#000000}
     * </PRE>
     * </p>
     *
     * @return String - the HTML
     */
    String getAlignmentStyles() {

	 return oStyler.getAlignmentStyles();
     }

    /**
     * Takes three sequences ( of the same length ) as strings and
     * returns three strings with colour markup in HTML.
     *
     * @param poAlignment - three strings representing a pairwise alignment
     *                      query, hit, consensus
     */
    void alignment2HTML( String[] poAlignment ) {

	if ( poAlignment == null || poAlignment.length != 3
	     || poAlignment[0].length() != poAlignment[1].length()
	     || poAlignment[0].length() != poAlignment[2].length() ) {

	    System.err.println( "-->" + poAlignment[0] + "<--" );
	    System.err.println( "-->" + poAlignment[1] + "<--" );
	    System.err.println( "-->" + poAlignment[2] + "<--" );

	    throw new IllegalArgumentException
		( "Only accept array of three strings, all of same length" );
	}

	if ( oStyler == null ) { // then no styles
	    return;
	}

	// Initialise

	markedUp[0].setLength(0);
	markedUp[1].setLength(0);
	markedUp[2].setLength(0);

	oCurrentStyle[0] = null;
	oCurrentStyle[1] = null;
	oCurrentStyle[2] = null;

	// For each position

	for( int i= 0, n = poAlignment[0].length(); i < n ; i++) {

	    String oFirst  = String.valueOf( poAlignment[0].charAt( i ) );
	    String oSecond = String.valueOf( poAlignment[1].charAt( i ) );

	    // Decide whether to apply colours.
	    if ( oColourCommand.isColoured( oFirst, oSecond ) ) {
		oStyler.getStyle( oFirst, oSecond,
				  oNewStyle );

	    } else {
		// put in a loop if want to change the number
		oNewStyle[0] = null;
		oNewStyle[1] = null;
		oNewStyle[2] = null;
	    }

	    this.applyStyles( oCurrentStyle, oNewStyle, markedUp );

	    markedUp[0].append( oFirst );
	    markedUp[1].append( oSecond );
	    markedUp[2].append( poAlignment[2].charAt( i ) );

	    System.arraycopy( oNewStyle, 0, oCurrentStyle,
			      0, oNewStyle.length );

	} // for each char in sequence

	this.flushStyles( oCurrentStyle, markedUp );

	poAlignment[0] = markedUp[0].substring(0);
	poAlignment[1] = markedUp[1].substring(0);
	poAlignment[2] = markedUp[2].substring(0);
    }


    /**
     * Simple utility function to call applyStyle on seq1, markup &
     * seq2.
     *
     * @param poCurrentStyle a <code>String[]</code>
     * @param poNewStyle a <code>String[]</code>
     * @param poMarkedUp a <code>StringBuffer[]</code>
     */
    private void applyStyles( String[] poCurrentStyle,
			      String[] poNewStyle,
			      StringBuffer[] poMarkedUp ) {

	for ( int i=0, n = poCurrentStyle.length; i < n ; i++ ) {
	    this.applyStyle( poCurrentStyle[i], poNewStyle[i], poMarkedUp[i] );
	}
    }

    /**
     * <p>
     * Apply the new style to the output.
     * </p>
     *
     * <p>
     * Takes care of runs of the same style and
     * changing styles.
     * </p>
     *
     * @param poCurrentStyle  <code>String</code> - the current style
     * @param poNewStyle      <code>String</code> - the new style
     * @param poOutput        <code>StringBuffer</code> - the styled output
     */
    private void applyStyle( String poCurrentStyle,
			     String poNewStyle,
			     StringBuffer poOutput ) {
	//
	// If new style is null we just need to close the old style
	// if it was different.
	//
	//   else if the new style is different from the old
	//     start a new style.
	//
	if ( poNewStyle == null ) {

	    if ( poCurrentStyle != null ) {
		poOutput.append( "</FONT>" );
	    }
	} else if ( !poNewStyle.equals( poCurrentStyle ) ) {
	    if ( poCurrentStyle != null ) {
		poOutput.append( "</FONT>" );
	    }
	    // start a new one
	    poOutput.append("<FONT CLASS=" + poNewStyle + ">");
	}
    }

    /**
     * If the last style is not null then close it.
     *
     */
    private void flushStyles( String[] poCurrentStyle,
                              StringBuffer[] poMarkedUp ) {

	for ( int i=0, n = poCurrentStyle.length; i < n ; i++ ) {
	    if ( poCurrentStyle[i] != null ) { //then close
		poMarkedUp[i].append( "</FONT>" );
	    }
	}
    }
}


