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

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * <p>
 * Abstract implementation of <code>AlignmentStyler</code>, contains
 * utility methods for generating a set of HTML styles from a list of
 * RGB colours.
 * </p>
 *
 * <p>
 * Thus <code>getAlignmentStyles()</code> is implemented and all that
 * remains to be implemented is the <code>getStyle</code> method.
 *
 * <pre>
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
 * </pre>
 * </p>
 *
 * <p>
 * This code released to the biojava project, May 2001
 * under the LGPL license.
 * </p>
 *
 * @author Cambridge Antibody Technology Group plc
 * @author Greg Cox
 * @version 1.0
 *
 */
public abstract class AbstractAlignmentStyler implements AlignmentStyler {

    /**
     * Store the unique colours for markup.
     */
    protected HashSet oColourSet = new HashSet();


    /**
     * <p>
     * Stores mapping from a Colour to a FONT Class.
     * </p>
     *
     * <p>
     * For example:
     * <PRE>
     *
     * Key      Value
     * ---      -----
     * #000000  C1-S
     *
     * </PRE>
     * </p>
     */
    protected HashMap oColourClassMap = new HashMap();



    /**
     * The number of unique colours.
     */
    protected int iNumberOfColours = 0;

    /**
     * <p>
     * Map between Char and the Colour class.
     * </p>
     *
     * <p>
     * Eg.
     * <PRE>
     *
     * Key      Value
     * ---      -----
     * A        C1-S
     *
     * </PRE>
     * </p>
     */
    protected HashMap oColourMap = new HashMap();


    /**
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
    public String getAlignmentStyles() {

	StringBuffer sb = new StringBuffer();

	if ( oColourSet.size() == 0 ) return "";

	//	sb.append("<STYLE TYPE=\"text/css\">\n");
	//	sb.append("<!--\n");

	Iterator it = oColourSet.iterator();
	while ( it.hasNext() ) {

	    sb.append( (String)it.next() );
	}

	//	sb.append( "-->\n</STYLE>\n" );
	return sb.substring(0);
    }

    /**
     * <p>
     * Return the styles for the two aligned characters.
     * (in the form of predefined font classes).
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
    public abstract void getStyle( String poFirst, String poSecond,
				   String[] poStyleHolder );


    /**
     * Add a colour style to this Styler.
     *
     * @param poChar the char for which this colour applies.
     * @param poColour the color in hex eg 'FFA2A2' for a nice red
     *                 ( R = FF, G = A2 and B = A2 )
     */
    public void addStyle ( String poChar, String poColour ) {

	String oColourClass  = this.getColourClass
	    ( poColour );

	oColourMap.put( poChar, oColourClass );
    }

    /**
     * <p>Returns the colour class for the specified colour (in hex).
     * If one is not already defined for that colour then a new class
     * is created and returned.</p>
     *
     * <p>
     * Colour specification is R G B in hex ie
     * FF00FF is r = 255, g = 0, b = 255.
     * </p>
     *
     * @param poColour - a colour, eg 'C8FFC8'
     * @return String - the colour class, eg 'C1-S'
     */
    protected String getColourClass( String poColour ) {

	String oColourClass = (String)oColourClassMap.get( poColour );
	if ( oColourClass == null ) {
	    // otherwise create a new one
	    oColourClass  = "C" + iNumberOfColours + "-S" ;

	    StringBuffer sb = new StringBuffer( 50 );

	    sb.append ( "FONT." );
	    sb.append ( oColourClass );
	    sb.append ( "{background-color:#" );
	    sb.append ( poColour );
	    sb.append ( ";color:#000000}\n" );

	    oColourSet.add( sb.substring(0) );
	    oColourClassMap.put( poColour, oColourClass );

	    iNumberOfColours++;
	}
	return oColourClass;
    }

}
