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

import java.io.FileInputStream;
import java.util.Iterator;
import java.util.Map;
import java.util.Properties;

/**
 * Simple implementation for specifying markup styles.
 * Has 3 modes of operation: SHOW_ALL, SHOW_SAME & SHOW_DIFF.<p>
 *
 * SHOW_ALL  - returns the default style for all given residues.
 * SHOW_SAME - only returns a markup style if the <B>styles</B> for both
 *             characters are the same.
 * SHOW_DIFF - only returns a markup style if the <B>styles</B> for both
 *             are different.
 *
 * Styles can be easily defined in two ways.<BR>
 * 
 * 1. Add each style by calling <CODE>addStyle( poChar, poColour )</CODE>
 *    For example,
 * <CODE>
 *    	String oRed = "FFA2A2";
 *	oStyler.addStyle( "-", oRed );
 *	oStyler.addStyle( "N", oRed );
 *	oStyler.addStyle( "A", oRed );
 *	oStyler.addStyle( "T", oRed );
 *	oStyler.addStyle( "C", oRed );
 *	oStyler.addStyle( "G", oRed );
 * </CODE><p>
 *
 * 2. Alternatively the styles could be specified in a java properties file
 *    and loaded by calling <CODE>readColourMapFromProperties( poFilename )</CODE>,
 *    or <CODE>readColourMap()</CODE> and setting the system property 'colourMap'
 *    to the correct filename. <BR>
 *
 *    This file should be in java properties format, mapping 
 *    characters to colours, specified in HEX RGB.
 * 
 * For example:
 * <PRE>
 * # set everything red
 * - = FFA2A2
 * N = FFA2A2
 * A = FFA2A2
 * T = FFA2A2
 * C = FFA2A2
 * G = FFA2A2
 * </PRE> 
 *
 * Note this is simply character based, so if you want to colour gaps then
 * you need to specify a colour for the gap character.<p>
 *
 * If no colour is specified for a character then it is uncoloured.
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
 * @version 1.0
 *
 */
public class SimpleAlignmentStyler extends AbstractAlignmentStyler {

    /**
     * Return default styles
     */
    public static int SHOW_ALL = 0;
    /**
     * Only return if the two colour classes for
     * query and subject are the same
     */
    public static int SHOW_SAME = 1;
    /**
     * As NORMAL except only return if the two colour classes for
     * query and subject are the different
     */
    public static int SHOW_DIFF = 2;

    private int iStyle = 0;

    /**
     * Creates a new <CODE>SimpleAlignmentStyler</CODE> instance.<p>
     *
     * The int flag should be one of SimpleAlignmentStyler.SHOW_ALL,
     * SimpleAlignmentStyler.SHOW_SAME or
     * SimpleAlignmentStyler.SHOW_DIFF.
     *
     * @param piStyle (one of SimpleAlignmentStyler.SHOW_SAME or SimpleAlignmentStyler.SHOW_DIFF).
     * @throws IllegalArgumentException - if style not one of allowed values
     */
    public SimpleAlignmentStyler( int piStyle ) {

	if ( piStyle != SHOW_DIFF &&  piStyle != SHOW_ALL &&
	     piStyle != SHOW_SAME ) {
	    throw new IllegalArgumentException
		( "Style flag not one of SimpleAlignmentStyler.SHOW_DIFF, " +
		  " SHOW_ALL or SHOW_SAME" );
	}
	iStyle = piStyle;
    }

    /**
     * Setup styles from java property file.
     *
     * @param poFileName - the file name of the property file.
     */
    protected void readColourMapFromProperties( String poFileName ) {

	// load in properties
	Properties oColourProps = new Properties();
        
	try{
	    FileInputStream fis = new FileInputStream( poFileName );
	    oColourProps.load(fis);
	    fis.close();
	} catch (java.lang.Exception e) {
        
	    System.out.println("Failed to read properties file: " +
			       poFileName );
	    System.out.println( e.getMessage() );
	    e.printStackTrace();
	    return;
	}

	    for (Iterator i=oColourProps.entrySet().iterator(); 
		 i.hasNext(); ) {
		Map.Entry e = (Map.Entry) i.next();

		String oColourClass  = this.getColourClass
		    ( (String)e.getValue() );

		oColourMap.put( e.getKey(), oColourClass );
	    }
    }

    /**
     * Read the the properties file that specifies the character/colour mapping.
     * The location of the property file is specified by the system property
     * 'colourMap'.
     *
     */
    protected void readColourMap() {

	String oPropFileName = System.getProperty("colourMap");
    
	if( oPropFileName == null ) {
            System.err.println
		("No ColourMap preference file specified " +
		 "with -DcolourMap=<filename>" );

        } else {
	    this.readColourMapFromProperties( oPropFileName );
	}
    }

    /**
     * Returns the styles for the two aligned characters in the form
     * of predefined font classes.<p>
     *
     * Null is acceptable value for no style.
     *
     * @param poFirst - the first char in the alignment
     * @param poSecond - the second char in the alignment
     * @param poStyleHolder - an array to hold the styles, [0] = first etc
     */
    public void getStyle( String poFirst, String poSecond,
			  String[] poStyleHolder ) {

 	poStyleHolder[0] = (String)oColourMap.get( poFirst  );
 	poStyleHolder[1] = (String)oColourMap.get( poSecond );

	if ( iStyle == SimpleAlignmentStyler.SHOW_SAME ) {
	    
	    if ( poStyleHolder[0] != poStyleHolder[1] ) {
		poStyleHolder[0] = null;
		poStyleHolder[1] = null;
	    }
	} else if ( iStyle == SimpleAlignmentStyler.SHOW_DIFF) {

	    if ( !(poStyleHolder[0] != poStyleHolder[1]) ) {
		poStyleHolder[0] = null;
		poStyleHolder[1] = null;
	    }
	} 
    }

} // end class
