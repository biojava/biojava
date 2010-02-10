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
 * Interface for specifying both the overall set of possible
 * markup styles and the particular style for any pair
 * of residues/bases.
 * </p>
 *
 * <p>
 * This can be independent of the choice of whether to apply a style
 * to a pair, by using the <code>ColourCommand</code> interface.
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
 * @version 1.0
 *
 */
interface AlignmentStyler {


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
     * </PRE></p>
     *
     * @return String - the HTML
     */
    String getAlignmentStyles();

    /**
     * <p>
     * Modifies in place the styles for the two aligned characters and
     * the consensus ( in the form of predefined font classes ).
     * </p>
     *
     * <p>
     * Null is acceptable value for no style.
     * </p>
     *
     * @param poFirst - the first char in the alignment
     * @param poSecond - the second char in the alignment
     * @param poStyleHolder - an array to hold the styles, 
     *                       [0] = first, 
     *                       [1] = second,
     *                       [2] = markup
     *                        
     */
    void getStyle( String poFirst, String poSecond, String[] poStyleHolder );

}
