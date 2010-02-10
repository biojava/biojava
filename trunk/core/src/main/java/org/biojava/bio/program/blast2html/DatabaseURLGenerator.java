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

import java.util.Properties;
/**
 * Takes a database ID and some configuration properties
 * ( such as base URL ) and returns either a URL or 
 * a full anchor tag.
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
 * @version 1.0
 */
public interface DatabaseURLGenerator {

    /**
     * Returns a string representation of a URL to the 
     * specified ID. Used in summary section.
     *
     * @param poID a database ID
     * @param poOptions <code>Properties</code> - any options needed
     * @return a <code>String</code> value
     */
    String toURL( String poID, Properties poOptions );
    /**
     * Returns a full <a href=_____>retrieve item</a> anchor
     * for the given database id. Used in detail section.
     *
     * @param poID a <code>String</code> value
     * @param poOptions a <code>Properties</code> value
     * @return a <code>String</code> value
     */
    String toLink( String poID, Properties poOptions );
}
