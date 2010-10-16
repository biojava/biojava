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
 * Interface for specifying whether a particular pair
 * of residues/bases should be coloured.
 * </p>
 *
 * <p>
 * This can be independent of the method that
 * chooses which colour.
 * </p>
 *
 * <p>
 * Example usage: highlight mismatches only.
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
public interface ColourCommand {

    /**
     * Returns true if the alignment pair should be coloured
     * else false.
     *
     * @param poFirst  - the first character
     * @param poSecond - the second character.
     * @return boolean - true if colour, else false.
     */
    public boolean isColoured( String poFirst, String poSecond );
}
