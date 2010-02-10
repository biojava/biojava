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

import java.util.List;


/**
 * Interface for a factory method to return a list of
 * 0, 1 or more <code>DatabaseURLGenerator</code> instances.
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
public interface URLGeneratorFactory {

    /**
     * Returns a list of 0, 1 or more  <code>DatabaseURLGenerator</code>
     * The first is used to link from the id in the summary table, all are
     * used as a list of links in the detail section.
     *
     * @return a <code>List</code>
     */
    public List getDatabaseURLGenerators();
}
