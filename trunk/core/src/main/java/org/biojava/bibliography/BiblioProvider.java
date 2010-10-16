// BiblioProvider.java
//
//    senger@ebi.ac.uk
//    March 2001
//

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
package org.biojava.bibliography;

import java.util.Hashtable;

/**
 * <p>
 * This class and its sub-classes define active participants of the process
 * of creation and dissemination of the bibliographic resources. The most
 * obvious examples are authors, but it includes also publishers and other
 * contributors.
 * </p>
 *
 * <p>
 * The participants can be people, organizations, or even software services
 * (mainly used for new digital resources). This specification does not provide
 * any unique identifiers for these providers (but does not preclude having such
 * identifiers as {@link #properties dynamic properties}).
 * </p>
 *
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @version $Id$
 * @since 1.3
 */

public class BiblioProvider {

    /**
     * Additional properties used when the explicit attributes of the
     * sub-classes are not sufficient.
     */
    public Hashtable properties = new Hashtable();

}
