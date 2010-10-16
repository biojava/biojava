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
package org.biojava.bio.program.sax;



/**
 * Support for SAX2 configuration of namespace support
 * <p>
 * Copyright &copy; 2000 Cambridge Antibody Technology Group plc.
 *
 * <p>
 * Primary author -<ul>
 * <li>Simon Brocklehurst (CAT)
 * </ul>
 * Other authors  -<ul>
 * <li>Tim Dilks          (CAT)
 * <li>Colin Hardman      (CAT)
 * <li>Stuart Johnston    (CAT)
 *</ul>
 *
 * This code was first released to the biojava.org project, July 2000.
 *
 * @author Cambridge Antibody Technology Group plc (CAT)
 * @version 1.0
 *
 */
interface NamespaceConfigurationIF {

    /**
     * Support SAX2 configuration of namespace support of parser.
     */
    boolean getNamespaces();
    /**
     * Support SAX2 configuration of namespace support of parser.
     */
    boolean getNamespacePrefixes();
    /**
     * Gets the URI for a namespace prefix, given that prefix,
     * or null if the prefix is not recognised.
     *
     * @param poPrefix a <code>String</code> The namespace prefix.
     */
    String getURIFromPrefix(String poPrefix);


}

