// BiblioScope.java
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
 * It represent an extent or scope of the content of the cited resource.
 * It can include spatial location (a place name or geographic co-ordinates),
 * temporal period (a period label, date, or date range), or both.
 * Finally, it can have additional dynamic properties such as jurisdiction.
 * </p>
 *
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @version $Id$
 * @since 1.3
 */
public class BiblioScope {

    /**
     * It may contain additional properties related to the extend or
     * scope of the cited resource. For example:
     * <pre>
     *    jurisdiction => government
     * </pre>
     * As with all properties, it is recommended to use a controlled vocabulary.
     * Such vocabulary should be named {@link BibRefSupport#ATTR_SCOPE}.
     */
    public Hashtable properties = new Hashtable();

    /**
     * <p>
     * It defines a spatial location of the cited resource.
     * </p>
     *
     * <p>
     * This specification does not suggest any rules for representing geographical
     * names but implementations may consider 
     * <a href="http://shiva.pub.getty.edu/tgn_browser/">Getty Thesaurus of Geographic Names</a>, or
     * <a href="http://lcweb.loc.gov/marc/countries/">MARC lists of countries</a> and
     * <a href="http://lcweb.loc.gov/marc/geoareas/">MARC list of geographical areas</a>.
     * </p>
     */
    public String spatialLocation;

    /**
     * It defines temporal period of the cited resource.
     * It may be a period label, date, or date range.
     */
    public String temporalPeriod;

}

