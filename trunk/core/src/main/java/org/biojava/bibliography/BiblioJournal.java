// BiblioJournal.java
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
 * A class describing journals. The citations referring to the journal articles
 * have a reference to this class. There are only few explicit attributes defined,
 * the rest are accessible using {@link #properties dynamic properties}.
 * </p>
 *
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @version $Id$
 * @since 1.3
 */

public class BiblioJournal {

    /**
     * Additional properties used when the explicit attributes
     * are not sufficient.
     */
    public Hashtable properties = new Hashtable();

    /**
     * A journal title.
     * The list of available titles can be provided using a controlled vocabulary,
     * taken, for example, from
     * <a href="http://www.ncbi.nlm.nih.gov/PubMed/fulltext.html"> MEDLINE Journals</a>.
     * Such controlled vocabulary should be named {@link BibRefSupport#JOURNAL_TITLES}.
     */
    public String name;

    /**
     * <p>
     * A standard number for journals.
     * </p>
     *
     * <p>
     * Be aware, however, that in the real world even this attribute may change over time.
     * Therefore, it is not suitable as a primary unique identifier for journals.
     * </p>
     */
    public String issn;

    /**
     * <p>
     * An abbreviation of the journal title.
     * </p>
     *
     * <p>
     * Note that some repositories use more abbreviation systems. For such cases,
     * use {@link #properties dynamic properties} for additional abbreviations.
     * </p>
     *
     * <p>
     * An example for biological journals is in
     * <a href="http://arachne.prl.msu.edu/journams/">Biological Journals and Abbreviations</a>.
     * </p>
     *
     * A controlled vocabulary with abbreviation should be named
     * {@link BibRefSupport#JOURNAL_ABBREV}.
     */
    public String abbreviation;

}
