// BiblioSubject.java
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
 * It represents the topic of the content of the cited resource.
 * It can be expressed in one or more ways.
 * </p>
 *
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @version $Id$
 * @since 1.3
 */

public class BiblioSubject {

    /**
     * The keywords are usually (but not limited to) one word long.
     * They are not controlled by any vocabulary.
     */
    public Hashtable keywords = new Hashtable();

    /**
     * The subject headings usually come from standard lists
     * such as <em>Sears List of Subject Headings</em>,
     * or <em>Library of Congress Subject Headings (LCSH)</em>.
     * This specification does not specify what list to use but implementors
     * are advised to provide a controlled vocabulary for the list that is used,
     * and to specify the source of subject headings in
     * {@link #subjectHeadingsSource} field. The name of such vocabulary
     * should be equal to {@link BibRefSupport#SUBJECT_HEADINGS}.
     */
    public Hashtable subjectHeadings = new Hashtable();

    /**
     * The source of {@link #subjectHeadings subject headings}.
     * For example:
     * <dl>
     *   <dt> <tt>SEARS</tt>
     *   <dd> for <em>Sears List of Subject Headings</em>
     *   <dt> <tt>LCSH</tt>
     *   <dd> for <em>Library of Congress Subject Headings (LCSH)</em>
     *   <dt> <tt>MeSH</tt>
     *   <dd> for <em>MEDLINE Mesh Terms</em>
     * </dl>
     */
    public String subjectHeadingsSource;

    /**
     * <p>
     * Classification code (call number) is usually either Dewey decimal or
     * Congress classification. But this specification does not prescribe it.
     * </p>
     *
     * <p>
     * Note that the classification codes are unique (unlike some subject headings).
     * Therefore, they can be even expressed as identifiers using the same notation
     * as used for the {@link BibRef#identifier citation identifiers}.
     * </p>
     */
    public Hashtable codes = new Hashtable();

}
