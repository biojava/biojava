// BiblioJournalArticle.java
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

/**
 * <p>
 * It represents a journal article.
 * </p>
 *
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @version $Id$
 * @since 1.3
 */

public class BiblioJournalArticle
    extends BiblioArticle {

  /**
   * The journal volume.
   */
    public String volume;

  /**
   *  The issue of the journal.
   */
    public String issue;

  /**
   * Suplement.
   */
    public String issueSupplement;

  /**
   * The journal that this article is in.
   */
    public BiblioJournal fromJournal;

}
