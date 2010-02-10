// BiblioCriterion.java
//
//    senger@ebi.ac.uk
//    April 2001
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
 * The criteria define how the matching or ordering should be done
 * during queries.
 *
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @version $Id$
 * @since 1.3
 */

public class BiblioCriterion {
  /**
   * A query criterion.
   */
    public static final int QUERY_CRITERION = 0;

  /**
   * A sort criterion.
   */ 
    public static final int SORT_CRITERION  = 1;

    /**
     * <p>
     * Each Criterion is identified by its name.
     * A list of criteria names is used in methods for querying and sorting
     * (see {@link BibRefQuery} interface).
     * </p>
     *
     * <p>
     * The implementations are advised to use descriptive names.
     * For example, the names for matching can be:
     * <pre>
     *     match all words
     *     match any word
     *     case insensitive
     *     case sensitive
     *     partial word match
     *     full word match
     * </pre>
     * and the names for ordering can be:
     * <pre>
     *      ascending
     *      descending
     * </pre>
     * Another example of how to use Criteria is to allow regular expressions in queries.
     * Not every implementation is supposed to have the capability of matching by regular
     * expressions but those who have can specify (and document), for example, criterion
     * with name <tt>regular expression</tt>.
     * </p>
     */
    public String name;

    /**
     * The criteria can be used for defining rules for matching
     * (type {@link #QUERY_CRITERION}), or for ordering (type {@link #SORT_CRITERION}).
     */
    public int type = QUERY_CRITERION;

    /**
     * <p>
     * A list of other criteria names that this criterion is mutually exclusive with.
     * </p>
     *
     * <p>
     * For example, a sort criterion <tt>ascending</tt> will probably have
     * <tt>descending</tt> in this list.
     * </p>
     */
    public String[] mutuallyExclusiveWith;

    /**
     * A name of a repository subset which this criterion is valid/used for.
     * @see BiblioEntryStatus#repositorySubset
     */
    public String forSubset;

}
