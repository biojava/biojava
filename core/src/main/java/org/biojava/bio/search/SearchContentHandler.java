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

package org.biojava.bio.search;


/**
 * <code>SearchContentHandler</code> is a notification interface for
 * objects which listen to search stream parsers. This is applicable
 * to all types of search results which are represented by flat files
 * created by external programs e.g. Fasta, (T)BlastN/PX, EMBOSS
 * programs. This is not limited to sequence similarity searches, but
 * includes any format consisting of a header followed by hits, each
 * of which may, or may not, have subhits.
 *
 * @author Keith James
 * @since 1.1
 */
public interface SearchContentHandler
{
    /**
     * <code>getMoreSearches</code> returns the state of the
     * <code>SearchContentHandler</code> with respect to further
     * searches from its data source. Used for handling streams of
     * search results.
     *
     * @return a <code>boolean</code> value.
     */
    public boolean getMoreSearches();

    /**
     * <code>setMoreSearches</code> sets the state of the
     * <code>SearchContentHandler</code>'s expectation of receiving
     * more results. Used for handling streams of search results.
     *
     * @param value a <code>boolean</code> value.
     */
    public void setMoreSearches(boolean value);

    /**
     * The <code>startSearch</code> method indicates the start of
     * useful search information.
     */
    public void startSearch();

    /**
     * The <code>endSearch</code> method indicates the end of useful
     * search information.
     */
    public void endSearch();

    /**
     * The <code>startHeader</code> method indicates the start of a
     * formatted header. This usually contains information relevant to
     * the search as a whole.
     */
    public void startHeader();

    /**
     * The <code>endHeader</code> method indicates the end of a
     * formatted header.
     */
    public void endHeader();

    /**
     * The <code>startHit</code> method indicates the start of a
     * formatted hit. This could be a single line, or a block of
     * lines.
     */
    public void startHit();

    /**
     * The <code>endHit</code> method indicates the end of a formatted
     * hit.
     */
    public void endHit();

    /**
     * The <code>startSubHit</code> method indicates the start of a
     * formatted subhit. There may be zero or more of these per hit.
     */
    public void startSubHit();

    /**
     * The <code>endSubHit</code> method indicates the end of a
     * formatted subhit.
     */
    public void endSubHit();

    /**
     * The <code>addSearchProperty</code> method adds a key/value pair
     * containing some property of the overall search result.
     *
     * @param key an <code>Object</code>.
     * @param value an <code>Object</code>.
     */
    public void addSearchProperty(Object key, Object value);

    /**
     * The <code>addHitProperty</code> method adds a key/value pair
     * containing some property of a particular hit.
     *
     * @param key an <code>Object</code>.
     * @param value an <code>Object</code>.
     */
    public void addHitProperty(Object key, Object value);

    /**
     * The <code>addSubHitProperty</code> method adds a key/value pair
     * containing some property of a particular subhit.
     *
     * @param key an <code>Object</code>.
     * @param value an <code>Object</code>.
     */
    public void addSubHitProperty(Object key, Object value);

    /**
     * <code>setQueryID</code> identifies the query sequence by a
     * name, ID or URN.
     *
     * @param queryID a <code>String</code> which should be an unique
     * identifer for the sequence.
     */
    public void setQueryID(String queryID);

    /**
     * <code>setDatabaseID</code> identifies the database searched by
     * a name, ID or URN.
     *
     * @param databaseID a <code>String</code> which should be an unique
     * identifier for the database searched.
     */
    public void setDatabaseID(String databaseID);
}
