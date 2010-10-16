// BiblioDescription.java
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
 * It represents an account of the content of the cited resource.
 * It is either an abstract, or table of contents, or both.
 * It can be written in a language different from the language of the cited resource.
 * </p>
 *
 * <p>
 * Both abstract and table of contents can contain more than just a plain text,
 * typically they may be expressed in a <em>markup language</em>.
 * Their formats are defined according to the
 * <a href="http://www.ietf.org/rfc/rfc2045.txt">MIME</a> specification.
 * </p>
 *
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @version $Id$
 * @since 1.3
 */

public class BiblioDescription {

    /**
     * <p>
     * It is an abstract of the cited resource. It can be expressed as a plain text
     * or in a <em>markup language</em>.
     * </p>
     *
     * @see #abstractType
     */
    public String theAbstract;

    /**
     * <p>
     * It specifies how {@link #theAbstract} is coded.
     * </p>
     *
     * <p>
     * If it is empty then {@link #theAbstract} is coded as a plain text, using <tt>us-ascii</tt>
     * coding.
     * Otherwise, this attribute is equivalent to the  <em>Content-Type Header Field</em>
     * of the <a href="http://www.ietf.org/rfc/rfc2045.txt">MIME</a> specification,
     * with exclusion of the keyword <tt>Content-Type</tt>.
     * For example, it can contain <tt>text/html</tt>, or, using additional parameters,
     * <tt>text/plain; charset=us-ascii</tt>.
     * </p>
     *
     * <p>
     * Often abstracts are also available from the same or separate repository as URLs.
     * There are several ways to provide this information in the here described data model.
     * The implementations may choose their own way and still remain compliant with this
     * specification. However, the first approach, described below, is recommended to achieve
     * interoperability between implementations.
     * </p>
     *
     * <p>
     * <ul>
     *   <li> Use here <tt>text/url</tt>
     *        and put the URL into {@link #theAbstract} field.
     *   <li> Use here <tt>text/plain; url=xxxxx</tt>
     *        where xxxxx is a URL of the abstract
     *        (in this case {@link #theAbstract} may still have a full or partial text of the
     *        abstract as a plain text).
     *   <li> Use here a <em>multi-part</em>
     *        (see <a href="http://www.ietf.org/rfc/rfc2045.txt">MIME</a> specification).
     *        In such case {@link #theAbstract} will have both the full or partial abstract text,
     *        and a URL.
     *   <li> Put the URL into {@link BibRef#properties} using key <tt>abstractURL</tt>.
     * </ul>
     * </p>
     */
    public String abstractType;

    /**
     * <p>
     * It is a table of contents of the cited resource.
     * It can be expressed as a plain text or in a <em>markup language</em>.
     * </p>
     *
     * @see #tableOfContentsType
     */
    public String tableOfContents;

    /**
     * <p>
     * It specifies how {@link #tableOfContents} is coded.
     * </p>
     *
     * <p>
     * If it is empty then {@link #tableOfContents} is coded as a plain text, using <tt>us-ascii</tt>
     * coding.
     * Otherwise, this attribute is equivalent to the  <em>Content-Type Header Field</em>
     * of the <a href="http://www.ietf.org/rfc/rfc2045.txt">MIME</a> specification,
     * with exclusion of the keyword <tt>Content-Type</tt>.
     * </p>
     *
     * @see #abstractType abstractType for example
     */
    public String tableOfContentsType;

    /**
     * <p>
     * It defines a language used for {@link #theAbstract} and {@link #tableOfContents}.
     * The recommended values are discussed in {@link BibRef#language}.
     * </p>
     *
     * <p>
     * <em>Note that there is no mechanism how to specify different languages for
     * the abstract and table of contents for one citation.</em>
     * </p>
     */
    public String language;

}
