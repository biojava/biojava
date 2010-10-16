// BibRefQuery.java
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

import java.io.InputStream;
import java.util.Enumeration;
import java.util.Hashtable;

/**
 * <p>
 * The interface <em>BibRefQuery</em> is a fundamental part of the Bibliographic Query
 * Service. It allows searching for and retrieving citations from a bibliographic
 * repository. The result of the query methods is again of type BibRefQuery which
 * allows further to refine the query. When the caller is satisfied with the query
 * results, the retrieval methods can be used to get either a list of citations (of type
 * {@link BibRef}), or an XML document representing citations. 
 * </p>
 *
 * <p>
 * Squeezing all query and retrieval methods into one interface allows to build very
 * flexible systems, both distributed (where the client and repository parts are
 * executed on different computers) and stand-alone (both parts are linked
 * together into one process).
 * </p>
 *
 * <p>
 *<table border=0 cellpadding=10>
 *<tr><td valign="top">
 * For example, this picture shows a client linked together
 * with a repository implementation. These two parts communicate
 * entirely via <em>BibRefQuery</em> interface. Each of them can be replaced
 * without changing the other one.
 *</td><th>
 * <img src="doc-files/BibRefQuery_simple.jpg">
 *</th></tr>
 *<tr><td valign="top">
 * In this example, a client uses <em>BibRefQuery</em> interface to communicate only
 * with a local implementation of a distributed architecture (a CORBA communication
 * protocol in this case). The repository implementation is similarly shielded by the
 * same interface from the communication protocol layer.
 *</td><th>
 * <img src="doc-files/BibRefQuery_corba.jpg">
 *</th></tr>
 *<tr><td valign="top">
 * The last picture shows yet another example of a distributed
 * architecture showing the parts which must be changed when a
 * different communication protocol is used (the SOAP-HTTP in this
 * case). Note that both the client and repository implementation
 * remained untouched.
 *</td><th>
 * <img src="doc-files/BibRefQuery_soap.jpg" align="left">
 *</th></tr>
 *</table>
 * </p>
 *
 * <h3>The implementation is advised to used the following constructor</h3>
 *
 * <p>
 *<pre>
 *    public <em>NameOfAnImplementation</em> (String[] args, Hashtable props) {...}
 *</pre>
 *    where both <tt>args</tt> and <tt>props</tt> contain implementation
 *    specific parameters and properties. However, some properties are
 *    more probable to be used - the suggested names for them are defined
 *    either in this interface or in the "sister" interface {@link BibRefSupport}.
 * </p>
 *
 * <p>
 * The use of this constructor makes easier to load dynamically different
 * implementations.
 * </p>
 *
 * <p>
 * The methods of the <em>BibRefQuery</em> interface can be divided into three groups.
 * The first group deals with connections to bibliographic repositories - here are
 * methods {@link #connect connect}, {@link #disconnect disconnect}, {@link #destroy destroy},
 * {@link #getCollectionId getCollectionId}, and {@link #isReady isReady}.
 * </p>
 *
 * <p>
 * The second and the most interesting group contains the query methods.
 * As mentioned above, these methods (mostly) return an another query collection
 * which is again query-able. Here belong methods {@link #find find}, {@link #findByAuthor findByAuthor},
 * {@link #findById findById}, {@link #query query}, {@link #getBibRefCount getBibRefCount}, and {@link #sort sort}.
 * </p>
 *
 * <p>
 * The last group has methods for retrieving citations from the resulting query collection.
 * The retrieval methods also allow to return citations not fully populated with all
 * available attribute data (for example, the long abstracts can be asked for only later).
 * Here belong methods {@link #getAllBibRefs getAllBibRefs}, {@link #getBibRefs getBibRefs},
 * {@link #getAllIDs getAllIDs}, {@link #getAllBibRefsAsXML getAllBibRefsAsXML},
 * {@link #getBibRefsAsXML getBibRefsAsXML}, and {@link #getBibRefAsXML getBibRefAsXML}. 
 * </p>
 *
 * <a name="how_attr_names">
 *<h3>Simple and Qualified Attribute Names</h3>
 * <p>
 * There are several places where method arguments represent attribute names:
 * </p>
 *
 * <p>
 *<ul>
 *   <li> In query methods, such as {@link #find find}, a list of attributes that should be
 *        searched.
 *   <li> The query results are citations represented as {@link BibRef} instances
 *        but not necessarily fully populated - they may contain only a subset of attributes,
 *        the <em>excluded</em> attribute lists used in several methods.
 *   <li> The results may be ordered by one or more attributes (method {@link #sort sort}).
 *</ul>
 * </p>
 *
 * <p>
 * Therefore, this interface defines several rules how to specify attribute names whenever
 * they have to be expressed as strings. The existence of these rules will make the
 * implementations interoperable. But, of course, they can be ignored if the interoperability
 * is not an issue.
 * </p>
 *
 * <p>
 * The following rules define how to create stringified names for individual attributes.
 *<ol>
 *  <li> The best recommended practice is to find attribute names from a controlled
 *       vocabulary - see details in {@link BibRefSupport} interface.
 *  <li> The stringified names of attributes of class {@link BibRef}  are equal to the
 *       member names of this class. For example, <tt>identifier</tt>, <tt>type</tt>,
 *       <tt>title</tt>, <tt>authors</tt>.
 *  <li> The stringified names of attributes of sub-classes derived from class {@link BibRef},
 *       and of attributes of other classes, are also equal to the member names but additionally
 *       they must be qualified by the resource type using two underscores ( __ ).
 *       For example,  <tt>Book__isbn</tt>,   <tt>JournalArticle__from_journal</tt>,
 *       <tt>Journal__name</tt>.
 *       <blockquote><em>
 *           The somewhat unusual  double underscore  is suggested here because in some
 *           query languages (where the stringified attribute names can be used as
 *           variables) is an underscore the only non-alphabetic character allowed for
 *           variables.
 *       </em></blockquote>
 *  <li> The qualification part of the stringified name (together with underscores) can
 *       be omitted if there is no ambiguity. For example, if an implementation does not
 *       use property name  <tt>isbn</tt>  anywhere else, the  <tt>Book__isbn</tt>  can
 *       be replaced by simple <tt>isbn</tt>.
 *       <blockquote><em>
 *       Be aware, however, that dropping the qualifier may compromise extendibility
 *       because a caller that expects a unique attribute name may break if another
 *       citation type is added with the same attribute name.
 *       </em></blockquote>
 *  <li> The stringified names of the attributes from {@link BibRef#properties dynamic properties}
 *       are equal to their property names, applying the rule about qualification as defined above.
 *       Thus, for example, an attribute  <tt>registry_number</tt>  hidden in member
 *       {@link BibRef#properties} will be stringified simply as <tt>registry_number</tt>,
 *       and an attribute  <tt>location</tt>  hidden in properties of a sub-class
 *       representing books will be stringified as  <tt>Book__location</tt>.
 *  <li> The stringified names of the attributes from {@link BibRef#properties dynamic properties}
 *       for instances without their own sub-class must be qualified (as described above) by
 *       the contents of their {@link BibRef#type}.  For example, a citation can be of type
 *       <tt>letter</tt>,  but there is no sub-class  <tt>Letter</tt>. Therefore, an attribute
 *       <tt>type</tt>  has value  <tt>letter</tt>.  This value is then used to create a qualified
 *       stringified name  <tt>letter__subject</tt>.
 *  <li> The stringified names should be considered case-insensitive. Thus,  <tt>book__location</tt>
 *       is the same as  <tt>Book__location</tt>,  and  <tt>journalarticle__issue</tt>  equals to
 *       <tt>JournalArticle__issue</tt>.
 *</ol>
 * </p>
 *
 *<a name="how_criteria">
 *<h3>Query Matching and Ordering Criteria</h3>
 * <p>
 * Several methods dealing with queries and sorting use a list of criteria.
 * The criteria define how the matching or ordering should be done.
 * </p>
 *
 * <p>
 * Each criterion is fully defined by an instance of {@link BiblioCriterion}.
 * Such definitions can be obtained from a controlled vocabulary - see
 * {@link BibRefSupport#getSupportedCriteria()}.
 * </p>
 *
 * <p>
 * Because each criterion is uniquely identifiable by its name, the querying and
 * sorting methods use only lists of criteria names, not lists of full criteria
 * definitions. 
 * </p>
 *
 *<a name="how_excluded">
 *<h3>Excluded and Only-included attributes</h3>
 * <p>
 * Several methods use parameter with <em>excluded</em> attributes, or a
 * parameter with <em>only-included</em> attributes.
 * There are two different meanings and uses of such attributes lists.
 * </p>
 *
 * <p>
 * The first meaning is used by the query methods. They return a new query collection.
 * From the practical and performance reasons it may be sometimes useful to define
 * <b>in advance</b> that the citations representing the resulting query collection
 * do not need to contain all attributes. The <em>excluded</em> list of attribute
 * names defines what attributes are not needed - typical use is to exclude
 * abstracts which may be quite long.
 * This, using the <em>excluded</em> list in the query method means that the
 * resulting query collection will never have all attributes fully filled with data
 * (unless, of course, the implementation ignores the <em>excluded</em> list).
 * </p>
 *
 * <p>
 * The second meaning is for the retrieval methods. They return citations from
 * a current query collection and can decide that only some attributes in the returned
 * citations are filled with data (such parameter list is always named <em>onlyAttrs</em>).
 * It may again mean that less data will be transferred
 * but it is a <b>post-act</b> decision because the query collection has already all
 * data and only does not return them now, but the next retrieval method (on the 
 * same collection) can retrieve them.
 * </p>
 *
 * <p>
 * The both uses may be applied in different scenarios, and their efficiency is
 * very dependent on the repository implementation. Sometimes the creation of a
 * query collection already includes heavy data manipulation - therefore, the
 * first usage may help with performance, But sometimes the resulting query
 * collection is more or less a virtual collection and the real data transfer
 * is applied only when the citations are being retrieved. In this case, the later
 * scenario may be more efficient.
 * </p>
 *
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @author Matthew Pocock
 * @version $Id$
 * @see BibRef
 * @see BibRefSupport
 * @see BiblioCriterion
 * @since 1.3
 */

public interface BibRefQuery {

    /**
     * <p>
     * A property name specifying a list of excluded attribute names
     * (the type of the property value should be <tt>String[]</tt>).
     * </p>
     *
     * <p>
     * The list is used to define attributes which are not returned in the
     * resulting citations (see discussion on
     * <a href="#how_excluded">excluded attributes</a>).
     * </p>
     *
     * @see #find find
     * @see #query query
     */
    static final String PROP_EXCLUDED_ATTRS = "excluded";

    /**
     * A property name specifying a list of searching and ordering criteria
     * names (type of the property value should be <tt>String[]</tt>). See
     * discussion on <a href="#how_criteria">criteria</a>.
     *
     * @see #find find
     * @see #query query
     */
    static final String PROP_CRITERIONS     = "criterions";

    /**************************************************************************
     * <p>
     * It returns an identification of the current query collection.
     * </p>
     *
     * <p>
     * At the beginning, the identification usually contain a bibliographic
     * repository name or its contents description. But later, usually after
     * {@link #connect} or after the first query, the identification may contain
     * information rich enough to be able to re-create the whole collection
     * (e.g. it can contain IDs of all records in the given collection).
     * </p>
     *
     * <p>
     * An implementation is not required to provide a persistent collection
     * identification. However, if it does provide, it should also be
     * able to accept the same identifier in the {@link #connect(byte[]) connect}
     * method, and to use it to re-create the same collection.
     * </p>
     *
     * @return an identification of the current collection (may be null)
     *************************************************************************/
    byte[] getCollectionId();

    /**************************************************************************
     * <p>
     * It creates a connection to a bibliographic repository, or/and it makes
     * all necessary initialization steps needed for further communication.
     * </p>
     *
     * <p>
     * However, there should be no need to call this method explicitly,
     * the other methods should do it automatically before they need something
     * from the repository.
     * </p>
     *
     * @throws BibRefException if the connection cannot be established
     *************************************************************************/
    void connect()
	throws BibRefException;

    /**************************************************************************
     * <p>
     * It creates a connection to a bibliographic repository, or/and it makes
     * all necessary initialization steps needed for further communication,
     * and it makes the collection described by <tt>collectionId</tt>
     * the current collection.
     * </p>
     *
     * @see #connect connect without parameters
     * @param collectionId a (usually persistent) token allowing to re-create
     *        a collection; the parameter is the same as an identifier returned
     *        earlier by method {@link #getCollectionId}
     * @throws BibRefException if the connection cannot be established, or if the
     *        collection with the given ID cannot be re-created
     *************************************************************************/
    void connect (byte[] collectionId)
	throws BibRefException;

    /**************************************************************************
     * It checks if the repository is available. The semantic of 
     * <em>available</em>depends on the implementation.
     *
     * @return true if it is ready
     *************************************************************************/
    boolean isReady();

    /**************************************************************************
     * <p>
     * It disconnects from the repository.
     * </p>
     *
     * <p>
     * The caller can use this method to announce that the current query
     * collection will not be needed soon. However, it may still be possible
     * to connect to it later again.
     * </p>
     *
     * @see #destroy destroy for more permanent action
     *************************************************************************/
    void disconnect();

    /*************************************************************************
     * <p>
     * It frees all resources related to this query collection.
     * </p>
     *
     * <p>
     * The caller explicitly announces no interest in the current
     * query collection at all. The existence of two separate
     * methods {@link #disconnect} and <tt>destroy</tt> allows more flexibility
     * for cases where an implementation deals with, for example,
     * temporary repositories. 
     * </p>
     *
     * @see #disconnect disconnect for less permanent action
     * @throws BibRefException if the connection to the repository is broken
     *************************************************************************/
    void destroy()
	throws BibRefException;

    /*************************************************************************
     * <p>
     * The easiest direct method for querying a repository.
     * </p>
     *
     * <p>
     * It is modeled on examples of web-based searches: A caller can specify
     * virtually anything in the list of keywords  and the implementation tries
     * to search for these in as many attributes as possible and reasonable,
     * applying logical AND between them. However, a caller can also specifically
     * limit the search only to attributes specified in the  searched  list.
     * </p>
     *
     * <p>
     * <em>Note that there is no real query language used by this method,
     * therefore, this method is not suitable for queries requiring
     * logical operators (others than AND).</em>
     * </p>
     *
     * <p>
     * The query result can be influenced by the additional properties:
     * <ul>
     *  <li> Property {@link #PROP_EXCLUDED_ATTRS} is of type <tt>String[]</tt>
     *       and contains list of attributes names which should not be
     *       included in the resulting query collection. See discussions on
     *       <a href="#how_excluded">excluded attributes</a> and on
     *       <a href="#how_attr_names">stringified attribute names</a>,
     *  <li> Property {@link #PROP_CRITERIONS} is also of type <tt>String[]</tt>
     *       and contains list of criteria names. The caller specifies here
     *       what criteria she wishes, and this method can change this property
     *       and return here the criteria really used for the query.
     *       See also discussion about <a href="#how_criteria">criteria</a>.
     * </ul>
     * </p>
     *
     * @param keywords keyword or phrases that are being looked for
     * @param attrs attributes names that should be searched; if this list is
     *              empty the implementation should search all reasonable
     *              attributes
     * @param properties specify attributes excluded from the results and
     *                   requested criteria for the query
     * @return a new query (and query-able) collection
     * @throws BibRefException if query failed (which can have many reasons :-))
     *         (note that an empty result does not cause an exception)
     *************************************************************************/
    BibRefQuery find (String[] keywords, String[] attrs, Hashtable properties)
	throws BibRefException;

    /*************************************************************************
     * <p>
     * This is a convenient method for a common query.
     * </p>
     *
     * <p>
     * The search is done only for attributes having non empty values in
     * parameter  <tt>author</tt>.  For example, a search for citations written
     * by authors with surname  <tt>Doe</tt>  can be specified by sending an
     * instance of <tt>BiblioPerson</tt> with <tt>surname</tt> filled with
     * <tt>Doe</tt>  and with other attributes empty. Or, a search for
     * institution  <tt>EBI</tt>  can be specified by sending an instance of
     * <tt>BiblioOrganization</tt> with <tt>name</tt> containing  <tt>EBI</tt>. 
     * </p>
     *
     * <p>
     * The query result can be influenced by the additional properties:
     * <ul>
     *  <li> Property {@link #PROP_EXCLUDED_ATTRS} is of type <tt>String[]</tt>
     *       and contains list of attributes names which should not be
     *       included in the resulting query collection. See discussions on
     *       <a href="#how_excluded">excluded attributes</a> and on
     *       <a href="#how_attr_names">stringified attribute names</a>,
     *  <li> Property {@link #PROP_CRITERIONS} is also of type <tt>String[]</tt>
     *       and contains list of criteria names. The caller specifies here
     *       what criteria she wishes, and this method can change this property
     *       and return here the criteria really used for the query.
     *       See also discussion about <a href="#how_criteria">criteria</a>.
     * </ul>
     * </p>
     *
     * @see #find find
     * @see BiblioPerson
     * @see BiblioOrganisation
     * @see BiblioService
     * @param author contains one or more attributes that are being search for
     * @param properties specify attributes excluded from the results and
     *                   requested criteria for the query
     * @return a new query (and query-able) collection
     * @throws BibRefException if query failed (which can have many reasons :-))
     *         (note that an empty result does not cause an exception)
     *************************************************************************/
    BibRefQuery findByAuthor (BiblioProvider author, Hashtable properties)
	throws BibRefException;

    /*************************************************************************
     * <p>
     * This is a convenient method returning just one citation.
     * </p>
     *
     * <p>
     * It queries the current collection in order to find and to retrieve
     * a citation with the given identifier. It depends on the implementation
     * what could be used as an identifier - see {@link BibRef#identifier}.
     * </p>
     *
     * @see #findById(String,String[]) findById with limited returned attributes
     * @param bibRefId an identifier of a citation that is being looked for
     * @return a found bibliographic reference (citation)
     * @throws BibRefException if such citation was not found (or something else
     *                     bad happened)
     *************************************************************************/
    BibRef findById (String bibRefId)
	throws BibRefException;

    /*************************************************************************
     * <p>
     * This is a convenient method returning just one citation, perhaps with
     * a limited number of attributes.
     * </p>
     *
     * <p>
     * It queries the current collection in order to find and to retrieve
     * a citation with the given identifier. It depends on the implementation
     * what could be used as an identifier - see {@link BibRef#identifier}.
     * </p>
     *
     * <p>
     * The returned citation will contain at least attributes whose names are
     * specified by the parameter <tt>onlyAttrs</tt> (see discussion on
     * <a href="#how_excluded">only-included attributes</a>.
     * </p>
     *
     * <p>
     * It is meant to provide more lightweight citation. The
     * implementation may provide more attributes than specified in
     * <tt>onlyAttrs</tt> (e.g. it is always recommended to include an
     * attribute representing a unique identifier of the citation even
     * if it is not asked for).
     * </p>
     *
     * <p>
     * Note that one can ask only for attributes that are available in the
     * current collection. If the collection was already created
     * <em>without</em> some attributes (using property
     * {@link #PROP_EXCLUDED_ATTRS}, e.g in method {@link #find find}) one cannot
     * expect to get them even if they are asked for by the parameter
     * <tt>onlyAttrs</tt>.
     * </p>
     *
     * @see #findById(String) findById
     * @param bibRefId an identifier of a citation that is being looked for
     * @param onlyAttrs a list of attribute names; at least these attributes
     *                  will be included in the returned citation
     * @return a found bibliographic reference (citation)
     * @throws BibRefException if such citation was not found (or something else
     *                     bad happened)
     *************************************************************************/
    BibRef findById (String bibRefId, String[] onlyAttrs)
	throws BibRefException;

    /*************************************************************************
     * <p>
     * It queries the current collection using a query language.
     * </p>
     *
     * <p>
     * Use this method when the simple {@link #find find} method is not sufficient.
     * For example, when more logical or relational operators are needed
     * to express the query,
     * </p>
     *
     * <p>
     * This specification does not propose any specific query language
     * to use (but may in the future). Roughly speaking, the query
     * method takes a <tt>query</tt> string and passes it to the repository
     * implementation, and if the implementation understands the query
     * the world is saved.
     * </p>
     *
     * <p>
     * Again, the query result can be influenced by the additional properties:
     * <ul>
     *  <li> Property {@link #PROP_EXCLUDED_ATTRS} is of type <tt>String[]</tt>
     *       and contains list of attributes names which should not be
     *       included in the resulting query collection. See discussions on
     *       <a href="#how_excluded">excluded attributes</a> and on
     *       <a href="#how_attr_names">stringified attribute names</a>,
     *  <li> Property {@link #PROP_CRITERIONS} is also of type <tt>String[]</tt>
     *       and contains list of criteria names. The caller specifies here
     *       what criteria she wishes, and this method can change this property
     *       and return here the criteria really used for the query.
     *       See also discussion about <a href="#how_criteria">criteria</a>.
     * </ul>
     * </p>
     *
     * @see #find find
     * @param query an expression in a query language
     * @param properties specify attributes excluded from the results and
     *                   requested criteria for the query
     * @return a new query (and query-able) collection
     * @throws BibRefException if query failed (which can have many reasons :-))
     *         (note that an empty result does not cause an exception)
     *************************************************************************/
    BibRefQuery query (String query, Hashtable properties)
	throws BibRefException;

    /*************************************************************************
     * <p>
     * It returns the number of citations in the current collection.
     * </p>
     *
     * @return the size of this collection
     * @throws BibRefException if a connection with the repository is broken
     *************************************************************************/
    int getBibRefCount()
	throws BibRefException;

    /*************************************************************************
     * <p>
     * It sorts the current collection and returns another collection which is
     * a sorted copy of the current collection.
     * </p>
     *
     * <p>
     * This is not strictly speaking a query method but it also returns
     * a query collection.
     * </p>
     *
     * <p>
     * The sorting result can be influenced by an additional property
     * {@link #PROP_CRITERIONS} (of type <tt>String[]</tt>) containing
     * a list of sorting criteria names. The caller specifies here
     * what criteria she wishes, and this method can change this property
     * and return here the criteria really used for sorting.
     * </p>
     *
     * @param orderedBy a list of attribute names that the collection should
     *                  be sorted by
     * @param properties FIXME: docs & params out of sync
     * @return a sorted collection
     * @throws BibRefException if sorting failed (which may also happen when
     *         the collection is too large)
     *************************************************************************/
    BibRefQuery sort (String[] orderedBy, Hashtable properties)
	throws BibRefException;

    /*************************************************************************
     * <p>
     * It returns all citations from the current collection as a
     * (possibly big) array. Obviously, the repository implementation
     * may limit the number of returned records.
     * </p>
     *
     * <p>
     * Some attributes may be missing (empty) if the property
     * {@link #PROP_EXCLUDED_ATTRS} was used for creating the current
     * collection. See discussion on
     * <a href="#how_excluded">excluded attributes</a>.
     * </p>
     *
     * @see #getAllBibRefs(String[]) getAllBibRefs with limited returned attributes
     * @return all citations from the current collection
     * @throws BibRefException if the collection is too large, or if the connection
     *        to the repository is broken
     *************************************************************************/
    BibRef[] getAllBibRefs()
	throws BibRefException;

    /*************************************************************************
     * <p>
     * It returns all citations from the current collection as a
     * (possibly big) array, perhaps with a limited number of attributes.
     * </p>
     *
     * <p>
     * The returned citations will contain at least attributes whose names are
     * specified by the parameter <tt>onlyAttrs</tt>. It is meant to provide
     * more lightweight citations. The implementation may provide more
     * attributes than specified in <tt>onlyAttrs</tt> (e.g. it may be always
     * good to include an attribute representing a unique identifier of a
     * citation even if it is not asked for). See discussion on
     * <a href="#how_excluded">only-included attributes</a>.
     * </p>
     *
     * <p>
     * Note that one can ask only for attributes that are available in the
     * current collection. If the collection was already created
     * <em>without</em> some attributes (using property
     * {@link #PROP_EXCLUDED_ATTRS}, e.g in method {@link #find find}) one
     * cannot expect to get them even if they are asked for by the parameter
     * <tt>onlyAttrs</tt>.
     * </p>
     *
     * @see #getAllBibRefs getAllBibRefs with all attributes
     * @see #getAllIDs getAllIDs
     *
     * @param onlyAttrs  attributes to attempt to include
     * @return all citations from the current collection
     * @throws BibRefException if the collection is too large, or if the connection
     *        to the repository is broken
     *************************************************************************/
    BibRef[] getAllBibRefs (String[] onlyAttrs)
	throws BibRefException;

    /*************************************************************************
     * <p>
     * A convenient method returning just identifiers of all current citations.
     * </p>
     *
     * @return a list of all identifiers
     * @throws BibRefException if the collection is too large, or if the connection
     *        to the repository is broken
     *************************************************************************/
    String[] getAllIDs()
	throws BibRefException;

    /*************************************************************************
     * <p>
     * It returns an enumeration of all citations from the current collection.
     * The type of elements in the enumeration is {@link BibRef} (or of its
     * sub-classes).
     * </p>
     *
     * <p>
     * Some attributes may be missing (empty) if the property
     * {@link #PROP_EXCLUDED_ATTRS} was used for creating the current
     * collection.
     * </p>
     *
     * @see #getAllBibRefs getAllBibRefs
     * @return an iterator over all citations
     * @throws BibRefException if the connection to the repository is broken
     *************************************************************************/
    Enumeration getBibRefs()
	throws BibRefException;

    /*************************************************************************
     * <p>
     * It returns an enumeration of all citations from the current collection,
     * perhaps with a limited number of attributes.
     * The type of elements in the enumeration is {@link BibRef} (or of its
     * sub-classes).
     * </p>
     *
     * <p>
     * The citations available through the enumeration will contain at least
     * attributes whose names are specified by the parameter <tt>onlyAttrs</tt>.
     * It is meant to provide more lightweight citations. The implementation
     * may provide more attributes than specified in <tt>onlyAttrs</tt> (e.g.
     * it may be always good to include an attribute representing a unique
     * identifier of a citation even if it is not asked for).
     * </p>
     *
     * <p>
     * Note that one can ask only for attributes that are available in the
     * current collection. If the collection was already created
     * <em>without</em> some attributes (using property
     * {@link #PROP_EXCLUDED_ATTRS}, e.g in method {@link #find find}) one cannot
     * expect to get them even if they are asked for by the parameter
     * <tt>onlyAttrs</tt>.
     * </p>
     *
     * @see #getAllBibRefs getAllBibRefs
     * @see #getBibRefs getBibRefs with all attributes
     *
     * @param onlyAttrs attributes to attempt to fetch
     * @return an iterator over all citations
     * @throws BibRefException if the connection to the repository is broken
     *************************************************************************/
    Enumeration getBibRefs (String[] onlyAttrs)
	throws BibRefException;

    /*************************************************************************
     * <p>
     * It returns all citations from the current collection as an XML stream.
     * The contents of such XML stream is widely repository dependent.
     * </p>
     *
     * <p>
     * Some attributes may be missing (empty) if the property
     * {@link #PROP_EXCLUDED_ATTRS} was used for creating the current
     * collection.
     * </p>
     *
     * @see #getAllBibRefs getAllBibRefs
     * @return an XML data stream containing all citations from the current
     *         collection
     * @throws BibRefException if the collection is too large, or if the connection
     *        to the repository is broken
     *************************************************************************/
    InputStream getAllBibRefsAsXML()
	throws BibRefException;

    /*************************************************************************
     * <p>
     * It returns an enumeration of all citations from the current collection.
     * The type of elements in the enumeration is <tt>String</tt>.
     * Each element represents one citation as an XML string.
     * The contents of such XML string is widely repository dependent.
     * </p>
     *
     * <p>
     * Some attributes may be missing (empty) if the property
     * {@link #PROP_EXCLUDED_ATTRS} was used for creating the current
     * collection.
     * </p>
     *
     * @see #getBibRefs getBibRefs
     * @see #getAllBibRefsAsXML getAllBibRefsAsXML
     * @return an iterator over all citations
     * @throws BibRefException if the connection to the repository is broken
     *************************************************************************/
    Enumeration getBibRefsAsXML()
	throws BibRefException;

    /*************************************************************************
     * <p>
     * A convenient utility method converting a given citation to its
     * XML representation. It is useful, for example, in cases when a
     * program annotates citations on-the-fly and needs them in the
     * same XML format.
     * </p>
     *
     * <p>
     * The XML format depends on the repository where the citation comes from.
     * </p>
     *
     * @param bibRef a citation being converted into an XML format
     * @return an XML representation of <tt>bibRef</tt>
     * @throws BibRefException if the implementation needs it :-)
     *************************************************************************/
    String getBibRefAsXML (BibRef bibRef)
	throws BibRefException;


}
