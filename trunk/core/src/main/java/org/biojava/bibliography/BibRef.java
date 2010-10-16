// BibRef.java
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
 * This class is a core class of the bibliographic data model - it
 * represents a bibliographic reference, a citation.
 * It is a super-class for all specialized citation types, but it
 * can also be instantiated and represent an additional specialized
 * citation type.
 * </p>
 *
 * <p>
 * The <em>BibRef</em> class has several explicit attributes, which are
 * reasonably general and which originate from the
 * <a href="http://dublincore.org">Dublin Core Metadata</a>, and a
 * hashtable that can hold any number of additional attributes. The same pattern is
 * repeatedly used on several other places of the data model. It achieves
 * extendibility without losing interoperability if the following rules are obeyed:
 * <ul>
 *   <li> The implementation must be prepared for cases when the explicitly
 *        defined attributes are empty (containing <tt>null</tt> value, or,
 *        in case of arrays, an empty list of elements). 
 *   <li> The names of additional properties (keys of the hashtable) must be
 *        obtainable and/or checkable using controlled vocabularies
 *        (see interface {@link BibRefSupport} for details).
 *   <li> The values stored in that hashtable should be of "reasonable" types. Any
 *        implementation should understand at least basic Java types (in most cases
 *        the <tt>String</tt> type is the best choice). The more exotic types
 *        are used, the less interoperability between implementations is likely.
 * </ul>
 * </p>
 *
 * <p>
 * The <em>BibRef</em> class is a parent class for derived classes representing bibliographic
 * references to specialized bibliographic resources. The following classes are
 * defined explicitly:
 *   {@link BiblioBook},
 *   {@link BiblioArticle},
 *   {@link BiblioBookArticle},
 *   {@link BiblioJournalArticle},
 *   {@link BiblioPatent},
 *   {@link BiblioThesis},
 *   {@link BiblioProceeding},
 *   {@link BiblioTechReport}, and
 *   {@link BiblioWebResource}
 * </p>
 *
 * <p>
 * The active participants of the process of creation and dissemination of the
 * bibliographic resources are defined by the class {@link BiblioProvider} and its
 * sub-classes. The participants can be people, organizations, or even software
 * services (mainly used for new digital resources). The most obvious examples
 * are authors, but it includes also publishers and other contributors.
 * </p>
 *
 * <p>
 * And finally, there is a class {@link BiblioJournal} describing journals.
 * The citations referring to the journal articles have a reference to this class. 
 * </p>
 *
 * <p>
 * This is an overview of all participating classes and their attributes:
 * <img src="doc-files/bibobjects_java.jpg">.
 * </p>
 *
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @version $Id$
 * @since 1.3
 */

public class BibRef {

    /**
     * Additional attributes of this citation.
     */
    public Hashtable properties = new Hashtable();

    /** It is an unambiguous reference to this citation "within the world".
     *  It is a string conforming to an identification system. An example
     *  of such system can be a combination of a well-known repository name
     *  and a unique identifier defined within this repository, such as
     *  <tt>MEDLINE/20000003</tt>.
     */
    public String identifier;

    /**
     * <p>
     * It defines the nature or genre of the cited resource.
     * </p>
     *
     * <p>
     * A recommended best practice is to use only values from a controlled vocabulary
     * named as defined in {@link BibRefSupport#RESOURCE_TYPES}.
     * Syntactically, and because of making query navigation easier, the value of this
     * attribute should be equal to a constant predefined in {@link BibRefSupport},
     * such as {@link BibRefSupport#TYPE_BOOK} for books, or
     * {@link BibRefSupport#TYPE_JOURNAL_ARTICLE} for journal articles.
     * However, there may be bibliographic resources, which are
     * not defined by specialized sub-classes (for example,  letters, practical guideline,
     * or archives), and therefore they do not have predefined names in
     * {@link BibRefSupport} interface.
     * </p>
     *
     * <p>
     * <em>Note that for the description of the physical or digital manifestation of the
     * cited resource there is an attribute {@link #format}.</em>
     * </p>
     */
    public String type;

    /**
     * <p>
     * It is an array of identifiers, all of them pointing to <em>the same cited source</em>
     * but usually stored in different bibliographic repositories.
     * </p>
     *
     * <p>
     * <em>Note that this attribute is not for referencing citations to other documents
     * that are related to the cited document.</em>
     * </p>
     */
    public String[] crossReferences;

    /**
     * A title given to the cited resource (a name by which the resource is formally known).
     */
    public String title;

    /**
     * It defines the topic of the content of the cited resource.
     */
    public BiblioSubject subject;

    /**
     * An account of the content of the cited resource.
     * It is either an abstract, or table of contents, or both.
     * It can be written in a language different from the language of the cited resource.
     */
    public BiblioDescription description;

    /**
     * It defines an extent or scope of the content of the cited resource.
     * It can include spatial location (a place name or geographic co-ordinates),
     * temporal period (a period label, date, or date range), or both.
     */
    public BiblioScope coverage;

    /**
     * <p>
     * The authors and contributors are responsible for creating the contents of the cited resource.
     * There is no formal definition of how this responsibility is divided between them. However,
     * the authors are usually primary creators while contributors may be illustrators, translators,
     * or other creative providers.
     * </p>
     *
     * <p>
     * The authors are in an ordered array (to be able to find the first author).
     * </p>
     */
    public BiblioProvider[] authors;

    /**
     * <p>
     * The authors and contributors are responsible for creating the contents of the cited resource.
     * There is no formal definition of how this responsibility is divided between them. However,
     * the authors are usually primary creators while contributors may be illustrators, translators,
     * or other creative providers.
     * </p>
     *
     * <p>
     * The contributors are in an ordered array (to be able to find the first contributor).
     * </p>
     */
    public BiblioProvider[] contributors;

    /** 
     * A publisher is responsible for making the resource available.
     */
    public BiblioProvider publisher;

    /**
     * <p>
     * It specifies information about rights over the cited resource.
     * Typically, it contains a rights management statement for the resource, or it refers
     * to a service providing such information. Rights information often encompasses
     * <a href="http://www.itds.treas.gov/ITDS/ITTA/ipr.html">Intellectual Property Rights</a>,
     * Copyright, and various Property Rights.
     * </p>
     *
     * <p>
     * If the attribute is empty, no assumptions can be made about the status of these and
     * other rights with respect to the cited resource.
     * </p>
     */
    public String rights;

    /**
     * <p>
     * Defines a date associated with an event in the life cycle of the cited resource
     * when this resource became available. Usually, it is a date of publishing. However,
     * for not yet published resources, it can be a date of creation.
     * </p>
     *
     * <p>
     * The suggested encoding is as defined in a W3C NOTE  
     * <a href="http://www.w3.org/TR/NOTE-datetime">Date and Time Formats</a>.
     * This NOTE defines a profile of <a href="http://www.iso.ch/markete/8601.pdf">ISO8601 standard</a>.
     * ISO8601 describes a large number of date/time formats and the NOTE reduces the scope and restricts
     * the supported formats to a small number. The profile offers a number of options from which this
     * attribute should contain/permit only the following ones:
     * </p>
     *
     * <p>
     * <dl>
     *   <dt> Year
     *   <dd> YYYY (e.g., 2000)
     *   <dt> Year and month
     *   <dd> YYYY-MM (e.g., 2000-12)
     *   <dt> Complete date
     *   <dd> YYYY-MM-DD (e.g., 2000-12-31)
     *   <dt> Complete date plus hours, minutes, and seconds
     *   <dd> YYYY-MM-DDThh:mm:ssZ (e.g., 2000-12-31T23:59:59Z)
     * </dl>
     * </p>
     *
     * <p>
     * Exactly the components shown here must be present, with exactly this punctuation.
     * Note that the <b>T</b> appears literally in the string, to indicate the beginning
     * of the time element, as specified in ISO 8601.
     * </p>
     *
     * <p>
     * Times are expressed in UTC (Coordinated Universal Time), with a special UTC designator
     * (<em>Z</em>), again as specified in ISO 8601.
     * </p>
     *
     * <p>
     * For query purposes, the format with fewer details is considered
     * as having all possible values in place of missing details. Thus, YYYY-MM would mean
     * all dates and times in the given month.
     * </p>
     */
    public String date;

    /**
     * <p>
     * It defines a language of the intellectual contents of the cited resource.
     * The recommendation is to use values as defined by
     * <a href="http://www.ietf.org/rfc/rfc1766.txt">RFC1766</a> which includes a two-letter
     * <em>Language Code</em> (taken from the ISO639 standard), followed optionally by a two-letter
     * <em>Country Code</em> (taken from the ISO3166 standard).
     * </p>
     *
     * <p>
     * For example, <tt>en</tt>  for English,
     * <tt>fr</tt> for French, or <tt>en-uk</tt>  for English used in the United Kingdom.
     * </p>
     *
     * <p>
     * Another possibility is to use
     * <a href="http://lcweb.loc.gov/marc/languages">MARC List of Languages</a>.
     * </p>
     *
     * <p>
     * In any case, the name of the used controlled vocabulary should be equal to
     * {@link BibRefSupport#LANGUAGES}.
     * </p>
     */
    public String language;

    /**
     * It describes the physical or digital manifestation of the cited
     * resource.  It can have very different content depending on the
     * citation type. Therefore, it is highly recommended to use a
     * controlled vocabulary to fill this attribute. The name of such
     * vocabulary should be equal to the type of the reosurce type
     * followed by {@link BibRefSupport#ATTR_FORMAT}. For example:
     * {@link BibRefSupport#TYPE_BOOK}/{@link BibRefSupport#ATTR_FORMAT}.
     */
    public String format;

    /**
     * It defines information related to the citation itself rather than to the cited resource.
     */
    public BiblioEntryStatus entryStatus;

}
