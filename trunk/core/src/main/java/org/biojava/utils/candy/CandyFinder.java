// CandyFinder.java
//
//    senger@ebi.ac.uk
//    February 2001
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
package org.biojava.utils.candy;

import java.beans.PropertyChangeListener;

/**
 * <p>
 * This interface is a main entry point to a set of controlled
 * vocabularies.
 * </p>
 *
 * <p>
 * The implementation is supposed to behave as a Java bean (usually an
 * invisible bean unless it implements some additional GUI methods
 * which are not defined in this interface).
 * </p>
 *
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @author Matthew Pocock
 * @version $Id$
 */

public interface CandyFinder
    extends PropertyChangeListener {

    /**************************************************************************
     * <p>
     * It creates a connection to an object representing a vocabulary
     * finder, or/and it makes all necessary initialization steps
     * needed for further communication.
     * </p>
     *
     * <p>
     * However, there should be no need to call this method
     * explicitly, the other methods should do it automatically before
     * they need to use the finder.
     * </p>
     *
     * @throws CandyException if the connection/initialization cannot
     *         be established
     *************************************************************************/
    void connect() throws CandyException;

    /**************************************************************************
     * It checks if a vocabulary finder object is available. The semantic of 
     * <em>available</em>depends on the implementation.
     *************************************************************************/
    boolean isReady();

    /**************************************************************************
     * It closes connection with the finder object. Implementations may
     * choose to use this method for freeing resources.
     *************************************************************************/
    void disconnect();

    /*************************************************************************
     * <p>
     * It returns names of all vocabularies known to this vocabulary
     * finder. Any of the returned names can be later used in the method
     * {@link #getVocabularyByName getVocabularyByName}.
     * </p>
     *
     * @return a list of available vocabulary names
     * @throws CandyException if the finder fails to communicate
     *         with its vocabularies
     *************************************************************************/
    String[] getAllVocabularyNames() throws CandyException;

    /*************************************************************************
     * It returns a selected vocabulary.
     *
     * @see #getAllVocabularyNames
     * @param name a name of a vocabulary to be returned
     * @return a selected vocabulary
     * @throws CandyException when the vocabulary cannot be found (likely the
     *         given name is wrong)
     *************************************************************************/
    CandyVocabulary getVocabularyByName (String name)
      throws CandyException;

    /*************************************************************************
     * It returns all available vocabularies.
     *
     * @return all available vocabularies
     * @throws CandyException if the finder fails to communicate
     *         with its vocabularies
     *************************************************************************/
    CandyVocabulary[] getAllVocabularies() throws CandyException;

    /*************************************************************************
     * It returns the number of available vocabularies.
     *
     * @throws CandyException if the finder fails to communicate
     *         with its vocabularies
     *************************************************************************/
    int getNumCount() throws CandyException;

    /*************************************************************************
     * 
     *    P r o p e r t i e s
     *
     *************************************************************************/

    /** A default name of this (and any) finder.
     *  It is used when no other name was given by the finder
     *  implementation.
     */
    static final String DEFAULT_FINDER_NAME = "Default vocabulary finder";

    //
    // Property names
    //

    /**
     * <p>
     * A property name.
     * </p>
     *
     * <p>
     * Its value is of type {@link CandyVocabulary}.
     * It this property is set a given vocabulary becomes part of
     * this finder.
     * </p>
     */
    static final String PROP_VOCABULARY = CandyVocabulary.PROP_VOCABULARY;

    //
    // Property access methods
    //

    /*************************************************************************
     * It returns a name of this vocabulary finder.
     *
     * @throws CandyException if the finder fails to return its name
     *************************************************************************/
    String getFinderName() throws CandyException;

}
