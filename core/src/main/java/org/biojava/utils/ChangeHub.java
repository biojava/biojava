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


package org.biojava.utils;


/**
 * Interface implemented by ChangeHubs, i.e.
 * classes that handle behaviour for
 * multiple instances of Changeable classes.
 * <p>
 * Listeners are indexed with a key and when
 * an event is fired, only listeners with the same
 * key are invoked.  The class manages the mapping
 * between key and listener.  It is the users responsibility
 * to compute the key.
 *
 * @author Thomas Down (original implementation)
 * @author David Huen (refactoring)
 * @since 1.3
 */
public interface ChangeHub
{
    /**
     * add a ChangeListener associated with given key.
     */
    public void addListener(Object key, ChangeListener listener, ChangeType ct);

    /**
     * remove a ChangeListener associated with given key.
     */
    public void removeListener(Object key, ChangeListener listener, ChangeType ct);

    /**
     * invoke the firePreChangeEvent on all ChangeListeners associated with
     * a specific key.
     */
    public void firePreChange(Object key, ChangeEvent cev) throws ChangeVetoException;

    /**
     * invoke the firePostChangeEvent on all ChangeListeners associated with
     * a specific key.
     */
    public void firePostChange(Object key, ChangeEvent cev);
}

