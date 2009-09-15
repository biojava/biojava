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
package org.biojava3.core.io;

import java.io.IOException;

/**
 * Receives data via the {@code set()} methods.
 * @author Richard Holland
 * @since 3.0
 */
public interface ThingReceiver {

    /**
     * Indicates that all subsequent things belong together. Will be called
     * even if only one thing is in the group. Remains in effect until
     * {@link #endGroupOfThings()} is called.
     * @throws IOException if IO problems occurred.
     */
    public void startGroupOfThings() throws IOException;

    /**
     * Indicates the end of a group of things. The builder must revert to the
     * state it was in before {@link #startGroupOfThings()} was called.
     * @exception IllegalStateException if called before 
     * {@link #finishThing()} has been called.
     * @throws IOException if IO problems occurred.
     */
    public void endGroupOfThings() throws IOException;

    /**
     * Inform the builder that it should start expecting data for a new Thing.
     * @exception IllegalStateException if called before 
     * {@link #startGroupOfThings()} has been called.
     * @throws IOException if IO problems occurred.
     */
    public void startThing() throws IOException;

    /**
     * Inform the builder that all data for the Thing it is currently building
     * has now been sent.
     * @exception IllegalStateException if called before 
     * {@link #startThing()} has been called.
     * @throws IOException if IO problems occurred.
     */
    public void finishThing() throws IOException;
    
    /**
     * If the receiver had to open any resources, close them now.
     * @throws IOException if it could not.
     */
    public void close() throws IOException;
}
