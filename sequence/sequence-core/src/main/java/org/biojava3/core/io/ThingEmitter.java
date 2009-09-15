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

import java.io.Serializable;
import java.util.Iterator;

/**
 * Takes a thing and emits its data as if it were being parsed from a data
 * source.
 * @author Richard Holland
 * @since 3.0
 */
public interface ThingEmitter<T extends Serializable> extends ThingReader {

    /**
     * Sets the data source for {@link #canReadNextThing()} and friends to
     * read things from.
     * @param things an iterator of things to emit.
     */
    public void setThingSource(Iterator<T> things);
}
