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
 * Can supply data.
 * @author Richard Holland
 * @since 3.0
 */
public interface ThingReader extends ThingParserMember {

    /**
     * Checks to see if the reader could read a record from its data source
     * if asked to.
     * @return {@code true} if it can.
     * @throws IOException if the read operation fails.
     */
    public boolean canReadNextThing() throws IOException;

    /**
     * Instructs the reader to read a record from its data source. Should only
     * be called if {@link #canReadOne()} returns {@code true}.
     * @throws IOException if the read operation fails.
     */
    public void readNextThing() throws IOException;
}
