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
 * Receives data about things and writes them to some form of output.
 * @author Richard Holland
 * @since 3.0
 */
public interface ThingWriter extends ThingReceiver {
    
    /**
     * If the writer had to open any resources, close them now.
     * @throws IOException if it could not.
     */
    public void close() throws IOException;
}