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
 * Created on 01-21-2010
 */
package org.biojava3.core.sequence.io;

import org.biojava3.core.sequence.io.template.GenbankHeaderParserInterface;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

public class GenericGenbankHeaderParser<S extends AbstractSequence<C>, C extends Compound> implements GenbankHeaderParserInterface<S,C> {

    /**
     * Parse out the components where some have a | and others do not
     * @param header
     * @return
     */
    private String[] getHeaderValues(String header) {
        String[] data = new String[0];
        return data;
    }

    /**
     * Parse the header and set the values in the sequence
     * @param header
     * @param sequence
     */
    public void parseHeader(String header, S sequence) {
        sequence.setOriginalHeader(header);
    }

    /**
     * 
     * @param args
     */
    public static void main(String[] args) {

    }
}
