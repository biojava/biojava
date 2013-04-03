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

package org.biojava3.core.sequence.io.template;

import java.io.DataInput;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public interface SequenceParserInterface {
/**
 *
 * @param dataInput
 * @param sequenceLength
 * @return
 * @throws Exception
 */
    public String getSequence(DataInput dataInput,int sequenceLength) throws Exception;
}
