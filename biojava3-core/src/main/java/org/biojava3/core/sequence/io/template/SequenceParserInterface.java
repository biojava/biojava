/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.io.template;

import java.io.DataInput;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public interface SequenceParserInterface {

    public String getSequence(DataInput dataInput,int sequenceLength) throws Exception;
}
