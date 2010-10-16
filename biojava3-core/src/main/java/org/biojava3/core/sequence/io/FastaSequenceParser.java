/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import java.io.DataInput;
import org.biojava3.core.sequence.io.template.SequenceParserInterface;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FastaSequenceParser implements SequenceParserInterface {

    public String getSequence(DataInput dataInput, int sequenceLength) throws Exception {
        StringBuilder sb;
        if (sequenceLength != -1) {
            sb = new StringBuilder(sequenceLength);
        } else {
            sb = new StringBuilder();
        }
        boolean keepGoing = true;
        while (keepGoing) {
            String line = dataInput.readLine();
            if (line == null || line.startsWith(">")) {
                break;
            }
            sb.append(line.trim());
        }
        return sb.toString();
    }
}
