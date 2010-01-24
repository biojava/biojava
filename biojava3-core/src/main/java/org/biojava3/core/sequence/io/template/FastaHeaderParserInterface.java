/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.io.template;

import org.biojava3.core.sequence.template.AbstractSequence;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public interface FastaHeaderParserInterface<S extends AbstractSequence> {

    public void parseHeader(String header,S sequence);
}
