/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence;

import org.biojava3.core.sequence.template.AbstractSequence;

/**
 *
 * @author Scooter
 */
public class SequenceLengthError extends Error{

    public SequenceLengthError(String message){
        super(message);
    }
}
