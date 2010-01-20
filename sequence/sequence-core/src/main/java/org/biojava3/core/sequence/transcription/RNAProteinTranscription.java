/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.transcription;

import org.biojava3.core.sequence.RNASequence;

/**
 *
 * @author Scooter
 */
public interface RNAProteinTranscription {


    public String translate(RNASequence rnaCodingSequence) throws Exception;

}
