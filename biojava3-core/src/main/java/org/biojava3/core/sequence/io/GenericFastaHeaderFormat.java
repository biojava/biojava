/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.io;

import org.biojava3.core.sequence.io.template.FastaHeaderFormatInterface;
import org.biojava3.core.sequence.template.AbstractSequence;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GenericFastaHeaderFormat implements FastaHeaderFormatInterface {

    public String getHeader(AbstractSequence sequence) {
        String header = "";

        if(sequence.getOriginalHeader() != null && sequence.getOriginalHeader().length() > 0){
            header = sequence.getOriginalHeader();
        }else{
            header = sequence.getAccession().getID();
        }

        return header;

    }

}
