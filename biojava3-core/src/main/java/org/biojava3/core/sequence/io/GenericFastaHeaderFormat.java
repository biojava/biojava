/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.io.template.FastaHeaderFormatInterface;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GenericFastaHeaderFormat<S extends AbstractSequence<?>, C extends Compound> implements FastaHeaderFormatInterface<S, C> {

    public String getHeader(S sequence) {
        String header = "";

        if (sequence.getOriginalHeader() != null && sequence.getOriginalHeader().length() > 0) {
            header = sequence.getOriginalHeader();
        } else {
            AccessionID accessionID = sequence.getAccession();
            if (accessionID != null) {
                header = accessionID.getID();
            }
        }

        return header;
    }
}
