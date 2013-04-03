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

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.io.template.FastaHeaderFormatInterface;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

/**
 * We store the original header if the sequence is parsed from a fasta file and will use that exact
 * sequence if we write out the sequences to a fasta file. If we don't have an orginal header then
 * use the accession id. This allows the implementation by the user to write out complex header
 * with id notes etc without rewriting the fasta writer
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
