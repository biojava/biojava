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
 */
package org.biojava3.sequencing.io.fastq;

import java.io.IOException;

/**
 * Writer for {@link FastqVariant#FASTQ_SANGER} formatted sequences.
 *
 * @since 3.0.3
 */
public final class SangerFastqWriter
    extends AbstractFastqWriter
{

    @Override
    protected void validate(final Fastq fastq) throws IOException
    {
        if (fastq == null)
        {
            return;
        }
        if (!fastq.getVariant().isSanger())
        {
            throw new IOException("sequence " + fastq.getDescription()
                                  + " not fastq-sanger format, was " + fastq.getVariant().lowercaseName());
        }
    }
}