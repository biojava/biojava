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
 * Unit test for SangerFastqWriter.
 */
public final class SangerFastqWriterTest
    extends AbstractFastqWriterTest
{

    @Override
    public FastqWriter createFastqWriter()
    {
        return new SangerFastqWriter();
    }

    @Override
    public Fastq createFastq()
    {
        return new FastqBuilder()
            .withDescription("description")
            .withSequence("sequence")
            .withQuality("quality_")
            .withVariant(FastqVariant.FASTQ_SANGER)
            .build();
    }

    public void testValidateNotSangerVariant()
    {
        SangerFastqWriter writer = new SangerFastqWriter();
        Appendable appendable = new StringBuilder();
        Fastq invalid = new FastqBuilder()
            .withDescription("description")
            .withSequence("sequence")
            .withQuality("quality_")
            .withVariant(FastqVariant.FASTQ_SOLEXA)
            .build();
        try
        {
            writer.append(appendable, invalid);
            fail("validate not fastq-sanger variant expected IOException");
        }
        catch (IOException e)
        {
            // expected
        }
    }
}