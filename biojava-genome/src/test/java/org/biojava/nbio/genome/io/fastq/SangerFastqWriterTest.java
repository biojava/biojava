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
package org.biojava.nbio.genome.io.fastq;


import org.junit.Test;

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

	@Test
	public void testConvertNotSangerVariant() throws Exception
	{
		SangerFastqWriter writer = new SangerFastqWriter();
		Appendable appendable = new StringBuilder();
		Fastq invalid = new FastqBuilder()
			.withDescription("description")
			.withSequence("sequence")
			.withQuality("quality_")
			.withVariant(FastqVariant.FASTQ_SOLEXA)
			.build();

		writer.append(appendable, invalid);
	}
}
