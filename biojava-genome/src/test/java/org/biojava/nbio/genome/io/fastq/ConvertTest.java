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

import java.io.File;
import java.io.FileWriter;
import java.util.List;
import java.util.Map;

import org.junit.Test;
import static org.junit.Assert.*;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * Round trip conversion functional tests.
 */
public final class ConvertTest {

	@Test
	public void testConvert() throws Exception
	{
		Map<FastqVariant, FastqReader> readers = Maps.newHashMap();
		readers.put(FastqVariant.FASTQ_SANGER, new SangerFastqReader());
		readers.put(FastqVariant.FASTQ_SOLEXA, new SolexaFastqReader());
		readers.put(FastqVariant.FASTQ_ILLUMINA, new IlluminaFastqReader());

		Map<FastqVariant, FastqWriter> writers = Maps.newHashMap();
		writers.put(FastqVariant.FASTQ_SANGER, new SangerFastqWriter());
		writers.put(FastqVariant.FASTQ_SOLEXA, new SolexaFastqWriter());
		writers.put(FastqVariant.FASTQ_ILLUMINA, new IlluminaFastqWriter());

		Map<FastqVariant, String> inputFileNames = Maps.newHashMap();
		inputFileNames.put(FastqVariant.FASTQ_SANGER, "sanger_full_range_as_sanger.fastq");
		inputFileNames.put(FastqVariant.FASTQ_SOLEXA, "solexa_full_range_as_solexa.fastq");
		inputFileNames.put(FastqVariant.FASTQ_ILLUMINA, "illumina_full_range_as_illumina.fastq");

		Map<FastqVariantPair, String> expectedFileNames = Maps.newHashMap();
		expectedFileNames.put(new FastqVariantPair(FastqVariant.FASTQ_SANGER, FastqVariant.FASTQ_SANGER), "sanger_full_range_as_sanger.fastq");
		expectedFileNames.put(new FastqVariantPair(FastqVariant.FASTQ_SANGER, FastqVariant.FASTQ_SOLEXA), "sanger_full_range_as_solexa.fastq");
		expectedFileNames.put(new FastqVariantPair(FastqVariant.FASTQ_SANGER, FastqVariant.FASTQ_ILLUMINA), "sanger_full_range_as_illumina.fastq");
		expectedFileNames.put(new FastqVariantPair(FastqVariant.FASTQ_SOLEXA, FastqVariant.FASTQ_SANGER), "solexa_full_range_as_sanger.fastq");
		expectedFileNames.put(new FastqVariantPair(FastqVariant.FASTQ_SOLEXA, FastqVariant.FASTQ_SOLEXA), "solexa_full_range_as_solexa.fastq");
		expectedFileNames.put(new FastqVariantPair(FastqVariant.FASTQ_SOLEXA, FastqVariant.FASTQ_ILLUMINA), "solexa_full_range_as_illumina.fastq");
		expectedFileNames.put(new FastqVariantPair(FastqVariant.FASTQ_ILLUMINA, FastqVariant.FASTQ_SANGER), "illumina_full_range_as_sanger.fastq");
		expectedFileNames.put(new FastqVariantPair(FastqVariant.FASTQ_ILLUMINA, FastqVariant.FASTQ_SOLEXA), "illumina_full_range_as_solexa.fastq");
		expectedFileNames.put(new FastqVariantPair(FastqVariant.FASTQ_ILLUMINA, FastqVariant.FASTQ_ILLUMINA), "illumina_full_range_as_illumina.fastq");

		for (FastqVariant variant1 : FastqVariant.values())
		{
			FastqReader reader = readers.get(variant1);
			String inputFileName = inputFileNames.get(variant1);
			for (FastqVariant variant2 : FastqVariant.values())
			{
				FastqWriter writer = writers.get(variant2);
				String expectedFileName = expectedFileNames.get(new FastqVariantPair(variant1, variant2));

				File tmp = File.createTempFile("convertTest", "fastq");
				FileWriter fileWriter = new FileWriter(tmp);

				for (Fastq fastq : reader.read(getClass().getResource(inputFileName))) {
					writer.append(fileWriter, fastq);
				}

				fileWriter.close();

				FastqReader resultReader = readers.get(variant2);
				List<Fastq> observed = Lists.newArrayList(resultReader.read(tmp));
				List<Fastq> expected = Lists.newArrayList(resultReader.read(getClass().getResource(expectedFileName)));

				assertEquals(expected.size(), observed.size());
				for (int i = 0; i < expected.size(); i++)
				{
					assertEquals(expected.get(i).getDescription(), observed.get(i).getDescription());
					assertEquals(expected.get(i).getSequence(), observed.get(i).getSequence());
					assertEquals(expected.get(i).getQuality(), observed.get(i).getQuality());
					assertEquals(expected.get(i).getVariant(), observed.get(i).getVariant());
				}
			}
		}
	}

	private final class FastqVariantPair
	{
		final FastqVariant variant1;
		final FastqVariant variant2;

		FastqVariantPair(final FastqVariant variant1, final FastqVariant variant2)
		{
			this.variant1 = variant1;
			this.variant2 = variant2;
		}

		@Override
		public int hashCode()
		{
			int result = 47;
			result = 31 * result + variant1.hashCode();
			result = 31 * result + variant2.hashCode();
			return result;
		}

		@Override
		public boolean equals(final Object o)
		{
			if (o == this)
			{
				return true;
			}
			if (!(o instanceof FastqVariantPair))
			{
				return false;
			}
			FastqVariantPair pair = (FastqVariantPair) o;
			return variant1.equals(pair.variant1) && variant2.equals(pair.variant2);
		}
	}
}
