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

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import java.io.*;
import java.net.URL;
import java.nio.charset.Charset;
import java.util.List;

/**
 * Abstract reader implementation for FASTQ formatted sequences.
 *
 * @since 3.0.3
 */
abstract class AbstractFastqReader
	implements FastqReader
{

	/**
	 * Return the FASTQ sequence format variant for this reader.
	 *
	 * @return the FASTQ sequence format variant for this reader
	 */
	protected abstract FastqVariant getVariant();

	@Override
	public final void parse(final Readable readable, final ParseListener listener)
		throws IOException
	{
		FastqParser.parse(readable, listener);
	}

	@Override
	public final void stream(final Readable readable, final StreamListener listener)
		throws IOException
	{
		StreamingFastqParser.stream(readable, getVariant(), listener);
	}

	@Override
	public final Iterable<Fastq> read(final File file) throws IOException
	{
		if (file == null)
		{
			throw new IllegalArgumentException("file must not be null");
		}

		BufferedReader reader = null;
		Collect collect = new Collect();
		try
		{
			reader = new BufferedReader(new FileReader(file));
			stream(reader, collect);
		}
		finally
		{
			try
			{
				reader.close();
			}
			catch (Exception e)
			{
				// ignore
			}
		}
		return collect.getResult();
	}

	@Override
	public final Iterable<Fastq> read(final URL url) throws IOException
	{
		if (url == null)
		{
			throw new IllegalArgumentException("url must not be null");
		}

		BufferedReader reader = null;
		Collect collect = new Collect();
		try
		{
			reader = Resources.asCharSource(url, Charset.forName("UTF-8")).openBufferedStream();
			stream(reader, collect);
		}
		finally
		{
			try
			{
				reader.close();
			}
			catch (Exception e)
			{
				// ignore
			}
		}
		return collect.getResult();
	}

	@Override
	public final Iterable<Fastq> read(final InputStream inputStream) throws IOException
	{
		if (inputStream == null)
		{
			throw new IllegalArgumentException("inputStream must not be null");
		}

		BufferedReader reader = null;
		Collect collect = new Collect();
		try
		{
			reader = new BufferedReader(new InputStreamReader(inputStream));
			stream(reader, collect);
		}
		finally
		{
			try
			{
				reader.close();
			}
			catch (Exception e)
			{
				// ignore
			}
		}
		return collect.getResult();
	}

	/**
	 * Collect FASTQ formatted sequences in a list.
	 */
	private static final class Collect implements StreamListener
	{
		/** List of FASTQ formatted sequences. */
		private final List<Fastq> result = Lists.newLinkedList();

		@Override
		public void fastq(final Fastq fastq)
		{
			result.add(fastq);
		}

		/**
		 * Return an unmodifiable iterable over the FASTQ formatted sequences collected by this stream listener.
		 *
		 * @return an unmodifiable iterable over the FASTQ formatted sequences collected by this stream listener
		 */
		public Iterable<Fastq> getResult()
		{
			return Iterables.unmodifiableIterable(result);
		}
	}
}
