package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava.nbio.core.sequence.io.template.SequenceHeaderParserInterface;
import org.biojava.nbio.core.util.InputStreamProvider;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Optional;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.Consumer;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Read from a FASTA file (or gzipped FASTA file) and create a Java stream of {@link ProteinSequence} objects
 * for use in a functional programming paradigm.
 *
 * @author Gary Murphy
 * @since 7.0.3
 */
public class FastaStreamer {

	private final Path path;
	private int batchSize = 1_000;
	private SequenceHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser;
	private SequenceCreatorInterface<AminoAcidCompound> sequenceCreator;
	private LinkedHashMap<String, ProteinSequence> chunk = new LinkedHashMap<>();
	private Iterator<Map.Entry<String, ProteinSequence>> iterator = Collections.emptyIterator();
	private boolean closed = false;

	/**
	 * The constructor is private.  Created via the <tt>from(...)</tt> static factory method
	 *
	 * @param path the path to the file containing the FASTA content (possibly GZipped)
	 */
	private FastaStreamer(final Path path) {
		this.path = path;
	}

	public static FastaStreamer from(final Path path) {
		return new FastaStreamer(path);
	}

	public static FastaStreamer from(File file) {
		return from(file.toPath());
	}

	public FastaStreamer withHeaderParser(SequenceHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser) {
		this.headerParser = headerParser;
		return this;
	}

	public FastaStreamer withSequenceCreator(SequenceCreatorInterface<AminoAcidCompound> sequenceCreator) {
		this.sequenceCreator = sequenceCreator;
		return this;
	}

	public FastaStreamer batchSize(int size) {
		this.batchSize = size;
		return this;
	}

	/**
	 * Create a stream of protein sequences from the contents of the path
	 * @return the stream
	 * @throws IOException if there is an error opening the file
	 */
	public Stream<ProteinSequence> stream() throws IOException {
		InputStreamProvider provider = new InputStreamProvider();
		InputStream input = provider.getInputStream(getPath().toFile());
		FastaReader<ProteinSequence, AminoAcidCompound> reader = new FastaReader<>(input, getHeaderParser(), getSequenceCreator());
		Spliterator<ProteinSequence> source = new Spliterators.AbstractSpliterator<>(Integer.MAX_VALUE, Spliterator.IMMUTABLE | Spliterator.NONNULL) {
			@Override
			public boolean tryAdvance(Consumer<? super ProteinSequence> action) {
				if (closed) {
					return false;
				}
				ProteinSequence protein = next(reader);
				if (null == protein) {
					return false;
				}
				action.accept(protein);
				return true;
			}

			/**
			 * Fetch the next header/protein tuple from the cache.  If the cache is empty, fetch another
			 * batch from the source file
			 *
			 * @param reader
			 * 		the input stream from which the FASTA content is read
			 * @return the protein sequence
			 */
			private ProteinSequence next(FastaReader<ProteinSequence, AminoAcidCompound> reader) {
				try {
					if (!iterator.hasNext()) {
						chunk = reader.process(getBatchSize());
						if (null == chunk) {
							closed = true;
							reader.close();
							return null;
						}
						iterator = chunk.entrySet().iterator();
					}
					if (iterator.hasNext()) {
						Map.Entry<String, ProteinSequence> entry = iterator.next();
						return createSequence(entry.getKey(), entry.getValue());
					}
					closed = true;
					reader.close();
				} catch (IOException exception) {
					throw new RuntimeException(String.format("I/O error reading the FASTA file from '%s'", getPath()));
				}
				return null;
			}
		}; // Spliterator
		return StreamSupport.stream(source, false);
	}

	/**
	 * Create the sequence with the information from the header.  This implementation return the sequence as-is, but
	 * this is an opportunity for the implementer to build specific information into the user collection space
	 * of the sequence
	 *
	 * @param header the original header
	 * @param sequence the protein sequence
	 * @return the sequence
	 */
	protected ProteinSequence createSequence(String header, ProteinSequence sequence) {
		return sequence;
	}

	protected Path getPath() {
		return path;
	}

	protected int getBatchSize() {
		return batchSize;
	}

	protected SequenceHeaderParserInterface<ProteinSequence, AminoAcidCompound> getHeaderParser() {
		return Optional.ofNullable(headerParser).orElse(new GenericFastaHeaderParser<>());
	}

	public SequenceCreatorInterface<AminoAcidCompound> getSequenceCreator() {
		return Optional.ofNullable(sequenceCreator).orElse(new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
	}
}
