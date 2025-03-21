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
 * @author Scooter Willis ;lt;willishf at gmail dot com&gt;
 * @author Karl Nicholas &lt;github:karlnicholas&gt;
 * @author Paolo Pavan
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 01-21-2010
 */
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.DataSource;
import org.biojava.nbio.core.sequence.TaxonomyID;
import org.biojava.nbio.core.sequence.features.DBReferenceInfo;
import org.biojava.nbio.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava.nbio.core.sequence.io.template.SequenceHeaderParserInterface;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.Compound;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Use {@link GenbankReaderHelper} as an example of how to use this class where {@link GenbankReaderHelper} should be the
 * primary class used to read Genbank files
 * @param <S> the sequence type
 * @param <C> the compound type
 */
public class GenbankReader<S extends AbstractSequence<C>, C extends Compound> {

	private SequenceCreatorInterface<C> sequenceCreator;
	private GenbankSequenceParser<S,C> genbankParser;
	private BufferedReader bufferedReader;
	private boolean closed;
	private final Logger logger = LoggerFactory.getLogger(this.getClass());

	public boolean isClosed() {
		return closed;
	}

	/**
	 * If you are going to use {@link FileProxyProteinSequenceCreator} then do not use this constructor because we need details about
	 * local file offsets for quick reads. {@link InputStream} does not give you the name of the stream to access quickly via file seek. A seek in
	 * an {@link InputStream} is forced to read all the data so you don't gain anything.
	 * @param is
	 * @param headerParser
	 * @param sequenceCreator
	 */
	public GenbankReader(final InputStream is, final SequenceHeaderParserInterface<S,C> headerParser,
						 final SequenceCreatorInterface<C> sequenceCreator) {
		this.sequenceCreator = sequenceCreator;
		bufferedReader = new BufferedReader(new InputStreamReader(is));
		genbankParser = new GenbankSequenceParser<>();
		closed = false;
	}

	/**
	 * If you are going to use the FileProxyProteinSequenceCreator then you
	 * need to use this constructor because we need details about
	 * the location of the file.
	 * @param file
	 * @param headerParser
	 * @param sequenceCreator
	 * @throws FileNotFoundException if the file does not exist, is a directory
	 * 	rather than a regular file, or for some other reason cannot be opened
	 * 	for reading.
	 * @throws SecurityException if a security manager exists and its checkRead
	 * 	method denies read access to the file.
	 */
	public GenbankReader(
			final File file,
			final SequenceHeaderParserInterface<S,C> headerParser,
			final SequenceCreatorInterface<C> sequenceCreator
			) throws FileNotFoundException {

		this.bufferedReader = new BufferedReader(new FileReader(file));
		this.sequenceCreator = sequenceCreator;
		genbankParser = new GenbankSequenceParser<>();
	}

	/**
	 * The parsing is done in this method.<br>
	 * This method will return all the available Genbank records
	 * in the File or InputStream, closes the underlying resource,
	 * and return the results in {@link LinkedHashMap}.<br>
	 * You don't need to call {@link GenbankReader#close()} after calling this method.
	 * @see #process(int)
	 * @return {@link HashMap} containing all the parsed Genbank records
	 * present, starting current fileIndex onwards.
	 * @throws IOException
	 * @throws CompoundNotFoundException
	 * @throws OutOfMemoryError if the input resource is larger than the allocated heap.
	 */
	public Map<String, S> process() throws IOException, CompoundNotFoundException {
		Map<String, S> result = process(-1);
		close();
		return result;
	}

	/**
	 * This method tries to parse maximum <code>max</code> records from
	 * the open File or InputStream, and leaves the underlying resource open.<br>
	 *
	 * Subsequent calls to the same method continue parsing the rest of the file.<br>
	 * This is particularly useful when dealing with very big data files,
	 * (e.g. NCBI nr database), which can't fit into memory and will take long
	 * time before the first result is available.<br>
	 * <b>N.B.</b>
	 * <ul>
	 * <li>This method can't be called after calling its NO-ARGUMENT twin.</li>
	 * <li>remember to close the underlying resource when you are done.</li>
	 * </ul>
	 * @see #process()
	 * @author Amr ALHOSSARY
	 * @since 3.0.6
	 * @param max maximum number of records to return.
	 * @return {@link HashMap} containing maximum <code>max</code> parsed Genbank records
	 * present, starting current fileIndex onwards.
	 * @throws IOException
	 * @throws CompoundNotFoundException
	 */
	public Map<String, S> process(final int max) throws IOException, CompoundNotFoundException {

		if(closed){
			throw new IOException("Cannot perform action: resource has been closed.");
		}

		Map<String, S> sequences = new LinkedHashMap<>();
		int i=0;
		while(true) {
			if(max>0 && i>=max) break;
			i++;
			String seqString = genbankParser.getSequence(bufferedReader, 0);
			//reached end of file?
			if(seqString==null) break;
			@SuppressWarnings("unchecked")
			S sequence = (S) sequenceCreator.getSequence(seqString, 0);
			GenericGenbankHeaderParser<S, C> genbankHeaderParser = genbankParser.getSequenceHeaderParser();			
			genbankHeaderParser.parseHeader(genbankParser.getHeader(), sequence);			
			String id = genbankHeaderParser.getAccession();
			int version = genbankHeaderParser.getVersion();
			String identifier = genbankHeaderParser.getIdentifier();
			AccessionID accession = new AccessionID(id , DataSource.GENBANK, version, identifier);
			sequence.setAccession(accession);
			
			// add features to new sequence
			genbankParser.getFeatures().values().stream()
			.flatMap(List::stream)
			.forEach(sequence::addFeature);

			// add taxonomy ID to new sequence
			List<DBReferenceInfo> dbQualifier = genbankParser.getDatabaseReferences().get("db_xref");
			if (dbQualifier != null){
				DBReferenceInfo q = dbQualifier.get(0);
				sequence.setTaxonomy(new TaxonomyID(q.getDatabase()+":"+q.getId(), DataSource.GENBANK));
			}

			sequences.put(sequence.getAccession().getID(), sequence);
		}

		return sequences;
	}

	public void close() {
		try {
			bufferedReader.close();
			this.closed = true;
		} catch (IOException e) {
			logger.error("Couldn't close the reader.", e);
			this.closed = false;
		}
	}
}

