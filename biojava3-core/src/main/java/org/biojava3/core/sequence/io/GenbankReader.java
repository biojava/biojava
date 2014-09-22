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
 * @author Karl Nicholas <github:karlnicholas>
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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedHashMap;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.io.template.SequenceHeaderParserInterface;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Use GenbankReaderHelper as an example of how to use this class where GenbankReaderHelper should be the
 * primary class used to read Genbank files
 *
 */
public class GenbankReader<S extends AbstractSequence<C>, C extends Compound> {

	private final static Logger logger = LoggerFactory.getLogger(GenbankReader.class);

    private SequenceCreatorInterface<C> sequenceCreator;
    private GenbankSequenceParser<S,C> genbankParser;
    private InputStream inputStream;
    
    /**
     * If you are going to use FileProxyProteinSequenceCreator then do not use this constructor because we need details about
     * local file offsets for quick reads. InputStreams does not give you the name of the stream to access quickly via file seek. A seek in
     * an inputstream is forced to read all the data so you don't gain anything.
     * @param br
     * @param headerParser
     * @param sequenceCreator
     */
    public GenbankReader(InputStream is, SequenceHeaderParserInterface<S,C> headerParser, SequenceCreatorInterface<C> sequenceCreator) {
        this.sequenceCreator = sequenceCreator;
        this.inputStream = is;
    	genbankParser = new GenbankSequenceParser<S,C>();
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
    		File file, 
    		SequenceHeaderParserInterface<S,C> headerParser, 
    		SequenceCreatorInterface<C> sequenceCreator
    		) throws FileNotFoundException {
    	
        inputStream = new FileInputStream(file);
        this.sequenceCreator = sequenceCreator;
    	genbankParser = new GenbankSequenceParser<S,C>();
    }

    /**
     * The parsing is done in this method.<br>
     * This method tries to process all the available Genbank records 
     * in the File or InputStream, closes the underlying resource, 
     * and return the results in {@link LinkedHashMap}.<br>
     * You don't need to call {@link #close()} after calling this method.
     * @see #process(int)
     * @return {@link HashMap} containing all the parsed Genbank records 
     * present, starting current fileIndex onwards.
     * @throws Exception 
     */
    public LinkedHashMap<String,S> process() throws Exception {
    	LinkedHashMap<String,S> sequences = process(-1);
    	return sequences;
    }

    /**
     * This method tries to parse maximum <code>max</code> records from
     * the open File or InputStream, and leaves the underlying resource open.<br>
     * Subsequent calls to the same method continue parsing the rest of the file.<br>
     * This is particularly useful when dealing with very big data files,
     * (e.g. NCBI nr database), which can't fit into memory and will take long
     * time before the first result is available.<br>
     * <b>N.B.</b>
     * <ul>
     * <li>This method ca't be called after calling its NO-ARGUMENT twin.</li> 
     * <li>remember to close the underlying resource when you are done.</li> 
     * </ul>
     * @see #process()
     * @author Amr AL-Hossary
     * @since 3.0.6
     * @param max maximum number of records to return, <code>-1</code> for infinity.
     * @return {@link HashMap} containing maximum <code>max</code> parsed Genbank records 
     * present, starting current fileIndex onwards.
     * @throws Exception 
     */
    public LinkedHashMap<String,S> process(int max) throws Exception {
        LinkedHashMap<String,S> sequences = new LinkedHashMap<String,S>();
		@SuppressWarnings("unchecked")
		S sequence = (S) sequenceCreator.getSequence(genbankParser.getSequence(new BufferedReader(new InputStreamReader(inputStream)), 0), 0);
		genbankParser.getSequenceHeaderParser().parseHeader(genbankParser.getHeader(), sequence);
    	sequences.put(sequence.getAccession().getID(), sequence);
    	close();
        return sequences;
    }

	public void close() throws IOException {
		inputStream.close();
	}

    public static void main(String[] args) throws Exception {
        String proteinFile = "src/test/resources/BondFeature.gb";
        FileInputStream is = new FileInputStream(proteinFile);

        GenbankReader<ProteinSequence, AminoAcidCompound> proteinReader = new GenbankReader<ProteinSequence, AminoAcidCompound>(is, new GenericGenbankHeaderParser<ProteinSequence,AminoAcidCompound>(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
        LinkedHashMap<String,ProteinSequence> proteinSequences = proteinReader.process();
        logger.info("Protein Sequences: {}", proteinSequences);

        String inputFile = "src/test/resources/NM_000266.gb";
        is = new FileInputStream(inputFile);
        GenbankReader<DNASequence, NucleotideCompound> dnaReader = new GenbankReader<DNASequence, NucleotideCompound>(is, new GenericGenbankHeaderParser<DNASequence,NucleotideCompound>(), new DNASequenceCreator(DNACompoundSet.getDNACompoundSet()));
        LinkedHashMap<String,DNASequence> dnaSequences = dnaReader.process();
        is.close();
        logger.info("DNA Sequences: {}", dnaSequences);
    }

}

