package org.biojava.bio.program.ssaha;

import java.io.File;
import java.io.IOException;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Packing;

/**
 * <p>
 * Builder for a data store.
 * </p>
 *
 * @author Matthew Pocock
 */
public interface DataStoreFactory {
  /**
   * Get a pre-built data store associated with a file.
   *
   * @param storeFile  the File to map in as a data store
   * @return the DataStore made by mapping the file
   *
   * @throws IOException if the file could not be mapped
   */
  public DataStore getDataStore(File storeFile)
  throws IOException;
  
  /**
   * Build a new DataStore.
   *
   * @param storeFile  the file to store the data store
   * @param seqDB  the SequenceDB to store in the data store
   * @param packing  the Packing used to bit-encode the sequences
   * @param wordLength the number of symbols per word
   * @param threshold  the number of times a word must appear to be ignored
   *
   * @throws IllegalAlphabetException if the packing does not agree with
   *         the sequences
   * @throws BioException if there is a problem building the data store
   */
  public DataStore buildDataStore(
    File storeFile,
    SequenceDB seqDB,
    Packing packing,
    int wordLength,
    int threshold
  ) throws
    IllegalAlphabetException,
    IOException,
    BioException;
}
