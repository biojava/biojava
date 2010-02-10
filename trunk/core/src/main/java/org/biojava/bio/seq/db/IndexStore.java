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

package org.biojava.bio.seq.db;

import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.SequenceBuilderFactory;
import org.biojava.bio.seq.io.SequenceFormat;
import org.biojava.bio.seq.io.SymbolTokenization;

/**
 * This defines the objects that IndexedSequenceDB uses to store all of the
 * database state, such as name, format, sequence builder and the actual file
 * offsets.
 *
 * You should only ever have to instantiate an implementation of this. Don't
 * call any of the inerface methods directly. They are intended for
 * IndexedSequenceDB, not you.
 *
 * In general, these objects should be transactional. Calls to store should add
 * the index to temporary storage. When commit is invoked, these indices should
 * all be added to the permanent storage. When rollback is invoked, these
 * indexes should be discarded. If commit fails for any reason, it should leave
 * the permanent storage in the pre-commit status.

 * @author Matthew Pocock
 */
public interface IndexStore {
  /**
   * Add the Index to the store.
   * <p>
   * This method should be transactional. If the store fails, the IndexStore
   * should be left in its original state.
   * <p>
   * If the file of the Index is not known yet, it is the
   * responsibility of the IndexStore to add it to the set returned by
   * getFiles.
   *
   * @param indx the Index to add
   * @throws IllegalIDException if the index has an invalid ID field
   * @throws BioException if the store failed
   */
  void store(Index indx) throws IllegalIDException, BioException;
  
  /**
   * Commit the stored indices to permanent storage.
   *
   * @throws BioException if for any reason the commit fails
   */
  void commit() throws BioException;
  
  /**
   * Discard all uncommited changes.
   */
  void rollback();
  
  /**
   * Fetch an Index based upon an ID.
   *
   * @param id  The ID of the sequence Index to retrieve
   * @throws IllegalIDException if the ID couldn't be found
   * @throws BioException if the fetch fails in the underlying storage mechanism
   */
  Index fetch(String id) throws IllegalIDException, BioException;
  
  /**
   * Retrieve the name of this store. This will be reflected as the name of the
   * IndexedSequenceDB.
   *
   * @return the String name of the index
   */
  String getName();
  
  /**
   * Retrieve the set of all current IDs.
   * <p>
   * This set should either be immutable, or modifiable totally
   * separately from the IndexStore.
   *
   * @return a Set of all legal IDs
   */
  Set getIDs();
  
  /**
   * Retrieve the Set of files that are currently indexed.
   */
  Set getFiles();
  
  /**
   * Retrieve the format of the index file.
   * <p>
   * This set should either be immutable, or modifiable totally
   * separately from the IndexStore.
   *
   * @return a Set of all indexed files
   */
  SequenceFormat getFormat();
  
  /**
   * Retrieve the SequenceBuilderFactory used to build Sequence instances.
   *
   * @return the associated SequenceBuilderFactory
   */
  SequenceBuilderFactory getSBFactory();
  
  /**
   * Retrieve the symbol parser used to turn the sequence characters
   * into Symobl objects.
   *
   * @return the associated SymbolParser
   */
  SymbolTokenization getSymbolParser();
}
