package org.biojava.bio.program.ssaha;

import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.SymbolList;

/**
 * A repository that can be searched with a sequence.
 *
 * @author Matthew Pocock
 */
public interface DataStore {
  /**
   * The alphabet of symbol lists that can be searched against this
   * DataStore.
   *
   * @return a FiniteAlphabet search types of SymbolList
   */
  public FiniteAlphabet getAlphabet();
  
  /**
   * Search the DataStore with a symbol list.
   *
   * @param id  the ID to report the symbol list by e.g. 'test' or 'foo1'
   * @param symList  the symbol list to search with
   * @param listener  the listener to inform of hits
   *
   * @throws IllegalAlphabetException if the symbol list is of a type that
   *         is not compatible with this data store
   */
  public void search(String id, SymbolList symList, SearchListener listener)
  throws IllegalAlphabetException, SearchException;
  
  /**
   * Resolve an ID to a sequence name.
   *
   * @param id  the int number of the sequence name to resolve
   * @return the name of that sequence as a String
   * @throws IndexOutOfBoundsException if id is negative or too large
   */
  public String seqNameForID(int id)
  throws IndexOutOfBoundsException, SearchException;
}
