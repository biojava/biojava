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

package org.biojava.bio.program.tagvalue;

import java.util.Map;
import java.util.Set;

import org.biojava.utils.SmallMap;

public class Index2Model {
  private final Map keys;
  private String primaryKeyName;
  
  public Index2Model() {
    keys = new SmallMap();
  }
  
  /**
   * <p>
   * Set the tag to use as a primary key in the index.
   * </p>
   *
   * <p>
   * Whenever a value for the primary key tag is seen, this is passed to the
   * indexer as the primary key for indexing.
   * </p>
   *
   * <p>
   * Primary keys must be unique between entries, and each entry must provide
   * exactly one primary key value.
   * </p>
   *
   * @param primaryKeyName the tag to use as primary key
   */
  public void setPrimaryKeyName(String primaryKeyName) {
    this.primaryKeyName = primaryKeyName;
  }
  
  /**
   * Retrieve the tag currently used as primary key.
   *
   * @return a String representing the primary key name
   */
  public String getPrimaryKeyName() {
    return primaryKeyName;
  }
  
  /**
   * <p>
   * Add a key and a path to that key in the tag-value hierachy.
   * </p>
   *
   * <p>
   * Secondary keys are potentialy non-unique properties of the entries being
   * indexed. Multiple records can use the same secondary key values, and a
   * single record can have multiple values for a secondary key. However, the
   * primary key must be unique.
   * </p>
   *
   * @param keyName  the name of the secondary key to add
   * @param path  the names of each tag to follow to reach the value of the key
   */
  public void addKeyPath(String keyName, Object[] path) {
    keys.put(keyName, path);
  }
  
  /**
   * Remove a key.
   *
   * @param keyName  the name of the key to remove
   */
  public void removeKeyPath(String keyName) {
    keys.remove(keyName);
  }
  
  public Object[] getKeyPath(String keyName) {
    return (Object []) keys.get(keyName);
  }

  public Set getKeys() {
    return keys.keySet();
  }
}
