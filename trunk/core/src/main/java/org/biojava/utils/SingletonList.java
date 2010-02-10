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
package org.biojava.utils;

import java.io.Serializable;
import java.util.AbstractList;

/**
 * @author Matthew Pocock
 */
public class SingletonList extends AbstractList implements Serializable {
  private final Object obj;
  
  public SingletonList(Object obj) {
    this.obj = obj;
  }
  
  public int size() {
    return 1;
  }
  
  public Object get(int i) throws IndexOutOfBoundsException {
    if(i == 0) {
      return obj;
    } else {
      throw new IndexOutOfBoundsException("Can't access item " + i + " of 1");
    }
  }
}
