package org.biojava.utils;

import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;

/**
 *
 *
 * @author Matthew Pocock
 */
public class MergingIterator
        implements Iterator
{
  private final Iterator sourceIt;
  private Iterator currentIt;
  private Object nextVal;

  public MergingIterator(Iterator sourceIt) {
    this.sourceIt = sourceIt;
    currentIt = Collections.EMPTY_SET.iterator();
    nextVal = findNextVal();
  }

  public boolean hasNext() {
    return nextVal != null;
  }

  public Object next() {
    Object val = nextVal;
    nextVal = findNextVal();
    return val;
  }

  public void remove() {
    throw new UnsupportedOperationException();
  }

  private Object findNextVal() {
    while(true) {
      if(!currentIt.hasNext()) {
        if(!sourceIt.hasNext()) {
          return null;
        }

        currentIt = ((Collection) (sourceIt.next())).iterator();

        continue;
      } else {
        return currentIt.next();
      }
    }
  }
}
