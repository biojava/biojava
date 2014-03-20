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
 */

package org.biojava3.ontology.utils;

import java.lang.ref.Reference;
import java.lang.ref.WeakReference;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * <p>
 * A utility class to provide management for informing ChangeListeners of
 * ChangeEvents.
 * </p>
 *
 * <p>
 * This is loosely modelled after the standard PropertyChangeEvent objects.
 * </p>
 *
 * <p>
 * For an object to correctly fire these events, they must follow a broad
 * outline like this:
 * <code><pre>
 * public void mutator(foo arg) throw ChangeVetoException {
 *   ChangeEvent cevt = new ChangeEvent(this, SOME_EVENT_TYPE, arg);
 *   synchronized(changeSupport) {
 *     changeSupport.firePreChangeEvent(cevt);
 *     // update our state using arg
 *     // ...
 *     changeSupport.firePostChangeEvent(cevt);
 *   }
 * }
 * </pre></code>
 * </p>
 *
 * <p>
 * The methods that delegate adding and removing listeners to a ChangeSupport
 * must take responsibility for synchronizing on the delegate.
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Keith James (docs)
 * @author Kalle Naslund (tiny bugfix)
 * @since 1.1
 */

public class ChangeSupport {
  private int listenerCount;
  private int delta;
  private Set unchanging;
  private Reference[] listeners;
  private ChangeType[] types;

  /**
   * Return true if we have any listeners registered at all.
   * 
   * @return true if there are listeners
   */
  public boolean hasListeners() {
      return (listenerCount > 0);
  }
  
  /**
   * Return true if we have listeners registered for a particular change type.
   * 
   * @param ct	the ChangeType to check
   * @return	true if there are listeners for this type
   */
  public boolean hasListeners(ChangeType ct)
  {
  	for(int i = 0; i < listenerCount ; i++ ) {
  	  ChangeType type = (ChangeType) types[i];
  	  if(ct.isMatchingType(type)) {
  	  	return true;
  	  }
  	}
  	
  	return false;
  }

  /**
   * Generate a new ChangeSupport instance.
   */
  public ChangeSupport() {
    this(5);
  }

  /**
   * Generate a new ChangeSupport instance which has room for initialSize
   * listeners before it needs to grow any resources.
   *
   * @param initialSize  the number of listeners that can be added before this
   *                     needs to grow for the first time
   */
  public ChangeSupport(int initialSize) {
    this(initialSize, 5);
  }

  /**
   * Generate a new ChangeSupport instance which has room for initialSize
   * listeners before it needs to grow any resources, and which will grow by
   * delta each time.
   *
   * @param initialSize  the number of listeners that can be added before this
   *                     needs to grow for the first time
   * @param delta  the number of listener slots that this will grow by each time
   *               it needs to
   */
  public ChangeSupport(int initialSize, int delta) {
    this(Collections.EMPTY_SET, initialSize, delta);
  }

  public ChangeSupport(Set unchanging) {
    this(unchanging, 0, 5);
  }

  /**
   * Generate a new ChangeSupport instance which has room for initialSize
   * listeners before it needs to grow any resources, and which will grow by
   * delta each time.
   *
   * @param unchanging Set of ChangeTypes that can never be fired
   * @param initialSize  the number of listeners that can be added before this
   *                     needs to grow for the first time
   * @param delta  the number of listener slots that this will grow by each time
   *               it needs to
   */
  public ChangeSupport(Set unchanging, int initialSize, int delta) {
    this.listenerCount = 0;
    this.listeners = new Reference[initialSize];
    this.types = new ChangeType[initialSize];

    this.delta = delta;
    this.unchanging = new HashSet(unchanging);
  }
  /**
   * Add a listener that will be informed of all changes.
   *
   * @param cl  the ChangeListener to add
   */
  public void addChangeListener(ChangeListener cl) {
    addChangeListener(cl, ChangeType.UNKNOWN);
  }

  /**
   * Add a listener that will be informed of changes of a given type (and it's subtypes)
   *
   * @param cl  the ChangeListener
   * @param ct  the ChangeType it is to be informed of
   */
  public void addChangeListener(ChangeListener cl, ChangeType ct) {
    if (ct == null) {
      throw new NullPointerException("Since 1.2, listeners registered for the null changetype are not meaningful.  Please register a listener for ChangeType.UNKNOWN instead");
    }

    if(isUnchanging(ct)) {
      return;
    }

    synchronized(this) {
      growIfNecessary();
      types[listenerCount] = ct;
      listeners[listenerCount] = new WeakReference(cl);
      listenerCount++;
    }
  }

  /**
   * Grows the internal resources if by adding one more listener they would be
   * full.
   */
  protected void growIfNecessary() {
    //try cleaning up first
    synchronized(this){
        reapGarbageListeners();
    }  
    if(listenerCount == listeners.length) {
      int newLength = listenerCount + delta;
      Reference[] newList = new Reference[newLength];
      ChangeType[] newTypes = new ChangeType[newLength];

      System.arraycopy(listeners, 0, newList, 0, listenerCount);
      System.arraycopy(types, 0, newTypes, 0, listenerCount);

      listeners = newList;
      types = newTypes;
    }
  }

  /**
   * Remove a listener that was interested in all types of changes.
   *
   * @param cl  a ChangeListener to remove
   */
  public void removeChangeListener(ChangeListener cl) {
    removeChangeListener(cl, ChangeType.UNKNOWN);
  }

  /**
   * Remove a listener that was interested in a specific types of changes.
   *
   * @param cl  a ChangeListener to remove
   * @param ct  the ChangeType that it was interested in
   */
  public void removeChangeListener(ChangeListener cl, ChangeType ct) {
    synchronized(this) {
      for(int i = 0; i < listenerCount; i++) {
        if( (listeners[i].get() == cl) && (types[i] == ct) ) {
          listenerCount--;
          System.arraycopy(listeners, i+1, listeners, i, (listenerCount - i));
          System.arraycopy(types, i+1, types, i, (listenerCount - i));
          return;
        }
      }
    }
  }

    /**
     * Remove all references to listeners which have been cleared by the
     * garbage collector.  This method should only be called when the
     * object is locked.
     */

    protected void reapGarbageListeners() {
	int pp = 0;
	for (int p = 0; p < listenerCount; ++p) {
	    Reference r = listeners[p];
	    if (r.get() != null) {
	        types[pp] = types[p];
		listeners[pp] = r;
		pp++;
	    }else{ //if it is null release the reference
                r = null;
            }
	}
	listenerCount = pp;
    }

  /**
   * <p>
   * Inform the listeners that a change is about to take place using their
   * firePreChangeEvent methods.
   * </p>
   *
   * <p>
   * Listeners will be informed if they were interested in all types of event,
   * or if ce.getType() is equal to the type they are registered for.
   * </p>
   *
   * <p>
   * This method must be called while the current thread holds the lock on this change support.
   * </p>
   * 
   * @param ce  the ChangeEvent to pass on
   * @throws ChangeVetoException if any of the listeners veto this change
   */
  public void firePreChangeEvent(ChangeEvent ce)
  throws ChangeVetoException {
     assert Thread.holdsLock(this)
            : "firePreChangeEvent must be called in a synchronized block locking the ChangeSupport";
    boolean needToReap = false;

    ChangeType ct = ce.getType();
    int listenerCount = this.listenerCount;
    ChangeType[] types = new ChangeType[listenerCount];
    System.arraycopy(this.types, 0, types, 0, listenerCount);

    Reference[] listeners = new Reference[listenerCount];
    System.arraycopy(this.listeners, 0, listeners, 0, listenerCount);

    for(int i = 0; i < listenerCount; i++) {
      ChangeType lt = types[i];
      if( ct.isMatchingType(lt)) {
        ChangeListener cl = (ChangeListener) listeners[i].get();
        if (cl != null) {
        	synchronized (cl) {
        		cl.preChange(ce);	
			}
          
        } else {
          needToReap = true;
        }
      }
    }

    if (needToReap) {
      reapGarbageListeners();
    }
  }

  /**
   * <p>
   * Inform the listeners that a change has taken place using their
   * firePostChangeEvent methods.
   * </p>
   *
   * <p>
   * Listeners will be informed if they were interested in all types of event,
   * or if ce.getType() is equal to the type they are registered for.
   * </p>
   *
   * <p>
   * This method must be called while the current thread holds the lock on this change support.
   * </p>
   *
   * @param ce  the ChangeEvent to pass on
   */

  public void firePostChangeEvent(ChangeEvent ce) {
    assert Thread.holdsLock(this)
            : "firePostChangeEvent must be called in a synchronized block locking the ChangeSupport";
    boolean needToReap = false;

    ChangeType ct = ce.getType();
    int listenerCount = this.listenerCount;
    ChangeType[] types = new ChangeType[listenerCount];
    System.arraycopy(this.types, 0, types, 0, listenerCount);

    Reference[] listeners = new Reference[listenerCount];
    System.arraycopy(this.listeners, 0, listeners, 0, listenerCount);

    for(int i = 0; i < listenerCount; i++) {
      ChangeType lt = types[i];
      if( ct.isMatchingType(lt) ) {
        ChangeListener cl = (ChangeListener) listeners[i].get();
        if (cl != null) {
          cl.postChange(ce);
        } else {
          needToReap = true;
        }
      }
    }

    if (needToReap) {
      reapGarbageListeners();
    }
  }

  public boolean isUnchanging(ChangeType ct) {
    if(unchanging == null) {
      return false;
    }

    for(Iterator i = ct.matchingTypes(); i.hasNext(); ) {
      if(unchanging.contains(i.next())) {
        return true;
      }
    }

    return false;
  }

  public String displayString()
  {
    StringBuffer sb = new StringBuffer();
    sb.append(this.toString());
    sb.append("\n");
    for(int i = 0; i < listenerCount; i++) {
      sb.append("\t");
      sb.append(listeners[i].get());
      sb.append("\t");
      sb.append(types[i]);
      sb.append("\n");
    }

    return sb.toString();
  }
}
