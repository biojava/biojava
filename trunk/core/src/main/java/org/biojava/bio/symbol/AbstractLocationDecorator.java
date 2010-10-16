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

package org.biojava.bio.symbol;

import java.util.Iterator;

/**
 * Abstract <code>Location</code> decorator (wrapper).
 *
 * <p>
 * These wrap up a normal <code>Location</code> object, and act as markers to indicate that
 * the location has some kind of special semantics.
 * </p>
 *
 * When implementing a new Location decorator.
 * @author Matthew Pocock
 */
public abstract class AbstractLocationDecorator implements Location {
  private final Location wrapped;

  /**
   * Construct a new decorator wrapping the specified Location.
   */

  protected AbstractLocationDecorator(Location wrapped) {
    this.wrapped = wrapped;
  }

  protected final Location getWrapped() {
    return wrapped;
  }

  protected abstract Location decorate(Location loc);

  public Location newInstance(Location loc) {
    if(loc instanceof AbstractLocationDecorator) {
      Location wrapped = ((AbstractLocationDecorator) loc).getWrapped();
      loc = wrapped.newInstance(wrapped);
    }
    return decorate(loc);
  }

  public Location getDecorator(Class decoratorClass) {
    if(decoratorClass.isInstance(this)) {
      return this;
    } else {
      return getWrapped().getDecorator(decoratorClass);
    }
  }

  public int getMin() {
    return getWrapped().getMin();
  }

  public int getMax() {
    return getWrapped().getMax();
  }

  public boolean overlaps(Location l) {
    return getWrapped().overlaps(l);
  }

  public boolean contains(Location l) {
    return getWrapped().contains(l);
  }

  public boolean contains(int p) {
    return getWrapped().contains(p);
  }

  public boolean equals(Object o) {
    return getWrapped().equals(o);
  }

  public Location intersection(Location l) {
    return getWrapped().intersection(l);
  }

  public Location union(Location l) {
    return getWrapped().union(l);
  }

  public SymbolList symbols(SymbolList seq) {
    return getWrapped().symbols(seq);
  }

  public Location translate(int dist) {
    return decorate(getWrapped().translate(dist));
  }

  public boolean isContiguous() {
    return getWrapped().isContiguous();
  }

  public Iterator blockIterator() {
    return getWrapped().blockIterator();
  }
}
