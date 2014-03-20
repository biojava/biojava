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

import java.util.Collections;
import java.util.Set;

/**
 * This is a utility implementation of Changeable that doesn't fire any events
 * or keep references to any listeners. Use this when you have a final immutable
 * class and can't be bothered to fill in all those method stubs.
 *
 * @author Matthew Pocock
 * @since 1.3
 */
public class Unchangeable
implements Changeable {
  public final void addChangeListener(ChangeListener cl) {}

  public final void addChangeListener(ChangeListener cl, ChangeType ct) {}

  public final Set getListeners(ChangeType ct) { return Collections.EMPTY_SET; }

  public final void removeChangeListener(ChangeListener cl) {}

  public final void removeChangeListener(ChangeListener cl, ChangeType ct) {}

  public final void addForwarder(ChangeForwarder cf, ChangeType ct) {}

  public final void removeForwarder(ChangeForwarder cf, ChangeType ct) {}

  public final Set getForwarders(ChangeType ct) { return Collections.EMPTY_SET; }

  public final boolean isUnchanging(ChangeType ct) { return true; }
}
