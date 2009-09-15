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
package org.biojava.bio.gui.sequence;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author Matthew Pocock
 */
public class SequenceViewerSupport {
  private List listeners = new ArrayList();
  
  public void addSequenceViewerListener(SequenceViewerListener svl) {
    synchronized(listeners) {
      listeners.add(svl);
    }
  }
  
  public void removeSequenceViewerListener(SequenceViewerListener svl) {
    synchronized(listeners) {
      listeners.remove(svl);
    }
  }
  
  public void fireMouseClicked(SequenceViewerEvent sve) {
    List l;
    synchronized(listeners) {
      l = new ArrayList(listeners);
    }
    for(Iterator i = l.iterator(); i.hasNext(); ) {
      SequenceViewerListener svl = (SequenceViewerListener) i.next();
      svl.mouseClicked(sve);
    }
  }
  
  public void fireMousePressed(SequenceViewerEvent sve) {
    List l;
    synchronized(listeners) {
      l = new ArrayList(listeners);
    }
    for(Iterator i = l.iterator(); i.hasNext(); ) {
      SequenceViewerListener svl = (SequenceViewerListener) i.next();
      svl.mousePressed(sve);
    }
  }
  
  public void fireMouseReleased(SequenceViewerEvent sve) {
    List l;
    synchronized(listeners) {
      l = new ArrayList(listeners);
    }
    for(Iterator i = l.iterator(); i.hasNext(); ) {
      SequenceViewerListener svl = (SequenceViewerListener) i.next();
      svl.mouseReleased(sve);
    }
  }
}
