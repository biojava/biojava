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

import java.awt.event.MouseEvent;
import java.util.EventObject;
import java.util.List;

/**
 * An event indicating that a mouse gesture was recognised within a widget that
 * renders sequences.
 *
 * @author Matthew Pocock
 * @since 1.2
 */
public class SequenceViewerEvent extends EventObject {
  private final List path;
  private final Object target;
  private final MouseEvent mouseEvent;
  private final int pos;
  
  /**
   * Construct a SequenceViewerEvent with the given source, target, mouseEvent
   * and path.
   *
   * @param source  the event source, presumably a GUI component
   * @param target  an Object that is the target of the gesture - a feature, or
   *                {alignment, label, index} or some other structure, or null
   *                if there is no obvious target
   * @param pos the position (offset) within the sequence
   * @param mouseEvent the MouseEvent that caused this event to be produced
   * @param path  a List of SequenceRenderer instances passed through to reach
   *              this event source
   */
  public SequenceViewerEvent(
    Object source,
    Object target,
    int pos,
    MouseEvent mouseEvent,
    List path
  ) {
    super(source);
    this.target = target;
    this.pos = pos;
    this.mouseEvent = mouseEvent;
    this.path = path;
  }
  
  /**
   * Get the list of SequenceRenderer instances that were passed through to
   * produce this event
   *
   * @return a List of SequenceRenderer instances
   */
  public List getPath() {
    return path;
  }
  
  /**
   * Get the Object that was the target of the mouse gesture or null if the
   * mouse is not gesturing over any recognizable rendered object.
   *
   * @return the Object gestured at by the mouse event
   */
  public Object getTarget() {
    return target;
  }
  
  /**
   * Get the offset within the sequence - the symbol index. This is not
   * guaranteed to be within the legal range of symbol indices.
   *
   * @return the position of the gesture in sequence coordinates
   */
  public int getPos() {
    return pos;
  }
  
  /**
   * Get the mouse event that caused this.
   *
   * @return the MouseEvent that caused this gesture to be noticed
   */
  public MouseEvent getMouseEvent() {
    return mouseEvent;
  }
  
  public String toString() {
    return "[" +
        super.toString() + "]" +
        "target: "+ target +
        "pos: " + pos +
        "mosuseEvent: " + mouseEvent;
  }
}
