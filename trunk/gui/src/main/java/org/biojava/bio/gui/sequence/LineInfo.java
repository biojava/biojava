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

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Encapsulates the rendering info for a single line of the display.
 * <p>
 * The single line of info may be divided into multiple regions, each rendered
 * by their own SequenceRenderer. It is the job of this class to cache the
 * information about how much space each one wants, and how much space they want
 * in total. A SequenceRenderer or SequencePanel that delegates rendering to
 * multiple child SequenceRenderer instances may want to use these objects
 * for storing this information about each row they are responsible for.
 *
 * @author Matthew Pocock
 */
public class LineInfo {
  private Map rendererToDepth = new HashMap();
  private double totalDepth = 0.0;
  
  public double getDepth(SequenceRenderer r) {
    Double depth = (Double) rendererToDepth.get(r);
    return depth.doubleValue();
  }
  
  public void setDepth(SequenceRenderer r, double depth) {
    totalDepth = Double.NaN;
    
    rendererToDepth.put(r, new Double(depth));
  }
  
  public double getTotalDepth() {
    if(Double.isNaN(totalDepth)) {
      totalDepth = 0.0;
      
      for(Iterator i = rendererToDepth.values().iterator(); i.hasNext(); ) {
        Double d = (Double) i.next(); // should never be null
        totalDepth += d.doubleValue();
      }
    }
    
    return totalDepth;
  }
}
