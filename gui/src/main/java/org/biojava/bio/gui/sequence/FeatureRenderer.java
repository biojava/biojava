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

import java.awt.Graphics2D;
import java.awt.event.MouseEvent;

import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;

/**
 * @author Thomas Down
 * @author Matthew Pocock
 */
public interface FeatureRenderer {
  void renderFeature(
    Graphics2D g2,
    Feature feat,
    SequenceRenderContext context
  );

  double getDepth(SequenceRenderContext src);

  public FeatureHolder processMouseEvent(
    FeatureHolder hits,
    SequenceRenderContext src,
    MouseEvent me
  );
}
