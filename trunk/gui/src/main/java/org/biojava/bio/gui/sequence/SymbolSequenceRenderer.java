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

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.event.MouseEvent;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.util.List;

import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SymbolList;

/**
 * <code>SymbolSequenceRenderer</code> renders symbols of a
 * <code>SymbolList</code>.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author David Huen
 * @author Keith James
 * @author Kalle Näslund
 */
public class SymbolSequenceRenderer implements SequenceRenderer
{
  private double depth = 25.0;
  private Paint  outline;

  public SymbolSequenceRenderer()
  {
    outline = Color.black;
  }

  public double getDepth(SequenceRenderContext context)
  {
    return depth + 1.0;
  }

  public double getMinimumLeader(SequenceRenderContext context)
  {
    return 0.0;
  }

  public double getMinimumTrailer(SequenceRenderContext context)
  {
    return 0.0;
  }

  public void paint(Graphics2D g2, SequenceRenderContext context)
  {
    Rectangle2D prevClip = g2.getClipBounds();
    AffineTransform prevTransform = g2.getTransform();

    g2.setPaint(outline);

    Font font = context.getFont();

    Rectangle2D maxCharBounds =
            font.getMaxCharBounds(g2.getFontRenderContext());

    double scale = context.getScale();

    if (scale >= (maxCharBounds.getWidth() * 0.3) &&
            scale >= (maxCharBounds.getHeight() * 0.3))
    {
      double xFontOffset = 0.0;
      double yFontOffset = 0.0;

      // These offsets are not set quite correctly yet. The
      // Rectangle2D from getMaxCharBounds() seems slightly
      // off. The "correct" application of translations based on
      // the Rectangle2D seem to give the wrong results. The
      // values below are mostly fudges.
      if (context.getDirection() == SequenceRenderContext.HORIZONTAL)
      {
        xFontOffset = maxCharBounds.getCenterX() * 0.25;
        yFontOffset = - maxCharBounds.getCenterY() + (depth * 0.5);
      }
      else
      {
        xFontOffset = - maxCharBounds.getCenterX() + (depth * 0.5);
        yFontOffset = - maxCharBounds.getCenterY() * 3.0;
      }

      SymbolList seq = context.getSymbols();
      SymbolTokenization toke = null;
      try {
        toke = seq.getAlphabet().getTokenization("token");
      } catch (Exception ex) {
        throw new BioRuntimeException(ex);
      }

      Location visible = GUITools.getVisibleRange(context, g2);
      for (int sPos = visible.getMin(); sPos <= visible.getMax(); sPos++)
      {
        double gPos = context.sequenceToGraphics(sPos);
        String s = "?";
        try {
          s = toke.tokenizeSymbol(seq.symbolAt(sPos));
        } catch (Exception ex) {
          // We'll ignore the case of not being able to tokenize it
        }

        if (context.getDirection() == SequenceRenderContext.HORIZONTAL)
        {
          g2.drawString(s,
                        (float) (gPos + xFontOffset),
                        (float) yFontOffset);
        }
        else
        {
          g2.drawString(s,
                        (float) xFontOffset,
                        (float) (gPos + yFontOffset));
        }
      }
    }

    g2.setClip(prevClip);
    g2.setTransform(prevTransform);
  }

  public SequenceViewerEvent processMouseEvent(SequenceRenderContext context,
                                               MouseEvent            me,
                                               List                  path)
  {
    path.add(this);
    int sPos = context.graphicsToSequence(me.getPoint());
    return new SequenceViewerEvent(this, null, sPos, me, path);
  }
}
