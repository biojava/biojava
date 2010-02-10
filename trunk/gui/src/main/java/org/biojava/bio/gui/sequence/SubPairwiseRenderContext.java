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

import java.awt.Font;
import java.awt.Point;
import java.awt.geom.Point2D;

import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SymbolList;

/**
 * <p><code>SubPairwiseRenderContext</code> is a rendering context
 * which wraps a delegate context and effectively hides some of the
 * delegate's properties with its own. If any of the
 * <code>SymbolList</code>, <code>FeatureHolder</code> or
 * <code>RangeLocation</code> arguments are not null, their values are
 * returned. Otherwise the delegate's method is called and its return
 * value is returned instead.</p>
 *
 * @author Keith James
 * @author Matthew Pocock
 * @since 1.2
 */
public class SubPairwiseRenderContext implements PairwiseRenderContext
{
    private final PairwiseRenderContext context;
    private final SymbolList            symbols;
    private final SymbolList            secondarySymbols;
    private final FeatureHolder         features;
    private final FeatureHolder         secondaryFeatures;
    private final RangeLocation         range;
    private final RangeLocation         secondaryRange;

    /**
     * Creates a new <code>SubPairwiseRenderContext</code>.
     *
     * @param context a <code>PairwiseRenderContext</code> to
     * wrap. This should not be null.
     * @param symbols a <code>SymbolList</code> to use instead of the
     * delegate's. May be null.
     * @param secondarySymbols a <code>SymbolList</code> to use
     * instead of the delegate's. May be null.
     * @param features a <code>FeatureHolder</code> to use instead of
     * the delegate's. May be null.
     * @param secondaryFeatures a <code>FeatureHolder</code> to use
     * instead of the delegate's. May be null.
     * @param range a <code>RangeLocation</code> to use instead of the
     * delegate's. May be null.
     * @param secondaryRange a <code>RangeLocation</code> to use
     * instead of the delegate's. May be null.
     */
    public SubPairwiseRenderContext(PairwiseRenderContext context,
                                    SymbolList            symbols,
                                    SymbolList            secondarySymbols,
                                    FeatureHolder         features,
                                    FeatureHolder         secondaryFeatures,
                                    RangeLocation         range,
                                    RangeLocation         secondaryRange)
    {
        this.context           = context;
        this.symbols           = symbols;
        this.secondarySymbols  = secondarySymbols;
        this.features          = features;
        this.secondaryFeatures = secondaryFeatures;
        this.range             = range;
        this.secondaryRange    = secondaryRange;
    }

    public SymbolList getSymbols()
    {
        if (symbols == null)
        {
            return context.getSymbols();
        }
        else
        {
            return symbols;
        }
    }

    public SymbolList getSecondarySymbols()
    {
        if (secondarySymbols == null)
        {
            return context.getSecondarySymbols();
        }
        else
        {
            return secondarySymbols;
        }
    }

    public FeatureHolder getFeatures()
    {
        if (features == null)
        {
            return context.getFeatures();
        }
        else
        {
            return features;
        }
    }

    public FeatureHolder getSecondaryFeatures()
    {
        if (secondaryFeatures == null)
        {
            return context.getSecondaryFeatures();
        }
        else
        {
            return secondaryFeatures;
        }
    }

    public RangeLocation getRange()
    {
        if (range == null)
        {
            return context.getRange();
        }
        else
        {
            return range;
        }
    }

    public RangeLocation getSecondaryRange()
    {
        if (secondaryRange == null)
        {
            return context.getSecondaryRange();
        }
        else
        {
            return secondaryRange;
        }
    }

    public int getDirection()
    {
        return context.getDirection();
    }

    public int getSecondaryDirection()
    {
        return context.getSecondaryDirection();
    }

    public double getScale()
    {
        return context.getScale();
    }

    public SequenceRenderContext.Border getLeadingBorder()
    {
        return context.getLeadingBorder();
    }

    public SequenceRenderContext.Border getTrailingBorder()
    {
        return context.getTrailingBorder();
    }

    public double sequenceToGraphics(int sequencePos)
    {
        return context.sequenceToGraphics(sequencePos);
    }

    public double secondarySequenceToGraphics(int sequencePos)
    {
        return context.secondarySequenceToGraphics(sequencePos);
    }

    public int graphicsToSequence(double graphicsPos)
    {
        return context.graphicsToSequence(graphicsPos);
    }

    public int graphicsToSequence(Point2D point)
    {
        return context.graphicsToSequence(point);
    }

    public int graphicsToSecondarySequence(double graphicsPos)
    {
        return context.graphicsToSecondarySequence(graphicsPos);
    }

    public int graphicsToSecondarySequence(Point point)
    {
        return context.graphicsToSecondarySequence(point);
    }

    public Font getFont()
    {
        return context.getFont();
    }
}
