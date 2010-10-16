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
import java.awt.Paint;

import org.biojava.bio.seq.FeatureFilter;
import org.biojava.utils.ChangeType;

/**
 * A concrete AbstractPeptideDigestRenderer. The features matching the given FeatureFilter
 * are rendered as blue arrows.
 *
 * @author Mark Southern
 * @since 1.5
 */
public class PeptideDigestRenderer extends AbstractPeptideDigestRenderer {
    public static final ChangeType DIGEST = new ChangeType("The peptide digest has changed",
        "org.biojava.bio.gui.sequence.PeptideDigestRenderer", "DIGEST",
        SequenceRenderContext.REPAINT
    );

    private Paint defaultPaint = Color.BLUE;
    
    public PeptideDigestRenderer(FeatureSource source) {
        super(source);
    }

    public PeptideDigestRenderer(FeatureSource source, FeatureFilter filter) {
        super(source,filter);
    }
    
    public PeptideDigestRenderer(FeatureSource source, FeatureFilter filter, int distanceBetweenFeatures)
    {
        super(source,filter,distanceBetweenFeatures);
    }

    public void setDefaultPaint(Paint p) {
        this.defaultPaint = p;
    }

    public Paint getDefaultPaint() {
        return this.defaultPaint;
    }
    
    public FeatureRenderer createRenderer(int lane)
    {
    	ArrowedFeatureRenderer fr = new ArrowedFeatureRenderer();
        fr.setFill(defaultPaint);
        fr.setOutline(defaultPaint);
        fr.setArrowHeadSize(5);
        fr.setArrowScoop(2);
        fr.setArrowSize(1);
        return fr;
    }
}
