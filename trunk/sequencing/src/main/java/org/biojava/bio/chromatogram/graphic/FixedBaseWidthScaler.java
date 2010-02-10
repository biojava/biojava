package org.biojava.bio.chromatogram.graphic;

import org.biojava.bio.BioError;
import org.biojava.bio.chromatogram.Chromatogram;
import org.biojava.bio.chromatogram.ChromatogramTools;

/**
 * A {@link ChromatogramNonlinearScaler} that scales all the 
 * base calls in a chromatogram to the same width in pixels,
 * optionally biasing the peak of the call to the center.
 *
 * @author Rhett Sutphin (<a href="http://genome.uiowa.edu/">UI CBCB</a>)
 * @author Matthew Pocock
 * @since 1.3
 */
public class FixedBaseWidthScaler implements ChromatogramNonlinearScaler {
    /** Set to true to get copious, cryptic debugging output on System.out */
    private static final boolean DEBUG = false;
    
    private Chromatogram lastC;
    private int lastSeqLength;
    private final float baseWidth;
    private final boolean centerPeaks;
    private float[] scales;

    /**    
     * Creates a new scaler that will scale bases to the specified width
     * without attempting to center their peaks.
     * @param width the desired call width in pixels
     */
    public FixedBaseWidthScaler(float width) {
        this(width, false);
    }

    /**
     * Creates a new scaler that will scale bases to the specified width
     * and may or may not bias the peaks to the center.
     * @param width the desired call width in pixels
     * @param centerPeaks if true, the scaler will try to put the peak of
     *        in the center of the scaled call.  Otherwise, the whole call
     *        will be scaled using the same factor.
     */
    public FixedBaseWidthScaler(float width, boolean centerPeaks) {
        baseWidth = width;
        this.centerPeaks = centerPeaks;
        lastC = null;
        scales = null;
        lastSeqLength = -1;
    }

    public float scale(Chromatogram c, int traceSampleIndex) throws IndexOutOfBoundsException {
        if (traceSampleIndex < 0 || traceSampleIndex >= c.getTraceLength())
            throw new IndexOutOfBoundsException("Requested a scale of a trace sample outside of the chromatogram");
        synchronized (this) {
            calcAllScales(c); 
            if (traceSampleIndex >= scales.length) 
                throw new BioError("This shouldn't happen: a valid trace sample is not included in the calculated scales array.  This is probably due to a failure in dirty checking.");
            return scales[traceSampleIndex];
        }
    }
    
    /**
     * Calculates all the scaled x-coordinates for the given chromatogram,
     * but only if it isn't the one that we already have the scales for.
     *
     * @param c  the Chromatogram to calculate the scale for
     */
    private synchronized void calcAllScales(Chromatogram c) {
        if (c == null) return;
        if (scales != null 
            && lastC != null && lastC == c 
            && scales.length == c.getTraceLength()
            && lastSeqLength == c.getSequenceLength()) 
            return;
        if (scales == null || scales.length != c.getTraceLength())
            scales = new float[c.getTraceLength()];
        if (DEBUG) System.out.println("bw=" + baseWidth + " cp="+centerPeaks);
        int[] peaks = ChromatogramTools.getTraceOffsetArray(c);
        int left = 0;
        int center = peaks[0];
        int right = peaks[0] + (int) Math.floor( ((double) (peaks[1] - peaks[0])) / 2 );
        if (DEBUG) System.out.print("b=1 ");
        scaleBase(left, center, right, 0.0f, false);
        for (int i = 1 ; i < peaks.length - 1 ; i++) {
            left = right + 1;
            center = peaks[i];
            right = peaks[i] + (int) Math.floor( ((double) (peaks[i+1] - peaks[i])) / 2 );
            if (DEBUG) System.out.print("b="+(i+1)+" ");
            scaleBase(left, center, right, i * baseWidth, false);
        }
        left = right + 1;
        center = peaks[peaks.length - 1];
        right = c.getTraceLength() - 1;
        if (DEBUG) System.out.print("b="+(peaks.length)+" ");
        scaleBase(left, center, right, (peaks.length - 1) * baseWidth, true);
        
        lastC = c;
        lastSeqLength = c.getSequenceLength();
    }

    /**
     * Calculate the scaled x coordinates for a base call.  The
     * base call consists of the samples in the range <code>[left, right]</code> 
     * and the peak is at <code>center</code> (which must be in that same range).
     * The parameter <code>widthSoFar</code> provides the x coordinate for
     * the sample at <code>left</code> and this method must return a similar
     * value for the base call immediately following this one (i.e., the one 
     * whose leftmost sample would be <code>right+1</code>).
     * @param left the leftmost sample index of the base (inclusive)
     * @param center the peak sample index
     * @param right the rightmost sample index of the base (inclusive)
     * @param nextX the total width of all the scaled bases up to this one
     * @param lastBase true if this is the last base, false otherwise.  The
     *        last base must be handled specially so that the last sample is 
     *        at baseWidth * c.sequenceLength.
     */
    private void scaleBase(int left, int center, int right, float nextX, boolean lastBase) {
        if (left > center)  throw new BioError("Assertion failure: left > center ; l="+left+"; c="+center+"; r="+right);
        if (center > right) throw new BioError("Assertion failure: center > right ; l="+left+"; c="+center+"; r="+right);
        if (left < 0)   throw new BioError("Assertion failure: left < 0 ; l="+left+"; c="+center+"; r="+right);
        if (center < 0) throw new BioError("Assertion failure: center < 0 ; l="+left+"; c="+center+"; r="+right);
        if (right < 0)  throw new BioError("Assertion failure: right < 0 ; l="+left+"; c="+center+"; r="+right);
        if (DEBUG) System.out.println("l="+left+" c="+center+" r="+right+" nX="+nextX);
        if (centerPeaks) {
            float leftIncrement = (baseWidth * 0.5f) / (center - left);
            if (DEBUG) System.out.println("  lside ("+leftIncrement+"): ");
            for (int i = left ; i < center ; i++) {
                scales[i] = nextX;
                if (DEBUG) System.out.println("    s["+i+"]="+scales[i]);
                nextX += leftIncrement;
            }
            float rightIncrement = (baseWidth * 0.5f) / (right - center + (lastBase?0:1));
            if (DEBUG) System.out.println("  rside ("+rightIncrement+"): ");
            for (int i = center ; i <= right ; i++) {
                scales[i] = nextX;
                if (DEBUG) System.out.println("    s["+i+"]="+scales[i]);
                nextX += rightIncrement;
            }
        }
        else {
            float increment = baseWidth / (right - left + (lastBase?0:1));
            if (DEBUG) System.out.println("  noside ("+increment+"): ");
            for (int i = left ; i <= right ; i++) {
                scales[i] = nextX;
                if (DEBUG) System.out.println("    "+scales[i]);
                nextX += increment;
            }
        }
    }
} 
