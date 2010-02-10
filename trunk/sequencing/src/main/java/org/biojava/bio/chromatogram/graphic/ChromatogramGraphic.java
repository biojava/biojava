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

package org.biojava.bio.chromatogram.graphic;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.BioError;
import org.biojava.bio.chromatogram.Chromatogram;
import org.biojava.bio.chromatogram.ChromatogramTools;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * Encapsulates a configurable method for drawing a {@link Chromatogram} 
 * into a graphics context.  
 *
 * @author Rhett Sutphin (<a href="http://genome.uiowa.edu/">UI CBCB</a>)
 * @author Matthew Pocock
 * @since 1.3
 */
public class ChromatogramGraphic implements Cloneable {
    /** A pseudo call list for use when a Chromatogram has no base calls. */
    private static SymbolList SINGLE_CALL = 
        new SimpleSymbolList(new Symbol[] { DNATools.getDNA().getGapSymbol() }, 
                             1, DNATools.getDNA());

    // array indices for subpaths, etc.
    private static final int A = 0, C = 1, G = 2, T = 3;
    
    private Chromatogram chromat;
    private float vertScale, horizScale;
    private int width, height;

  /** Flag for subpaths. */
    protected boolean subpathsValid;

  /** Flag for call boxes. */
  protected boolean callboxesValid;

  /** Flag for drawable call boxes. */
  protected boolean drawableCallboxesValid;
    
    private GeneralPath[][] subpaths;
    private Rectangle2D.Float[] callboxes;
    private Rectangle2D[] drawableCallboxes;
    private AffineTransform drawnCallboxesTx;

    /** The map containing the {@link Option}s and values for this instance. */
    protected Map options;
    /** The map containing the colors for drawing traces.  Keys are DNA Symbols. */
    protected Map colors;
    /** The map containing the fill colors for callboxes.  Keys are DNA Symbols. */
    protected Map fillColors;
    
    /** The default values for the {@link #colors} map. */
    private static final Map DEFAULT_COLORS = new HashMap();
    static {
        DEFAULT_COLORS.put(DNATools.a(),   Color.green);
        DEFAULT_COLORS.put(DNATools.c(),   Color.blue);
        DEFAULT_COLORS.put(DNATools.g(),   Color.black);
        DEFAULT_COLORS.put(DNATools.t(),   Color.red);
    }
    
    /** Default constructor with no Chromatogram. */
    public ChromatogramGraphic() {
        this(null);
    }
    
    /**
     * Creates a new <code>ChromatogramGraphic</code>, initially displaying
     * the given chromatogram.
     *
     * @param c the Chromomatogram to display
     */
    public ChromatogramGraphic(Chromatogram c) {
        options = new HashMap(Option.DEFAULTS);
        colors = new HashMap();
        fillColors = new HashMap();
        height = -1; width = -1;
        vertScale = -1.0f; horizScale = -1.0f;
        subpaths = new GeneralPath[4][];
        for (Iterator it = DEFAULT_COLORS.keySet().iterator() ; it.hasNext() ; ) {
            Symbol key = (Symbol) it.next();
            setBaseColor(key, (Color) DEFAULT_COLORS.get(key));
        }
        setChromatogram(c);
        //if (System.getProperty("os.name").equalsIgnoreCase("mac os x")) {
        //    setOption(Option.USE_PER_SHAPE_TRANSFORM, Boolean.TRUE);
        //}
    }

    /**
     * Precomputes the {@link java.awt.geom.GeneralPath}s used to draw the
     * traces.
     */
    protected synchronized void generateSubpaths() {
        if (chromat != null && !subpathsValid) {
            if (subpaths == null)
                subpaths = new GeneralPath[4][];
            int[][] samples = new int[4][];
            try {
                samples[A] = chromat.getTrace(DNATools.a());
                samples[C] = chromat.getTrace(DNATools.c());
                samples[G] = chromat.getTrace(DNATools.g());
                samples[T] = chromat.getTrace(DNATools.t());
            } catch (IllegalSymbolException ise) {
                throw new BioError("Can't happen");
            }
            ChromatogramNonlinearScaler horizScaler = (ChromatogramNonlinearScaler) getOption(Option.HORIZONTAL_NONLINEAR_SCALER);
            float max = chromat.getMax();
            int subpathLength = getIntOption(Option.SUBPATH_LENGTH);
            int countSubpaths = (int)Math.ceil( horizScaler.scale(chromat, samples[A].length-1) / subpathLength );
            int[] subpathTraceBounds = new int[countSubpaths+1];
            {
                subpathTraceBounds[0] = 0;
                int bound = 1;
                float scaled;
                for (int i = 0 ; i < samples[A].length ; i++) {
                    scaled = horizScaler.scale(chromat, i);
                    if (scaled >= bound*subpathLength) {
                        subpathTraceBounds[bound] = i;
                        bound++;
                    }
                }
                subpathTraceBounds[countSubpaths] = samples[A].length - 1;
            }
            int offsetIdx, thisPathLen;
            for (int i = 0 ; i < 4 ; i++) {
                subpaths[i] = new GeneralPath[countSubpaths];
                for (int j = 0 ; j < countSubpaths ; j++) {
                    thisPathLen = subpathTraceBounds[j+1] - subpathTraceBounds[j];
                    subpaths[i][j] = new GeneralPath(GeneralPath.WIND_EVEN_ODD, thisPathLen);
                    offsetIdx = subpathTraceBounds[j];
                    subpaths[i][j].moveTo(horizScaler.scale(chromat, offsetIdx), max - samples[i][offsetIdx]);
                    for (int k = offsetIdx+1 ; k < samples[i].length && k <= offsetIdx + thisPathLen ; k++) {
                        subpaths[i][j].lineTo(horizScaler.scale(chromat, k), max - samples[i][k]);
                    }
                }
            }
            subpathsValid = true;
        }
    }
    
    /**
     * Precomputes the {@link java.awt.geom.Rectangle2D}s that are the in-memory
     * representation of the callboxes.  These rectangles are used for drawing
     * (via generateDrawableCallboxes) as well as queries (e.g., 
     * {@link #getCallContaining}).
     */
    protected synchronized void generateCallboxes() {
        if (chromat != null && !callboxesValid) {
            if (chromat.getSequenceLength() < 2) {
                callboxes = new Rectangle2D.Float[1];
                callboxes[0] = new Rectangle2D.Float(0, 0, chromat.getTraceLength() - 1, chromat.getMax());
            }
            else {
                int[] bcOffsets = ChromatogramTools.getTraceOffsetArray(chromat);
                if (callboxes == null || callboxes.length != bcOffsets.length)
                    callboxes = new Rectangle2D.Float[bcOffsets.length];
                ChromatogramNonlinearScaler horizScaler = 
                    (ChromatogramNonlinearScaler) getOption(Option.HORIZONTAL_NONLINEAR_SCALER);
                float max = chromat.getMax();
                int leftIdx  = 0;
                int rightIdx = bcOffsets[0] + (int) Math.floor( ((double) (bcOffsets[1] - bcOffsets[0])) / 2) + 1;
                float left  = horizScaler.scale(chromat, leftIdx);
                float right = horizScaler.scale(chromat, rightIdx);
                callboxes[0] = new Rectangle2D.Float(left, 0, right - left, max);
                //System.out.println("[cg.gcb] cb[0] left="+leftIdx+":"+left+" center="+bcOffsets[0]+" right="+rightIdx+":"+right+" box="+callboxes[0]);
                int i = 1;
                while (i < bcOffsets.length - 1) {
                    leftIdx = rightIdx;
                    rightIdx = bcOffsets[i] + (int) Math.floor( ((double) (bcOffsets[i+1] - bcOffsets[i])) / 2) + 1;
                    left = right;
                    right = horizScaler.scale(chromat, rightIdx);
                    callboxes[i] = new Rectangle2D.Float(left, 0, right - left, max);
                    //System.out.println("[cg.gcb] cb["+i+"] left="+leftIdx+":"+left+" center="+bcOffsets[i]+" right="+rightIdx+":"+right+" box="+callboxes[i]);
                    i++;
                }
                left = right;
                right = horizScaler.scale(chromat, chromat.getTraceLength() - 1);
                callboxes[i] = new Rectangle2D.Float(left, 0, right - left, max);
                //System.out.println("[cg.gcb] cb["+i+"] left="+left+" center="+bcOffsets[i]+" right="+right+" box="+callboxes[i]);
            }
            callboxesValid = true;
            drawableCallboxesValid = false;
        }
        //System.out.println("[cg.gcb]");
        //for (int i = 0 ; i < callboxes.length ; i++)
        //    System.out.println("callboxes["+i+"]="+callboxes[i]);
    }
    
    /**
     * Precomputes the callboxes in screen coordinates.
     * @param shapeTx the transform to apply to the callboxes to move them into
     *        screen space.
     */
    protected synchronized void generateDrawableCallboxes(AffineTransform shapeTx) {
        if (!callboxesValid)
            generateCallboxes();
        if (!shapeTx.equals(drawnCallboxesTx) || !drawableCallboxesValid) {
            //System.out.println("Regen drawableCallboxes with " + shapeTx + " condition: tx="+shapeTx.equals(drawnCallboxesTx)+" valid="+drawableCallboxesValid);
            if (drawableCallboxes == null || drawableCallboxes.length != callboxes.length)
                drawableCallboxes = new Rectangle2D[callboxes.length];
            for (int i = 0 ; i < drawableCallboxes.length ; i++)
                drawableCallboxes[i] = shapeTx.createTransformedShape(callboxes[i]).getBounds2D();
            drawnCallboxesTx = (AffineTransform) shapeTx.clone();
            drawableCallboxesValid = true;
            //System.out.println("[cg.gdcb]");
            //for (int i = 0 ; i < drawableCallboxes.length ; i++)
            //    System.out.println("cb["+i+"]="+callboxes[i]+"\ndcb["+i+"]="+drawableCallboxes[i]);
        }
    }
    
    /**
     * Accessor for the in-use chromatogram.
     * @return the chromatogram that a call to {@link #drawTo} will draw
     */
    public Chromatogram getChromatogram() { return chromat; }
    /**
     * Sets the chromatogram to draw.
     * @param c the new chromatogram
     * @see Option#WIDTH_IS_AUTHORITATIVE
     * @see Option#HEIGHT_IS_AUTHORITATIVE
     */
    public synchronized void setChromatogram(Chromatogram c) {
        this.chromat = c;
        callboxesValid = false;
        subpathsValid = false;
        // set width, height, horizScale, vertScale for new chromat (even if null)
        if (optionIsTrue(Option.WIDTH_IS_AUTHORITATIVE)) 
            setWidth(width);
        else
            setHorizontalScale(horizScale);
        if (optionIsTrue(Option.HEIGHT_IS_AUTHORITATIVE))
            setHeight(height);
        else
            setVerticalScale(vertScale);
        // drawing bounds default to show the whole chromat
        setOption(Option.FROM_TRACE_SAMPLE, new Integer(0));
        if (c == null)
            setOption(Option.TO_TRACE_SAMPLE, new Integer(Integer.MAX_VALUE));
        else
            setOption(Option.TO_TRACE_SAMPLE, new Integer(c.getTraceLength() - 1));
        //System.out.println("chromatgfx[w=" + width + "; h=" + height + "; hs=" + horizScale + "; vs=" + vertScale + "]");
    }
    
    /**
     * Returns the width of the whole graphic (in pixels).
     *
     * @return the width
     */
    public int getWidth()  { return width;  }

    /**
     * Returns the height of the whole graphic (in pixels).
     *
     * @return the height
     */
    public int getHeight() { return height; }

    /**
     * Returns the in-use horizontal scale factor.
     * The "units" of this value are (trace samples) / pixel.
     * For example, a horizontal scale of 1.0 means that there will be one 
     * pixel horizontally for each trace sample.
     *
     * @return the horizontal scale
     */
    public float getHorizontalScale() { return horizScale; }

    /**
     * Returns the in use vertical scale factor.  
     * The "units" of this value are (trace value bins) / pixel.
     * For example, a vertical scale of 1.0 means that there will be one 
     * pixel vertically for each value in the range 
     * [0, <code>getChromatogram().getMax()</code>].
     *
     * @return the vertical scale
     */
    public float getVerticalScale()   { return vertScale;  }
    
    /** 
     * Returns the width of the graphic as it will be rendered.
     * This means that the {@link Option#FROM_TRACE_SAMPLE} and 
     * {@link Option#TO_TRACE_SAMPLE} bounds are taken into account.
     *
     * @return the rendered width
     */
    public int getRenderedWidth() { 
        return getRenderedWidth(horizScale);
    }
    
    /** 
     * Returns the width of the graphic as it would be rendered with
     * the specified horizontal scale. The {@link Option#FROM_TRACE_SAMPLE} and 
     * {@link Option#TO_TRACE_SAMPLE} bounds are taken into account.
     *
     * @param horizontalScale  the horizontal scale
     * @return the rendered width at that scale
     */
    public int getRenderedWidth(float horizontalScale) {
        return (int) Math.ceil(
                        horizontalScale *
                            (getFloatOption(Option.TO_TRACE_SAMPLE)   -
                             getFloatOption(Option.FROM_TRACE_SAMPLE) + 
                             1)
                     );
    }
    
    /**
     * Sets the height (in pixels).  This will also change the
     * vertical scale.
     *
     * @param h the desired height in pixels
     * @see Option#HEIGHT_IS_AUTHORITATIVE
     */
    public void setHeight(int h) {
        if (h != height) drawableCallboxesValid = false;
        height = h;
        if (chromat != null) 
            vertScale = ( (float) height ) / chromat.getMax();
        else
            vertScale = -1.0f;
    }
    
    /**
     * Sets the vertical scale (proportional).  This will also
     * change the height.
     *
     * @param vs the desired vertical scale.  See {@link #getVerticalScale}
     *        for semantics.
     * @see Option#HEIGHT_IS_AUTHORITATIVE
     */
    public void setVerticalScale(float vs) {
        if (vs != vertScale) drawableCallboxesValid = false;
        vertScale = vs;
        if (chromat != null)
            height = (int) (vertScale * chromat.getMax());
        else
            height = -1;
    }
    
    /**
     * Sets the width of the whole graphic (in pixels).  This will also change 
     * the horizontal scale.
     *
     * @param w the desired width in pixels
     * @see Option#WIDTH_IS_AUTHORITATIVE
     */
    public void setWidth(int w) {
        if (w != width) drawableCallboxesValid = false;
        width = w;
        ChromatogramNonlinearScaler horizScaler = 
            (ChromatogramNonlinearScaler) getOption(Option.HORIZONTAL_NONLINEAR_SCALER);
        if (chromat != null)
            horizScale = ( (float) width ) / horizScaler.scale(chromat, chromat.getTraceLength()-1);
        else
            horizScale = -1.0f;
    }
    
    /**
     * Sets the horizontal scale (proportional).  This will also
     * change the width.
     *
     * @param hs the desired vertical scale.  See {@link #getHorizontalScale}
     *        for semantics.
     * @see Option#WIDTH_IS_AUTHORITATIVE
     */
    public void setHorizontalScale(float hs) {
        if (hs != horizScale) drawableCallboxesValid = false;
        horizScale = hs;
        
        ChromatogramNonlinearScaler horizScaler = 
            (ChromatogramNonlinearScaler) getOption(Option.HORIZONTAL_NONLINEAR_SCALER);
        if (chromat != null)
            width = (int) (horizScale * horizScaler.scale(chromat, chromat.getTraceLength()-1));
        else
            width = -1;
    }

    /** 
     * Returns the color that will be used to draw the trace for the
     * given DNA symbol.
     * @param b the symbol
     * @return the color, or null if none is set
     */
    public Color getBaseColor(Symbol b)     { return (Color) colors.get(b); }
    /** 
     * Returns the color that will be used to fill in the callboxes for
     * calls with the given symbol.
     * @param b the symbol
     * @return the color, or null if none is set
     */
    public Color getBaseFillColor(Symbol b) { return (Color) fillColors.get(b); }
    /**
     * Maps a color to a DNA symbol.  The color as specified will be used for
     * to draw the trace for the symbol (if any).  The fill color for calls to
     * the symbol will be derived from the trace color.
     * @param b the symbol
     * @param c the color
     */
    public void setBaseColor(Symbol b, Color c) {
        colors.put(b, c);
        // fade color
        float[] hsb = Color.RGBtoHSB(c.getRed(), c.getGreen(), c.getBlue(), null);
        //System.out.println("Base: "+b+"; H="+hsb[0]+"; S="+hsb[1]+"; B="+hsb[2]);
        fillColors.put(b, Color.getHSBColor(hsb[0], hsb[1] * 0.09f, Math.max(hsb[2], 0.8f)));
    }
    
    /**
     * Returns the number of callboxes, regenerating them if necessary.  This 
     * should almost always equal 
     * <code>getChromatogram().getSequenceLength()</code>
     * @return the number of callboxes
     */
    public int getCallboxCount() {
        if (!callboxesValid)
            generateCallboxes();
        if (callboxesValid)
            return callboxes.length;
        else
            return 0;
    }
    
    /**
     * Returns the screen-coordinate bounds of the callbox for a given call.
     *
     * @param index the callbox for which to get the bounds <b>0-based</b>
     * @return a Rectangle2D giving the bounds of the call box
     */
    public Rectangle2D getCallboxBounds(int index) {
        return getCallboxBounds(index, true);
    }
    
    /**
     * Returns the bounds of the callbox for a given call.
     * @param index the callbox for which to get the bounds <b>0-based</b>.
     * @param boundsOnScreen determines the coordinate system of the returned
     *        bounds
     * @return the bounds of the callbox in screen coordinates if 
     *          <code>boundsOnScreen</code> is true, otherwise the bounds
     *          of the callbox in chromatogram coordinates
     */
    public Rectangle2D getCallboxBounds(int index, boolean boundsOnScreen) {
        if (chromat != null && index >= 0 && index < getCallboxCount()) {
            if (!callboxesValid)
                generateCallboxes();
            if (boundsOnScreen) {
                if (!drawableCallboxesValid)
                    generateDrawableCallboxes(getTransform());
                return drawableCallboxes[index].getBounds2D();
                //return getTransform().createTransformedShape(callboxes[index]).getBounds2D();
            }
            else {
                return callboxes[index].getBounds2D();
            }
        }
        else {
            return null;
        }
    }
    
    /**
     * Returns the <b>0-based</b> index of the call containing a given
     * point.  The point may be either in screen space or chromatogram
     * space, scale-wise.  If the point is in screen space, the caller must 
     * translate the point such that if it is, for instance, from a mouse 
     * click, a click on the upper-left corner of the graphic would be (0,0).
     * @param point the point to search for
     * @param pointOnScreen if true, the point will be treated as though it
     *        is in screen space.  Otherwise, it will be considered to be
     *        in chromatogram space.
     * @return the <b>0-based</b> index of the callbox which contains the point
     */
    public int getCallContaining(Point2D point, boolean pointOnScreen) {
        if (chromat != null) {
            if (!callboxesValid)
                generateCallboxes();
            Point2D trans = new Point2D.Double(point.getX(), point.getY());
            if (pointOnScreen)
                getInvTransform().transform(point, trans);
            int i = 0;
            // XXX: binary search is possible since callboxes is sorted
            while (i < callboxes.length && trans.getX() > callboxes[i].getMaxX())
                i++;
            return (i < callboxes.length) ? (i) : (i-1);
        }
        else {
            return 0;
        }
    }
    
    /** 
     * Synonym for {@link #getCallContaining(Point2D, boolean)} with
     * <code>pointOnScreen</code>=true.
     *
     * @param point  the Point2D to search with
     * @return the call containing this point
     */
    public int getCallContaining(Point2D point) {
        return getCallContaining(point, true);
    }
    
    /**
     * Same as {@link #getCallContaining(Point2D, boolean)}, except that
     * only the x-coordinate of the point is specified.
     *
     * @param x the x-coordinate to search for
     * @param xOnScreen whether the coordinate in screen space or chromatogram
     *        space
     * @return the index of the call containing the position x
     */
    public int getCallContaining(float x, boolean xOnScreen) {
        return getCallContaining(new Point2D.Float(x, 0), xOnScreen);
    }
    
    /** 
     * Synonym for {@link #getCallContaining(float, boolean)} with
     * <code>pointOnScreen</code>=true.
     *
     * @param x the x-coordinate to search for
     * @return the index of the call containing the position x
     */
    public int getCallContaining(float x) {
        return getCallContaining(x, true);
    }
    
    private int getSubpathContaining(float x) {
        // XXX: this could be a binary search
        Rectangle2D bounds;
        if (x < subpaths[A][0].getBounds2D().getX())
            return 0;
        for (int i = 0 ; i < subpaths[A].length ; i++) {
            bounds = subpaths[A][i].getBounds2D();
            if (bounds.getX() <= x && x <= bounds.getMaxX())
                return i;
        }
        return subpaths[A].length - 1;
    }
    
    /**
     * Returns a new AffineTransform describing the transformation
     * from chromatogram coordinates to output coordinates.
     */
    public AffineTransform getTransform() {
        AffineTransform at = new AffineTransform();
        getTransformAndConcat(at);
        return at;
    }
    
    /**
     * Concatenates the chromatogram-to-output transform to the 
     * provided given AffineTransform.
     */
    public void getTransformAndConcat(AffineTransform target) {
        target.scale(horizScale, vertScale);
        target.translate(-1.0 * getFloatOption(Option.FROM_TRACE_SAMPLE), 0.0);
    }
    
    /** 
     * Returns a new AffineTransform describing the transformation from
     * output space to chromatogram space.  Should be much more efficient
     * than <code>getTransform().createInverse()</code>
     */
    public AffineTransform getInvTransform() {
        AffineTransform at = new AffineTransform();
        at.translate(getFloatOption(Option.FROM_TRACE_SAMPLE), 0.0);
        at.scale(1.0 / horizScale, 1.0 / vertScale);
        return at;
    }
    
    /**
     * Draws the chromatogram onto the provided graphics context.
     *
     * @param g2  the Graphics2D to draw to
     */
    public void drawTo(Graphics2D g2) {
        //System.out.println("drawTo(" + g2 + ", " + fromTraceSample + ", " + toTraceSample + ")");
        AffineTransform origTx = g2.getTransform();
        //System.out.println("origTx:   " + origTx);
        Color origC = g2.getColor();
        Shape origClip = g2.getClip();
        //System.out.println("origClip: " + origClip);
        Stroke origStroke = g2.getStroke();
        
        Rectangle2D clip;
        if (origClip != null)
            clip = origClip.getBounds2D();
        else
            clip = new Rectangle(0, 0, getWidth(), getHeight());
        //System.out.println("clipping bounds: " + clip);

        // Allow user to decide whether AA should be on
        //g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        //g2.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_NORMALIZE);

        boolean usePerShpTx = optionIsTrue(Option.USE_PER_SHAPE_TRANSFORM);
        
        AffineTransform shapeTx = (AffineTransform) origTx.clone();
        getTransformAndConcat(shapeTx);
        if (usePerShpTx)
            g2.setTransform(new AffineTransform()); // set to identity, do tx per shape
        else
            g2.setTransform(shapeTx); // do tx in graphics2d object
        //System.out.println("using g2.Tx="+g2.getTransform());
        //System.out.println("shapeTx="+shapeTx);
        
        if (optionIsTrue(Option.USE_CUSTOM_STROKE))
            g2.setStroke((Stroke) getOption(Option.SEPARATOR_STROKE));
        if (!callboxesValid)
            generateCallboxes();
        if (usePerShpTx) {
            generateDrawableCallboxes(AffineTransform.getScaleInstance(shapeTx.getScaleX(), shapeTx.getScaleY()));
            g2.translate(shapeTx.getTranslateX(), shapeTx.getTranslateY());
        }
        int leftCbIdx = getCallContaining((float)clip.getX()); 
        int rightCbIdx = getCallContaining((float)clip.getMaxX());
        //System.out.println("Drawing calls with stroke=" + basicStrokeToString((BasicStroke)g2.getStroke()));
        //System.out.println("Drawing " + (rightCbIdx - leftCbIdx + 1) + " visible calls from " + leftCbIdx + " to " + rightCbIdx);
        SymbolList calls = ChromatogramTools.getDNASequence(chromat);
        // if the chromatogram has no calls, fake it
        if (calls.length() < 1)
            calls = SINGLE_CALL;
        boolean drawSep = optionIsTrue(Option.DRAW_CALL_SEPARATORS);
        Color fill, line;
        Line2D.Double sep = new Line2D.Double();
        int callIdx;
        for (int i = leftCbIdx ; i <= rightCbIdx ; i++) {
            callIdx = i + 1;
            if (doDrawCallbox(calls.symbolAt(callIdx))) {
                fill = (Color) fillColors.get(calls.symbolAt(callIdx));
                line = (Color) colors.get(calls.symbolAt(callIdx));
                if (line == null) {
                    line = Color.black;
                    fill = Color.white;
                }
                g2.setColor(fill);
                if (usePerShpTx) {
                    g2.fill(drawableCallboxes[i]);
                    //System.out.println("Drawing drawableCallboxes["+i+"]="+drawableCallboxes[i]);
                }
                else {
                    g2.fill(callboxes[i]);
                    //System.out.println("Drawing callboxes["+i+"]="+callboxes[i]);
                }
            }
            if (drawSep) {
                g2.setColor((Color) getOption(Option.SEPARATOR_COLOR));
                if (usePerShpTx)
                    sep.setLine(drawableCallboxes[i].getX(), 
                                drawableCallboxes[i].getY(), 
                                drawableCallboxes[i].getX(), 
                                drawableCallboxes[i].getMaxY());
                else
                    sep.setLine(callboxes[i].x, callboxes[i].y, callboxes[i].x, callboxes[i].y + callboxes[i].height);
                //System.out.println("sep["+i+"]=(" + sep.x1 + ", " + sep.y1 + ") -> (" + sep.x2 + ", " + sep.y2 + ")");
                //System.out.println("        " + sep.getBounds2D());
                g2.draw(sep);
            }
        }
        if (usePerShpTx)
            g2.translate(-1 * shapeTx.getTranslateX(), -1 * shapeTx.getTranslateY());

        
        if (optionIsTrue(Option.USE_CUSTOM_STROKE))
            g2.setStroke((Stroke) getOption(Option.TRACE_STROKE));
        if (!subpathsValid)
            generateSubpaths();
        float toTraceSample   = getFloatOption(Option.TO_TRACE_SAMPLE  );
        float fromTraceSample = getFloatOption(Option.FROM_TRACE_SAMPLE);
        //System.out.println("Drawing from trace sample " + fromTraceSample + " to " + toTraceSample);
        //System.out.println("Drawing traces with stroke=" + basicStrokeToString((BasicStroke)g2.getStroke()));
        int loSubpath, hiSubpath;
        //loSubpath = (int)Math.floor((clip.getX() / shapeTx.getScaleX() + fromTraceSample) / subpathLength);
        loSubpath = getSubpathContaining((float) (clip.getX() / shapeTx.getScaleX() + fromTraceSample));
        //hiSubpath = (int)Math.floor(Math.min(clip.getMaxX() / shapeTx.getScaleX() + fromTraceSample, toTraceSample) / subpathLength);
        hiSubpath = getSubpathContaining((float) Math.min(clip.getMaxX() / shapeTx.getScaleX() + fromTraceSample, toTraceSample));
        //System.out.println("Drawing subpaths ["+loSubpath+","+hiSubpath+"] of "+subpaths[0].length+" each "+subpathLength+" long");
        if (optionIsTrue(Option.DRAW_TRACE_A)) {
            g2.setColor((Color) colors.get(DNATools.a()));
            if (usePerShpTx)
                for (int j = loSubpath ; j <= hiSubpath ; j++)
                    g2.draw(shapeTx.createTransformedShape(subpaths[A][j]));
            else
                for (int j = loSubpath ; j <= hiSubpath ; j++)
                    g2.draw(subpaths[A][j]);
        }
        if (optionIsTrue(Option.DRAW_TRACE_C)) {
            g2.setColor((Color) colors.get(DNATools.c()));
            if (usePerShpTx)
                for (int j = loSubpath ; j <= hiSubpath ; j++)
                    g2.draw(shapeTx.createTransformedShape(subpaths[C][j]));
            else
                for (int j = loSubpath ; j <= hiSubpath ; j++)
                    g2.draw(subpaths[C][j]);
        }
        if (optionIsTrue(Option.DRAW_TRACE_G)) {
            g2.setColor((Color) colors.get(DNATools.g()));
            if (usePerShpTx)
                for (int j = loSubpath ; j <= hiSubpath ; j++)
                    g2.draw(shapeTx.createTransformedShape(subpaths[G][j]));
            else
                for (int j = loSubpath ; j <= hiSubpath ; j++)
                    g2.draw(subpaths[G][j]);
        }
        if (optionIsTrue(Option.DRAW_TRACE_T)) {
            g2.setColor((Color) colors.get(DNATools.t()));
            if (usePerShpTx)
                for (int j = loSubpath ; j <= hiSubpath ; j++)
                    g2.draw(shapeTx.createTransformedShape(subpaths[T][j]));
            else
                for (int j = loSubpath ; j <= hiSubpath ; j++)
                    g2.draw(subpaths[T][j]);
        }
        
        g2.setStroke(origStroke);
        g2.setTransform(origTx);
        g2.setColor(origC);
        g2.setClip(origClip);
    }

    /**
     * Sets a new value for the specified option.  Be sure that the
     * value is appropriate per the documentation, or you'll induce a 
     * ClassCastException somewhere else.
     *
     * @param opt Option to set
     * @param value new value for the option
     * @see Option
     */
    public void setOption(Option opt, Object value) {
        options.put(opt, value);
        if (opt == Option.SUBPATH_LENGTH) {
            subpathsValid = false;
        }
        if (opt == Option.HORIZONTAL_NONLINEAR_SCALER) {
            subpathsValid = false;
            callboxesValid = false;
            drawableCallboxesValid = false;
            if (optionIsTrue(Option.WIDTH_IS_AUTHORITATIVE))
                setWidth(width);
            else
                setHorizontalScale(horizScale);
        }
    }
    
    /**
     * Returns the current value for the specified option.
     * @see Option
     */
    public Object getOption(Option opt) {
        return options.get(opt);
    }
    
    /**
     * Helper method for converting a {@link java.lang.Boolean}-valued
     * option into a <code>boolean</code> primitive.
     *
     * @param opt the {@link Option} to convert
     * @return true if the option is enabled
     * @throws ClassCastException when the option isn't <code>Boolean</code>-valued
     */
    public boolean optionIsTrue(Option opt) throws ClassCastException {
        if (getOption(opt) instanceof Boolean)
            return ((Boolean)getOption(opt)).booleanValue();
        else
            throw new ClassCastException("Option \""+opt+"\" is not set to a Boolean value");
    }
    
    /**
     * Helper method for converting a {@link java.lang.Number}-valued
     * option into a <code>float</code> primitive.
     * @param opt the {@link Option} to convert
     * @throws ClassCastException when the option isn't <code>Number</code>-valued
     */
    public float getFloatOption(Option opt) throws ClassCastException {
        if (getOption(opt) instanceof Number)
            return ((Number) getOption(opt)).floatValue();
        else
            throw new ClassCastException("Option \""+opt+"\" is not set to a Number value");
    }
    
    /**
     * Helper method for converting a {@link java.lang.Number}-valued
     * option into an <code>int</code> primitive.
     * @param opt the {@link Option} to convert
     * @throws ClassCastException when the option isn't <code>Number</code>-valued
     */
    public int getIntOption(Option opt) throws ClassCastException {
        if (getOption(opt) instanceof Number)
            return ((Number) getOption(opt)).intValue();
        else
            throw new ClassCastException("Option \""+opt+"\" is not set to a Number value");
    }
    
    /**
     * Utility method for determining whether to draw a callbox for a particular called Symbol.
     *
     * @param bc Symbol to evaluate
     * @return true if this should be drawn
     */
    private boolean doDrawCallbox(Symbol bc) {
        if      (bc == DNATools.a()) return optionIsTrue(Option.DRAW_CALL_A);
        else if (bc == DNATools.c()) return optionIsTrue(Option.DRAW_CALL_C);
        else if (bc == DNATools.g()) return optionIsTrue(Option.DRAW_CALL_G);
        else if (bc == DNATools.t()) return optionIsTrue(Option.DRAW_CALL_T);
        else if (DNATools.getDNA().contains(bc)) return optionIsTrue(Option.DRAW_CALL_OTHER);
        else return false;
    }
    
    /**
     * Performs a partial deep copy and invalidates regenerable structures.
     *
     * @return an Object that is castable to ChromatogramGraphic
     */
    public Object clone() {
        ChromatogramGraphic copy = null;
        try {
            copy = (ChromatogramGraphic) super.clone();
            copy.callboxesValid = false;
            copy.drawableCallboxesValid = false;
            copy.drawableCallboxes = null;
            copy.subpathsValid = false;
            copy.subpaths = null;
            // copy options
            copy.options = new HashMap();
            for (Iterator it = this.options.keySet().iterator() ; it.hasNext() ; ) {
                Object next = it.next();
                copy.options.put(next, this.options.get(next));
            }
        } catch (CloneNotSupportedException e) {
            System.err.println(e);
            throw new BioError("Can't happen");
        }
        return copy;
    }
    
    /**
     * A typesafe enumeration of the options available for configuring
     * the behavior of a {@link ChromatogramGraphic} instance.
     * The semantics and expected values are described with the
     * enumerated options.
     *
     * @author Rhett Sutphin (<a href="http://genome.uiowa.edu/">UI CBCB</a>)
     * @since 1.3
     */
    public static class Option {
        private String desc;
        private static Map map = new HashMap();

      /**
       * Create a new Obtion.
       *
       * @param desc  option description
       * @param def   option default
       */
        private Option(String desc, Object def) {
            this.desc = desc;
            map.put(desc, this);
            DEFAULTS.put(this, def);
        }
        
        public String toString() { return desc; }
        
        /**
         * Looks up an <code>Option</code> instance based on its
         * string description.
         * @param desc the description of the desired <code>Option</code>
         * @return the <code>Option</code> with the specified description
         *          or null if there isn't one
         */
        public static final Option lookup(String desc) {
            return (Option) map.get(desc);
        }
        
        /**
         * Default values table
         */
        static final Map DEFAULTS = new HashMap();
        
        /**
         * Option indicating whether to fill in the callboxes for calls of
         * nucleotide A.  
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.TRUE</code>.
         * </p>
         */
        public static final Option DRAW_CALL_A   = new Option("draw-A-calls",     Boolean.TRUE);
        /**
         * Option indicating whether to fill in the callboxes for calls of
         * nucleotide C.  
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.TRUE</code>.
         * </p>
         */
        public static final Option DRAW_CALL_C   = new Option("draw-C-calls",     Boolean.TRUE);
        /**
         * Option indicating whether to fill in the callboxes for calls of
         * nucleotide G.  
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.TRUE</code>.
         * </p>
         */
        public static final Option DRAW_CALL_G   = new Option("draw-G-calls",     Boolean.TRUE);
        /**
         * Option indicating whether to fill in the callboxes for calls of
         * nucleotide T.  
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.TRUE</code>.
         * </p>
         */
        public static final Option DRAW_CALL_T   = new Option("draw-T-calls",     Boolean.TRUE);
        /**
         * Option indicating whether to fill in the callboxes for non-base calls
         * (gaps, ambiguities).
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.TRUE</code>.
         * </p>
         */
        public static final Option DRAW_CALL_OTHER = new Option("draw-other-calls", Boolean.TRUE);
        
        /**
         * Option indicating whether to draw the chromatogram trace for 
         * nucleotide A.
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.TRUE</code>.
         * </p>
         */
        public static final Option DRAW_TRACE_A = new Option("draw-A-trace", Boolean.TRUE);
        /**
         * Option indicating whether to draw the chromatogram trace for 
         * nucleotide C.
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.TRUE</code>.
         * </p>
         */
        public static final Option DRAW_TRACE_C = new Option("draw-C-trace", Boolean.TRUE);
        /**
         * Option indicating whether to draw the chromatogram trace for 
         * nucleotide G.
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.TRUE</code>.
         * </p>
         */
        public static final Option DRAW_TRACE_G = new Option("draw-G-trace", Boolean.TRUE);
        /**
         * Option indicating whether to draw the chromatogram trace for 
         * nucleotide T.
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.TRUE</code>.
         * </p>
         */
        public static final Option DRAW_TRACE_T = new Option("draw-T-trace", Boolean.TRUE);
        
        /**
         * Option indicating whether to draw vertical lines separating
         * the calls.
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.TRUE</code>.
         * </p>
         */
        public static final Option DRAW_CALL_SEPARATORS = 
            new Option("draw-call-separators", Boolean.TRUE);
        /**
         * Option indicating the color that the call separators
         * should be.
         * <p>
         * Value type: {@link java.awt.Color}.<br/>
         * Default value: <code>Color.lightGray</code>.
         * </p>
         */
        public static final Option SEPARATOR_COLOR = 
            new Option("separator-color", Color.lightGray);
        
        /**
         * Option indicating whether width or horizontal scale is
         * the authoritative measure.  If the value is true, then
         * when the Chromatogram displayed by the graphic is changed, the 
         * horizontal scale may be changed but the width will stay the same.
         * If the value is false, the width may change but the horizontal
         * scale will stay the same.
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.FALSE</code>.
         * </p>
         */
        public static final Option WIDTH_IS_AUTHORITATIVE  = 
            new Option("width-is-authoritative", Boolean.FALSE);
        /**
         * Option indicating whether height or vertical scale is
         * the authoritative measure.  If the value is true, then
         * when the Chromatogram displayed by the graphic is changed, the 
         * vertical scale may be changed but the height will stay the same.
         * If the value is false, the height may change but the vertical
         * scale will stay the same.
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.TRUE</code>.
         * </p>
         */
        public static final Option HEIGHT_IS_AUTHORITATIVE = 
            new Option("height-is-authoritative", Boolean.TRUE);
        
        /**
         * Option indicating whether to use custom strokes when
         * drawing traces and separators.
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.TRUE</code>.
         * </p>
         */
        public static final Option USE_CUSTOM_STROKE = 
            new Option("use-custom-stroke", Boolean.TRUE);
        /**
         * Option providing the the stroke to use for drawing
         * the chromatogram traces.
         * <p>
         * Value type: {@link java.awt.Stroke}.<br/>
         * Default value: {@link BasicStroke} with width 1.0, cap CAP_ROUND, join JOIN_ROUND.
         * </p>
         */
        public static final Option TRACE_STROKE = 
            new Option("trace-stroke", new BasicStroke(1.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        /**
         * Option providing the the stroke to use for drawing
         * call separators.
         * <p>
         * Value type: {@link java.awt.Stroke}.<br/>
         * Default value: {@link BasicStroke} with width 1.0, default cap & join.
         * </p>
         */
        public static final Option SEPARATOR_STROKE = 
            new Option("separator-stroke", new BasicStroke(1.0f));

        /**
         * Option indicating whether to apply scaling and translation
         * transforms to each shape individually or to apply a single
         * transform to the graphics context.  For putative performance
         * reasons, the latter is the default.  However, setting this 
         * property to true may result in more attractive output, particularly
         * when the horizontal and vertical scales are very different.
         * This value must also be set to true if using a custom stroke
         * while drawing into a Swing graphics context on JDK 1.3.1 on Mac OS X,
         * due to a nasty rendering bug on that platform.
         * <p>
         * Value type: {@link java.lang.Boolean}.<br/>
         * Default value: <code>Boolean.FALSE</code>.
         * </p>
         */
        public static final Option USE_PER_SHAPE_TRANSFORM =
            new Option("use-per-shape-transform", Boolean.FALSE);
        
        /**
         * To improve performance, the drawing objects for the chromatogram
         * traces are precomputed.  Specifically, the traces are stored as a set
         * of {@link java.awt.geom.GeneralPath}s. This option indicates how 
         * long (in trace samples) each one of these should be.  Ideally, this
         * value would be slightly more than the average number of trace samples
         * visible at once in the application using the graphic.  However,
         * constantly changing this value is counterproductive as it forces the
         * recalculation of the subpaths.  In general, having a value that is
         * too small should be preferred to one that is too large.
         * <p>
         * Value type: {@link java.lang.Integer}.<br/>
         * Default value: <code>250</code>.
         * </p>
         */
        public static final Option SUBPATH_LENGTH =
            new Option("subpath-length", new Integer(250));
        
        /**
         * Option indicating the lowest (leftmost) trace sample that should be
         * drawn.  The sample at this (0-based) index will be drawn at x=0 in the
         * output chromatogram.  Note that this option is reset to zero every time
         * {@link ChromatogramGraphic#setChromatogram} is called.
         * <p>
         * Value type: {@link java.lang.Integer}.<br/>
         * Default value: <code>0</code>.
         * </p>
         */
        public static final Option FROM_TRACE_SAMPLE = 
            new Option("from-trace-sample", new Integer(0));
        
        /**
         * Option indicating the highest (rightmost) trace sample that should be
         * drawn.  The sample at this (0-based) index will be the last drawn in the
         * output chromatogram.  Note that this option is reset to the length of the new
         * chromatogram every time {@link ChromatogramGraphic#setChromatogram} is called.
         * <p>
         * Value type: {@link java.lang.Integer}.<br/>
         * Default value: <code>Integer.MAX_VALUE</code>.
         * </p>
         */
        public static final Option TO_TRACE_SAMPLE =
            new Option("to-trace-sample",   new Integer(Integer.MAX_VALUE));
        
        /**
         * Option specifying the non-linear scaling function to apply, as
         * embodied in a {@link ChromatogramNonlinearScaler} object.
         * <p>
         * Value type: {@link ChromatogramNonlinearScaler}.<br/>
         * Default value: an instance of {@link ChromatogramNonlinearScaler.Identity}.
         */
        public static final Option HORIZONTAL_NONLINEAR_SCALER =
            new Option("horiz-nonlinear-scaler", ChromatogramNonlinearScaler.Identity.getInstance());
        
        //public static final Option EXAGGERATE_SAMPLE_POINTS =
        //    new Option("exaggerate-sample-points", Boolean.FALSE);
    }
    
    /*
    private final static String basicStrokeToString(BasicStroke bs) {
        StringBuffer sb = new StringBuffer(bs.toString());
        sb.append("[width=").append(bs.getLineWidth());
        sb.append("; EndCap=");
        switch (bs.getEndCap()) {
            case BasicStroke.CAP_BUTT:   sb.append("CAP_BUTT");   break;
            case BasicStroke.CAP_ROUND:  sb.append("CAP_ROUND");  break;
            case BasicStroke.CAP_SQUARE: sb.append("CAP_SQUARE"); break;
        }
        sb.append("; Join=");
        switch (bs.getLineJoin()) {
            case BasicStroke.JOIN_BEVEL: sb.append("JOIN_BEVEL"); break;
            case BasicStroke.JOIN_MITER: sb.append("JOIN_MITER"); break;
            case BasicStroke.JOIN_ROUND: sb.append("JOIN_ROUND"); break;
        }
        sb.append("; MiterLimit=").append(bs.getMiterLimit()).append(']');
        return sb.toString();
    }
     */       
}
