package org.biojava.bio.gui.sequence;

import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.symbol.SymbolList;


/**
 * A context providing information for rendering sequences into circular coordinate systems.
 *
 * <p>
 * <b>Note:</b> All angles are measured in radians, using the normal Java graphics concepts of
 * angles.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public interface CircularRendererContext {
  /**
   * Get the angle through which the origin of the sequence is rotated through.
   *
   * <p>
   * This is equivalent to adding the offset to all calculated angles.
   * </p>
   *
   * @return the rotation offset
   */
  public double getOffset();

  /**
   * Return the angle for an index into a sequence.
   *
   * @param indx  the sequence offset
   * @return  the angle this offset is to be rendered to
   */
  public double getAngle(int indx);

  /**
   * Calculate the position in the sequence relating to the angle.
   *
   * @param angle  the angle arround the circle
   * @return  the index of the symbol rendered at that angle
   */
  public int getIndex(double angle);

  /**
   * Get the current radius at which data should be rendered.
   *
   * @return the radius
   */
  public double getRadius();


  /**
   *  The SymbolList that is currently rendered by this context.
   *
   * @return    the Sequence value
   */
  SymbolList getSymbols();

  /**
   * The features to render.
   *
   * @return a FeatureHolder with the Features to render
   */
  FeatureHolder getFeatures();
}
