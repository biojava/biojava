package org.biojava.bio.seq.projection;

import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;

/**
 * A ProjectionContext that translates and optionaly flips features.
 *
 * <p>
 * Use this to 'reverse complement' a feature hierachy, or just to shift it
 * sideways a bit.
 * </p>
 *
 * <p>
 * If the flipping mode is dissabled, then all translated features are projected
 * as having locations equivalent to feat.getLocation().translate(translation).
 * If the flipping mode is enabled, then all features are flipped arround
 * translation so that translation-i becomes translation+i.
 * </p>
 *
 * @author Matthew Pocock
 */
public class TranslateFlipContext
extends ReparentContext {
  // Fixme: we should have a simple constructor that takes a
  // translation and a StrandedFeature.Strand and does the flip
  private final int translation;
  private final boolean oppositeStrand;

  /**
   * Create a new TranslateFlipContext with explicit translation and flip.
   *
   * <p>
   * Locations will be mapped according to the rules in @link ProjectionUtils.
   * </p>
   *
   * @param parent          the parent to graft all projected features onto
   * @param wrapped         the featurs to project
   * @param translate       the translation
   * @param oppositeStrand  wether or not to flip
   */
  public TranslateFlipContext(FeatureHolder parent,
                 FeatureHolder wrapped,
                 int translate,
                 boolean oppositeStrand)
  {
    super(parent, wrapped);
    this.translation = translate;
    this.oppositeStrand = oppositeStrand;
  }

  /**
   * Create a new TranslateFlipContext that flips all featurs arround min and
   * max.
   *
   * <p>
   * A Location at exactly min will become one at max, and a Location at exactly
   * max will become one at min.
   * </p>
   *
   * <p>
   * This is equivalent to
   * <code>TranslateFlipContext(parent, wrapped, min + max, true)</code> and is
   * provided to make client code more readable.
   * </p>
   *
   * @param parent  the parent to graft all projected features ont
   * @param wrapped the features to project
   * @param min     the lower position
   * @param max     the higher position
   */
  public TranslateFlipContext(FeatureHolder parent,
                              FeatureHolder wrapped,
                              int min,
                              int max)
  {
    super(parent, wrapped);

    if(min > max) {
      throw new IllegalArgumentException("Max must not be less than min: " +
                                         min + "," + max);
    }

    this.translation = min + max;
    this.oppositeStrand = true;
  }


  /**
   * Create a new TranslateFlipContext with translation only.
   *
   * <p>
   * This is equivalent to
   * <code>TranslateFlipContext(parent, wrapped, translation, false)</code> and
   * is provided to make client code more readable.
   * </p>
   *
   * @param parent          the parent to graft all projected features onto
   * @param wrapped         the featurs to project
   * @param translation       the translation
   */
  public TranslateFlipContext(FeatureHolder parent,
                              FeatureHolder wrapped,
                              int translation)
  {
    super(parent, wrapped);

    this.translation = translation;
    this.oppositeStrand = false;
  }

  public final int getTranslation() {
    return translation;
  }

  public final boolean isOppositeStrand() {
    return oppositeStrand;
  }

  public Location projectLocation(Location oldLoc) {
    return oldLoc.newInstance(ProjectionUtils.transformLocation(
            oldLoc, translation, oppositeStrand));
  }

  public final Location revertLocation(Location oldLoc) {
    return oldLoc.newInstance(ProjectionUtils.revertLocation(
            oldLoc, translation, oppositeStrand));
  }

  public final StrandedFeature.Strand projectStrand(StrandedFeature.Strand strand) {
      if (oppositeStrand) {
          return strand.flip();
      } else {
          return strand;
      }
  }

  public final StrandedFeature.Strand revertStrand(StrandedFeature.Strand strand) {
    return projectStrand(strand);
  }

  protected FilterUtils.FilterTransformer getTransformer() {
    final FilterUtils.FilterTransformer delegate = super.getTransformer();

    return new FilterUtils.FilterTransformer() {
      public FeatureFilter transform(FeatureFilter ff) {
        ff = delegate.transform(ff);

        if (ff instanceof FeatureFilter.OverlapsLocation) {
          return new FeatureFilter.OverlapsLocation(projectLocation(((FeatureFilter.OverlapsLocation) ff).getLocation()));
        } else if (ff instanceof FeatureFilter.ContainedByLocation) {
          return new FeatureFilter.ContainedByLocation(projectLocation(((FeatureFilter.ContainedByLocation) ff).getLocation()));
        } else if (ff instanceof FeatureFilter.StrandFilter) {
          return new FeatureFilter.StrandFilter(projectStrand(((FeatureFilter.StrandFilter) ff).getStrand()));
        } else {
          return ff;
        }
      }
    };
  }

  protected FilterUtils.FilterTransformer getReverter() {
    final FilterUtils.FilterTransformer delegate = super.getReverter();

    return new FilterUtils.FilterTransformer() {
      public FeatureFilter transform(FeatureFilter ff) {
        ff = delegate.transform(ff);

        if (ff instanceof FeatureFilter.OverlapsLocation) {
          return new FeatureFilter.OverlapsLocation(revertLocation(((FeatureFilter.OverlapsLocation) ff).getLocation()));
        } else if (ff instanceof FeatureFilter.ContainedByLocation) {
          return new FeatureFilter.ContainedByLocation(revertLocation(((FeatureFilter.ContainedByLocation) ff).getLocation()));
        } else if (ff instanceof FeatureFilter.StrandFilter) {
          return new FeatureFilter.StrandFilter(revertStrand(((FeatureFilter.StrandFilter) ff).getStrand()));
        } else {
          return ff;
        }
      }
    };
  }
}
