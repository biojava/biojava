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

package org.biojava.bio.seq.projection;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.RangeLocation;

/**
 * Some common things you want to do while projecting features.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.3
 */
public class ProjectionUtils {
  /**
   * Transform a location, translating and flipping as required.
   *
   * <p>
   * If oppositeStrand is false, this is equivalent to translating the location.
   * If it is true, this is equivalent to flipping it.
   * </p>
   *
   * @param oldLoc          the Location to transform
   * @param translation     the translation to apply
   * @param oppositeStrand  wether or not this is a flip
   * @return  the transformed location
   */
  public static Location transformLocation(Location oldLoc,
                                           int translation,
                                           boolean oppositeStrand)
  {
    if (oppositeStrand) {
      return flipLocation(oldLoc, translation);
    } else {
      return oldLoc.translate(translation);
    }
  }

  /**
   * Revert a location, translating and flipping as required.
   *
   * <p>
   * If oppositeStrand is false, this is equivalent to un-translating the
   * location. If it is true, this is equivalent to (un)flipping it.
   * </p>
   *
   * @param oldLoc          the Location to revert
   * @param translation     the translation to unapply
   * @param oppositeStrand  wether or not this is a flip
   * @return  the reverted location
   */
  public static Location revertLocation(Location oldLoc,
                                        int translation,
                                        boolean oppositeStrand)
  {
    if (oppositeStrand) {
      return flipLocation(oldLoc, translation);
    } else {
      return oldLoc.translate(-translation);
    }
  }

  /**
   * Flip a location.
   *
   * <p>
   * All points <code>p</code> map to <code>translation - p</code>. Clearly,
   * this mapping is its own inverse. If you wish to flip all locations between
   * 1 and length, you should use a translation of length + 1. In general, if
   * you wish to flip all features between x and y, you should use a translation
   * of x + y.
   * </p>
   *
   * @param oldLoc      the Location to flip
   * @param translation the translation to use
   * @return  the flipped Location
   */
  public static Location flipLocation(Location oldLoc, int translation) {
    if (oldLoc.isContiguous()) {
      if (oldLoc instanceof PointLocation) {
        return new PointLocation(translation - oldLoc.getMin());
      } else {
        return new RangeLocation(translation - oldLoc.getMax(),
                                 translation - oldLoc.getMin());
      }
    } else {
      Location compound;
      List locList = new ArrayList();
      for (Iterator i = oldLoc.blockIterator(); i.hasNext();) {
        Location oldBlock = (Location) i.next();
        locList.add(new RangeLocation(translation - oldBlock.getMax(),
                                      translation - oldBlock.getMin()));
      }
      compound = LocationTools.union(locList);
      return compound;
    }
  }

  public static StrandedFeature.Strand flipStrand(StrandedFeature.Strand s) {
    if (s == StrandedFeature.POSITIVE) {
      return StrandedFeature.NEGATIVE;
    } else if (s == StrandedFeature.NEGATIVE) {
      return StrandedFeature.POSITIVE;
    } else {
      return StrandedFeature.UNKNOWN;
    }
  }
}
