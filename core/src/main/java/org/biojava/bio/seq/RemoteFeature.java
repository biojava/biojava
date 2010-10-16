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

package org.biojava.bio.seq;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.symbol.Location;

/**
 * A feature that indicates that there is some remote feature that can't be
 * represented entirely on a single Sequence.
 * <p>
 * These are the sort of features that represent things like the horible Embl
 * remote feature spans (such as AL24199:100..200). The method getRemoteFeature
 * should return a Feature on another Sequence that properly represents the
 * Location. This Seqeunce will often be an Assembly of the parent Sequence to
 * this Feature, and the Sequence associated with the remote Location sequence
 * ID. Thse Features are also applicable to the case when a portion of a
 * Sequence is projected but some Feature overlaps the boundary of the projected
 * portion.
 *
 * @see org.biojavax.bio.seq.RichFeature
 * @author Matthew Pocock
 * @author Greg Cox
 * @since 1.2
 */
public interface RemoteFeature extends StrandedFeature {
  /**
   * Retrieve the list of Regions that locate this feature both localy and
   * remotely. Local Regions have a null sequence ID.
   *
   * @return an immutable List of Regions
   */
  List getRegions();

  /**
   * Retrieve the Feature on some assembly Sequence that can represent this
   * RemoteFeature properly.
   * <p>
   * This method should be equivalent to calling
   * <code>getResolver().resolve(this)</code>.
   *
   * @return the Feature on another Sequence that this is projected onto
   * @throws BioException if for any reason the remote Feature could not be
   *            constructed
   */
  Feature getRemoteFeature() throws BioException;

  Resolver getResolver();

  public class Template extends StrandedFeature.Template {
    public List regions;
    public Resolver resolver;

    public Template() {
      super();
      location = Location.empty;
      regions = new ArrayList();
      resolver = null;
    }

    /**
     * Creates a RemoteFeature.Template that has the same values as the
     * template passed in.  Fields that are in the template passed in but
     * not in RemoteFeature Templates will be silently discarded.  Regions is
     * set to an empty list and the resolver is set to null.
     *
     * @param theTemplate the template for this template.
     */
    public Template(Feature.Template theTemplate) {
      location = theTemplate.location;
      if (location == null) {
        location = Location.empty;
      }
      type = theTemplate.type;
      source = theTemplate.source;
      annotation = theTemplate.annotation;
      if (theTemplate instanceof StrandedFeature.Template) {
        strand = ((StrandedFeature.Template) theTemplate).strand;
      } else {
        strand = StrandedFeature.UNKNOWN;
      }
      resolver = null;
      regions = null;
    }

  }

  /**
   * The interface for objects that actually can take a RemoteFeature and
   * return a Sequence object with the feature resolved into a real feature.
   * <p>
   * An implementation may choose to create a new Assembly from Sequences
   * in a SequenceDB instance, or to return some existing larger Sequence that
   * the sequence parent of the Feature is part of. This interface should ensure
   * canonicalization of the returned Feature and the Sequence it resides on
   * (using, for instance, SoftReferenceCacheMap keyed by a set of sequence
   * IDs).
   *
   * @author Matthew Pocock
   * @since 1.1
   */
  public static interface Resolver {
    /**
     * Resolve rFeat.
     * <p>
     * This method returns a Feature that represents the RemoteFeature rFeat on
     * some Sequence where all of the Regions can be represented localy. This
     * may be an assembly of the parent sequence of rFeat and each of the
     * Sequences that have IDs listed in the Region list of rFeat.
     *
     * @param rFeat  the RemoteFeature to resolve
     * @return a Feature on some other Seqence where each Region of rFeat is
     *          resolved into a local Location
     */
    Feature resolve(RemoteFeature rFeat) throws IllegalIDException, BioException;
  }

  /**
   * A tuple of Location and sequence ID.
   * <p>
   * For local locations, the Region is just a wrapper for a Location. For
   * remote Regions, it also contains a String representing the Sequence ID of
   * the remote Location.
   *
   * @author Matthew Pocock
   * @since 1.1
   */
  public final static class Region {
    private final Location location;
    private final String seqID;
    private final boolean isRemote;

    /**
     * Create a new Region.
     *
     * @param location  the Location of the Region
     * @param seqID the ID of the Sequence containing the Location, or null if
     *          it is a local Region
     */
    public Region(Location location, String seqID, boolean isRemote) {
      this.location = location;
      this.seqID = seqID;
      this.isRemote = isRemote;
    }

    /**
     * Retrieve the Location of the Region.
     *
     * @return the Location of this Region
     */
    public Location getLocation() {
      return location;
    }

    /**
     * Return the remote Sequence ID if this Region is on another Sequence
     * (isRemote will return true), null otherwise.
     *
     * @return the ID of the remote Sequence containing this Region
     */
    public String getSeqID() {
      return seqID;
    }

    /**
     * Return whether this Region is remote or local.
     * <p>
     * If this returns true, getSeqID() will return the ID of the remote
     * sequence. Otherwise, getSeqID() will return null.
     *
     * @return true if this is a remote Region, false otherwise
     */
    public boolean isRemote() {
      return isRemote;
    }
  }
}
