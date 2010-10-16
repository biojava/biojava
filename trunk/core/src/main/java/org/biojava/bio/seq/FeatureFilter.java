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

import java.io.Serializable;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.AnnotationTools;
import org.biojava.bio.AnnotationType;
import org.biojava.bio.CardinalityConstraint;
import org.biojava.bio.CollectionConstraint;
import org.biojava.bio.PropertyConstraint;
import org.biojava.bio.seq.homol.SimilarityPairFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.walker.WalkerFactory;

/**
 * A filter for accepting or rejecting a feature.
 *
 * <p>
 * It is possible to write custom <code>FeatureFilter</code>s by implementing this
 * interface.  There are also a wide range of built-in features, and it is possible
 * to build complex queries using <code>FeatureFilter.And</code>, <code>FeatureFilter.Or</code>,
 * and <code>FeatureFilter.Not</code>.  Where possible, use of the built-in filters
 * is preferable to writing new filters, since the methods in the <code>FilterUtils</code>
 * class have access to special knowledge about the built-in filter types and how they
 * relate to one another.
 * </p>
 *
 * <p>
 * If the filter is to be used in a remote process, it is recognized that it may
 * be serialized and sent over to run remotely, rather than each feature being
 * retrieved locally.
 * </p>
 *
 * @since 1.0
 * @author Matthew Pocock
 * @author Thomas Down
 */

public interface FeatureFilter extends Serializable {
  /**
   * This method determines whether a feature is to be accepted.
   *
   * @param f the Feature to evaluate
   * @return  true if this feature is to be selected in, or false if it is to be ignored
   */
  boolean accept(Feature f);

  /**
   * All features are selected by this filter.
   */
  static final public FeatureFilter all = new AcceptAllFilter();

  /**
   * No features are selected by this filter.
   */
  static final public FeatureFilter none = new AcceptNoneFilter();


  /**
   *  A filter that returns all features not accepted by a child filter.
   *
   * @author Thomas Down
   * @author Matthew Pocock
   * @since 1.0
   */
  public final static class Not implements FeatureFilter {
    static { WalkerFactory.getInstance().addTypeWithParent(Not.class); }

    FeatureFilter child;

    public FeatureFilter getChild() {
      return child;
    }

    public Not(FeatureFilter child) {
        this.child = child;
    }

    public boolean accept(Feature f) {
        return !(child.accept(f));
    }

    public boolean equals(Object o) {
      return
        (o instanceof Not) &&
        (((Not) o).getChild().equals(this.getChild()));
    }

    public int hashCode() {
      return getChild().hashCode();
    }

    public String toString() {
      return "Not(" + child + ")";
    }
  }

  /**
   *  A filter that returns all features accepted by both child filter.
   *
   * @author Thomas Down
   * @author Matthew Pocock
   * @since 1.0
   */
  public final static class And implements FeatureFilter {
    static { WalkerFactory.getInstance().addTypeWithParent(And.class); }

    FeatureFilter c1, c2;

    public FeatureFilter getChild1() {
      return c1;
    }

    public FeatureFilter getChild2() {
      return c2;
    }

    public And(FeatureFilter c1, FeatureFilter c2) {
        this.c1 = c1;
        this.c2 = c2;
    }

    public boolean accept(Feature f) {
        return (c1.accept(f) && c2.accept(f));
    }

    public boolean equals(Object o) {
      if(o instanceof FeatureFilter) {
        return FilterUtils.areEqual(this, (FeatureFilter) o);
      } else {
        return false;
      }
    }

    public int hashCode() {
      return getChild1().hashCode() ^ getChild2().hashCode();
    }

      public String toString() {
	  return "And(" + c1 + " , " + c2 + ")";
       }
  }

  /**
   *  A filter that returns all features accepted by at least one child filter.
   *
   * @author Thomas Down
   * @author Matthew Pocock
   * @since 1.0
   */
  public final static class Or implements FeatureFilter {
    static { WalkerFactory.getInstance().addTypeWithParent(Or.class); }

    FeatureFilter c1, c2;

    public FeatureFilter getChild1() {
      return c1;
    }

    public FeatureFilter getChild2() {
      return c2;
    }

    public Or(FeatureFilter c1, FeatureFilter c2) {
        this.c1 = c1;
        this.c2 = c2;
    }

    public boolean accept(Feature f) {
        return (c1.accept(f) || c2.accept(f));
    }

    public boolean equals(Object o) {
      if(o instanceof FeatureFilter) {
        return FilterUtils.areEqual(this, (FeatureFilter) o);
      } else {
        return false;
      }
    }

    public int hashCode() {
      return getChild1().hashCode() ^ getChild2().hashCode();
    }

    public String toString() {
      return "Or(" + c1 + " , " + c2 + ")";
    }
  }

  /**
   * Construct one of these to filter features by type.
   *
   * @author Matthew Pocock
   * @since 1.0
   */
  final public static class ByType implements OptimizableFilter {
    private String type;

    public String getType() {
      return type;
    }

    /**
     * Create a ByType filter that filters in all features with type fields
     * equal to type.
     *
     * @param type  the String to match type fields against
     */
    public ByType(String type) {
        if (type == null)
            throw new NullPointerException("Type may not be null");
        this.type = type;
    }

    /**
     * Returns true if the feature has a matching type property.
     */
    public boolean accept(Feature f) {
      return type.equals(f.getType());
    }

    public boolean equals(Object o) {
      return
        (o instanceof ByType) &&
        (((ByType) o).getType().equals(this.getType()));
    }

    public int hashCode() {
      return getType().hashCode();
    }

    public boolean isProperSubset(FeatureFilter sup) {
      return this.equals(sup) || (sup instanceof AcceptAllFilter);
    }

    public boolean isDisjoint(FeatureFilter filt) {
      return (filt instanceof AcceptNoneFilter) || (
        (filt instanceof ByType) &&
        !getType().equals(((ByType) filt).getType())
      );
    }

    public String toString() {
      return "ByType(" + type + ")";
    }
  }

  /**
   * Construct one of these to filter features by source.
   *
   * @author Matthew Pocock
   * @since 1.0
   */
  public final static class BySource implements OptimizableFilter {
    private String source;

    public String getSource() {
      return source;
    }

    /**
     * Create a BySource filter that filters in all features which have sources
     * equal to source.
     *
     * @param source  the String to match source fields against
     */
    public BySource(String source) {
        if (source == null)
            throw new NullPointerException("Source may not be null");
        this.source = source;
    }

    public boolean accept(Feature f) { return source.equals(f.getSource()); }

    public boolean equals(Object o) {
      return
        (o instanceof BySource) &&
        (((BySource) o).getSource().equals(this.getSource()));
    }

    public boolean isProperSubset(FeatureFilter sup) {
      return this.equals(sup) || (sup instanceof AcceptAllFilter);
    }

    public int hashCode() {
      return getSource().hashCode();
    }

    public boolean isDisjoint(FeatureFilter filt) {
      return (filt instanceof AcceptNoneFilter) || (
        (filt instanceof BySource) &&
        !getSource().equals(((BySource) filt).getSource())
      );
    }

      public String toString() {
	  return "BySource(" + source + ")";
       }
  }

  /**
   * Filter which accepts only those filters which are an instance
   * of a specific Java class
   *
   * @author Thomas Down
   * @author Matthew Pocock
   * @since 1.1
   */

  public final static class ByClass implements OptimizableFilter {
    private Class clazz;

    public ByClass(Class clazz) {
        if (clazz == null) {
            throw new NullPointerException("Clazz may not be null");
        }
        if(!Feature.class.isAssignableFrom(clazz)) {
            throw new ClassCastException(
              "Filters by class must be over Feature classes: " +
              clazz
            );
        }
        this.clazz = clazz;
    }

    public boolean accept(Feature f) {
      return clazz.isInstance(f);
    }

    public Class getTestClass() {
      return clazz;
    }

    public boolean equals(Object o) {
      return
        (o instanceof ByClass) &&
        (((ByClass) o).getTestClass() == this.getTestClass());
    }

    public int hashCode() {
      return getTestClass().hashCode();
    }

    public boolean isProperSubset(FeatureFilter sup) {
      if(sup instanceof ByClass) {
        Class supC = ((ByClass) sup).getTestClass();
        return supC.isAssignableFrom(this.getTestClass());
      }
      return (sup instanceof AcceptAllFilter);
    }

    public boolean isDisjoint(FeatureFilter feat) {
      if(feat instanceof ByClass) {
        Class featC = ((ByClass) feat).getTestClass();
        return
          ! (featC.isAssignableFrom(getTestClass())) &&
          ! (getTestClass().isAssignableFrom(featC));
      } else if (feat instanceof ByComponentName) {
	  return !getTestClass().isAssignableFrom(ComponentFeature.class);
      }

      return (feat instanceof AcceptNoneFilter);
    }

    public String toString() {
      return "ByClass(" + clazz.getName() + ")";
    }
  }


  /**
   * Accept features with a given strandedness.
   *
   * @author Matthew Pocock
   * @since 1.1
   */
  public final static class StrandFilter implements OptimizableFilter {
    private StrandedFeature.Strand strand;

    /**
     * Build a new filter that matches all features of a given strand.
     *
     * @param strand the Strand to match
     */
    public StrandFilter(StrandedFeature.Strand strand) {
      this.strand = strand;
    }

    /**
     * Retrieve the strand this matches.
     *
     * @return the Strand matched
     */
    public StrandedFeature.Strand getStrand() {
      return strand;
    }

    /**
     * Accept the Feature if it is an instance of StrandedFeature and matches
     * the value of getStrand().
     *
     * @param f the Feature to check
     * @return true if the strand matches, or false otherwise
     */
    public boolean accept(Feature f) {
      if(f instanceof StrandedFeature) {
        StrandedFeature sf = (StrandedFeature) f;
        return sf.getStrand() == strand;
      } else {
        return strand == StrandedFeature.UNKNOWN;
      }
    }

    public boolean equals(Object o) {
      return
        (o instanceof StrandFilter) &&
        (((StrandFilter) o).getStrand() == this.getStrand());
    }

    public int hashCode() {
      return getStrand().hashCode();
    }

    public String toString() {
      return "StrandedFilter(" + strand + ")";
    }

    public boolean isProperSubset(FeatureFilter sup) {
      return this.equals(sup);
    }

    public boolean isDisjoint(FeatureFilter filt) {
      return (filt instanceof AcceptNoneFilter) || (
        (filt instanceof StrandFilter) &&
        ((StrandFilter) filt).getStrand() == getStrand()
      );
    }
  }

  /**
   * Accept features that reside on a sequence with a particular name.
   *
   * @author Matthew Pocock
   * @since 1.3
   */
  public final static class BySequenceName
  implements OptimizableFilter {
    private String seqName;

    public BySequenceName(String seqName) {
      this.seqName = seqName;
    }

    public String getSequenceName() {
      return seqName;
    }

    public boolean accept(Feature f) {
      return f.getSequence().getName().equals(seqName);
    }

    public boolean isProperSubset(FeatureFilter sup) {
      return equals(sup);
    }

    public boolean isDisjoint(FeatureFilter filt) {
        if (filt instanceof BySequenceName) {
            return !equals(this);
        } else {
            return false;
        }
    }

    public boolean equals(Object o) {
      return
       (o instanceof BySequenceName) &&
       ((BySequenceName) o).getSequenceName().equals(seqName);
    }

    public int hashCode() {
      return seqName.hashCode();
    }
  }

  /**
   *  A filter that returns all features contained within a location.
   *
   * @author Matthew Pocock
   * @since 1.0
   */
  public final static class ContainedByLocation implements OptimizableFilter {
    private Location loc;

    public Location getLocation() {
      return loc;
    }

    /**
     * Creates a filter that returns everything contained within loc.
     *
     * @param loc  the location that will contain the accepted features
     */
    public ContainedByLocation(Location loc) {
        if (loc == null) {
            throw new NullPointerException("Loc may not be null");
        }
        this.loc = loc;
    }

    /**
     * Returns true if the feature is within this filter's location.
     */
    public boolean accept(Feature f) {
      return loc.contains(f.getLocation());
    }

    public boolean equals(Object o) {
      return
        (o instanceof ContainedByLocation) &&
        (((ContainedByLocation) o).getLocation().equals(this.getLocation()));
    }

    public int hashCode() {
      return getLocation().hashCode();
    }

    public boolean isProperSubset(FeatureFilter sup) {
      if(sup instanceof ContainedByLocation) {
        Location supL = ((ContainedByLocation) sup).getLocation();
        return supL.contains(this.getLocation());
      } else if(sup instanceof OverlapsLocation) {
        Location supL = ((OverlapsLocation) sup).getLocation();
        return supL.contains(this.getLocation());
      } else if (sup instanceof ShadowOverlapsLocation) {
        Location supL = ((ShadowOverlapsLocation) sup).getLocation();
        return supL.contains(this.getLocation());
      } else if (sup instanceof ShadowContainedByLocation) {
        Location supL = ((ShadowContainedByLocation) sup).getLocation();
        return supL.contains(this.getLocation());
      }
      return (sup instanceof AcceptAllFilter);
    }

    public boolean isDisjoint(FeatureFilter filt) {
      if(filt instanceof ContainedByLocation) {
        Location loc = ((ContainedByLocation) filt).getLocation();
        return !getLocation().overlaps(loc);
      } else if (filt instanceof OverlapsLocation) {
	    Location filtL = ((OverlapsLocation) filt).getLocation();
	    return !filtL.overlaps(this.getLocation());
      } else if (filt instanceof ShadowOverlapsLocation) {
        Location filtL = ((ShadowOverlapsLocation) filt).getLocation();
        return filtL.getMax() < loc.getMin() || filtL.getMin() > loc.getMax();
      } else if (filt instanceof ShadowContainedByLocation) {
        Location filtL = ((ShadowContainedByLocation) filt).getLocation();
        return filtL.getMax() < loc.getMin() || filtL.getMin() > loc.getMax();
      }

      return (filt instanceof AcceptNoneFilter);
    }

    public String toString() {
      return "ContainedBy(" + loc + ")";
    }
  }

  /**
   *  A filter that returns all features overlapping a location.
   *
   * @author Matthew Pocock
   * @since 1.0
   */
  public final static class OverlapsLocation implements OptimizableFilter {
    private Location loc;

    public Location getLocation() {
      return loc;
    }

    /**
     * Creates a filter that returns everything overlapping loc.
     *
     * @param loc  the location that will overlap the accepted features
     */
    public OverlapsLocation(Location loc) {
        if (loc == null) {
            throw new NullPointerException("Loc may not be null");
        }
        this.loc = loc;
    }

    /**
     * Returns true if the feature overlaps this filter's location.
     */
    public boolean accept(Feature f) {
      return loc.overlaps(f.getLocation());
    }

    public boolean equals(Object o) {
      return
        (o instanceof OverlapsLocation) &&
        (((OverlapsLocation) o).getLocation().equals(this.getLocation()));
    }

    public int hashCode() {
      return getLocation().hashCode();
    }

    public boolean isProperSubset(FeatureFilter sup) {
      if(sup instanceof OverlapsLocation) {
          Location supL = ((OverlapsLocation) sup).getLocation();
          return supL.contains(this.getLocation());
      } else if (sup instanceof ShadowOverlapsLocation) {
          Location supL = ((ShadowOverlapsLocation) sup).getLocation();
          return supL.contains(this.getLocation());
      }
      return (sup instanceof AcceptAllFilter);
    }

    public boolean isDisjoint(FeatureFilter filt) {
      if (filt instanceof ContainedByLocation)  {
        Location loc = ((ContainedByLocation) filt).getLocation();
        return !getLocation().overlaps(loc);
      } else if (filt instanceof ShadowContainedByLocation) {
        Location loc = ((ShadowContainedByLocation) filt).getLocation();
        return getLocation().getMax() < loc.getMin() || getLocation().getMin() > loc.getMax();
      }
      return (filt instanceof AcceptNoneFilter);
    }

    public String toString() {
      return "Overlaps(" + loc + ")";
    }
  }

  /**
   *  A filter that accepts all features whose shadow overlaps a specified
   * <code>Location</code>.  The shadow is defined as the interval between the
   * minimum and maximum positions of the feature's location.  For features
   * with contiguous locations, this filter is equivalent to
   * <code>FeatureFilter.OverlapsLocation</code>..
   *
   * <p>
   * A typical use of this filter is in graphics code where you are rendering
   * features with non-contiguous locations in a `blocks and connectors' style,
   * and wish to draw the connector even when no blocks fall within the
   * selected field of view
   * </p>
   *
   * @author Thomas Down
   * @since 1.3
   */

  public final static class ShadowOverlapsLocation implements OptimizableFilter {
    private Location loc;

    public Location getLocation() {
      return loc;
    }

    /**
     * Creates a filter that returns everything overlapping loc.
     *
     * @param loc  the location that will overlap the accepted features
     */
    public ShadowOverlapsLocation(Location loc) {
        if (loc == null) {
            throw new NullPointerException("Loc may not be null");
        }
        this.loc = loc;
    }

    /**
     * Returns true if the feature overlaps this filter's location.
     */
    public boolean accept(Feature f) {
        Location floc = f.getLocation();
        if (!floc.isContiguous()) {
            floc = new RangeLocation(floc.getMin(), floc.getMax());
        }
        return loc.overlaps(floc);
    }

    public boolean equals(Object o) {
      return
        (o instanceof ShadowOverlapsLocation) &&
        (((ShadowOverlapsLocation) o).getLocation().equals(this.getLocation()));
    }

    public int hashCode() {
      return getLocation().hashCode() +77;
    }

    public boolean isProperSubset(FeatureFilter sup) {
      if(sup instanceof ShadowOverlapsLocation) {
        Location supL = ((ShadowOverlapsLocation) sup).getLocation();
        return supL.contains(this.getLocation());
      }
      return (sup instanceof AcceptAllFilter);
    }

    public boolean isDisjoint(FeatureFilter filt) {
      if (filt instanceof ShadowContainedByLocation)  {
          Location loc = ((ShadowContainedByLocation) filt).getLocation();
          return !getLocation().overlaps(loc);
      }  else if (filt instanceof ContainedByLocation) {
          Location loc = ((ContainedByLocation) filt).getLocation();
          return (loc.getMax() < getLocation().getMin() || loc.getMin() > getLocation().getMax());
      }
      return (filt instanceof AcceptNoneFilter);
    }

    public String toString() {
      return "ShadowOverlaps(" + loc + ")";
    }
  }

  /**
   *  A filter that accepts all features whose shadow is contained by a specified
   * <code>Location</code>.  The shadow is defined as the interval between the
   * minimum and maximum positions of the feature's location.  For features
   * with contiguous locations, this filter is equivalent to
   * <code>FeatureFilter.ContainedByLocation</code>.
   *
   * @author Thomas Down
   * @since 1.3
   */

  public final static class ShadowContainedByLocation implements OptimizableFilter {
    private Location loc;

    public Location getLocation() {
      return loc;
    }

    /**
     * Creates a filter that returns everything contained within loc.
     *
     * @param loc  the location that will contain the accepted features
     */
    public ShadowContainedByLocation(Location loc) {
        if (loc == null) {
            throw new NullPointerException("Loc may not be null");
        }
        this.loc = loc;
    }

    /**
     * Returns true if the feature is within this filter's location.
     */
    public boolean accept(Feature f) {
        Location floc = f.getLocation();
        if (!floc.isContiguous()) {
            floc = new RangeLocation(floc.getMin(), floc.getMax());
        }
        return loc.contains(floc);
    }

    public boolean equals(Object o) {
      return
        (o instanceof ShadowContainedByLocation) &&
        (((ShadowContainedByLocation) o).getLocation().equals(this.getLocation()));
    }

    public int hashCode() {
      return getLocation().hashCode() + 88;
    }

    public boolean isProperSubset(FeatureFilter sup) {
      if(sup instanceof ShadowContainedByLocation) {
        Location supL = ((ShadowContainedByLocation) sup).getLocation();
        return supL.contains(this.getLocation());
      } else if (sup instanceof ShadowOverlapsLocation) {
        Location supL = ((ShadowOverlapsLocation) sup).getLocation();
        return supL.contains(this.getLocation());
      }
      return (sup instanceof AcceptAllFilter);
    }

    public boolean isDisjoint(FeatureFilter filt) {
      if(filt instanceof ContainedByLocation) {
        Location filtL = ((ShadowContainedByLocation) filt).getLocation();
        return filtL.getMax() < loc.getMin() || filtL.getMin() > loc.getMax();
      } if(filt instanceof ShadowContainedByLocation) {
        Location loc = ((ShadowContainedByLocation) filt).getLocation();
        return !getLocation().overlaps(loc);
      } else if (filt instanceof OverlapsLocation) {
	    Location filtL = ((OverlapsLocation) filt).getLocation();
	    return filtL.getMax() < loc.getMin() || filtL.getMin() > loc.getMax();
      } else if (filt instanceof ShadowOverlapsLocation) {
        Location filtL = ((ShadowOverlapsLocation) filt).getLocation();
        return !filtL.overlaps(this.getLocation());
      }

      return (filt instanceof AcceptNoneFilter);
    }

    public String toString() {
      return "ShadowContainedBy(" + loc + ")";
    }
  }

  /**
   * A filter that returns all features that have an annotation bundle that is of a given
   * annotation type.
   *
   * @author Matthew Pocock
   * @author Thomas Down
   * @since 1.3
   */
  public static class ByAnnotationType
  implements OptimizableFilter {
    private AnnotationType type;

    protected ByAnnotationType() {
      this(AnnotationType.ANY);
    }

    public ByAnnotationType(AnnotationType type) {
      this.type = type;
    }

    public AnnotationType getType() {
      return type;
    }

    protected void setType(AnnotationType type) {
      this.type = type;
    }

    public boolean accept(Feature f) {
      return type.instanceOf(f.getAnnotation());
    }

    public boolean equals(Object o) {
      if(o instanceof ByAnnotationType) {
        ByAnnotationType that = (ByAnnotationType) o;
        return this.getType() == that.getType();
      }

      return false;
    }

    public int hashCode() {
      return getType().hashCode();
    }

    public boolean isDisjoint(FeatureFilter filter) {
      if(filter instanceof AcceptNoneFilter) {
        return true;
      } else if(filter instanceof ByAnnotationType) {
        // check for common property names
        ByAnnotationType that = (ByAnnotationType) filter;
        Set props = that.getType().getProperties();
        Set ourProps = this.getType().getProperties();
        Set allProps = new HashSet(props);
        allProps.addAll(ourProps);
        for(Iterator i = allProps.iterator(); i.hasNext(); ) {
          Object prop = i.next();

          CollectionConstraint thisC = this.getType().getConstraint(prop);
          CollectionConstraint thatC = that.getType().getConstraint(prop);
          if (AnnotationTools.intersection(thisC, thatC) == CollectionConstraint.NONE) {
              return true;
          }
        }
      }

      return false;
    }

    public boolean isProperSubset(FeatureFilter filter) {
      if(filter instanceof ByAnnotationType) {
        ByAnnotationType that = (ByAnnotationType) filter;

        Set thisProps = this.getType().getProperties();
        Set thatProps = that.getType().getProperties();
        for(Iterator i = thatProps.iterator(); i.hasNext(); ) {
          Object prop = i.next();

          if(!thisProps.contains(prop)) {
            return false;
          }

          CollectionConstraint thisP = this.getType().getConstraint(prop);
          CollectionConstraint thatP = that.getType().getConstraint(prop);

          if(
            !thatP.subConstraintOf(thisP)
          ) {
            return false;
          }
        }

        return true;
      }

      return false;
    }

    public String toString() {
      return "ByAnnotationType {" + type + "}";
    }
  }

  /**
   * Retrieve features that contain a given annotation with a given value.
   *
   * @author Matthew Pocock
   * @author Keith James
   * @since 1.1
   */
  public final static class ByAnnotation
  extends ByAnnotationType {
    private Object key;
    private Object value;

    /**
     * Make a new ByAnnotation that will accept features with an annotation
     * bundle containing 'value' associated with 'key'.
     *
     * @param key  the Object used as a key in the annotation
     * @param value the Object associated with key in the annotation
     */
    public ByAnnotation(Object key, Object value) {
      this.key = key;
      this.value = value;

      AnnotationType.Impl type = new AnnotationType.Impl();
      type.setConstraints(
        key,
        new PropertyConstraint.ExactValue(value),
        CardinalityConstraint.ONE
      );
      setType(type);
    }

    public Object getKey() {
      return key;
    }

    public Object getValue() {
      return value;
    }
  }

  /**
   * Retrieve features that contain a given annotation, and that the set of values
   * contains the value given.
   *
   * @author Thomas Down
   * @since 1.3
   */
  public final static class AnnotationContains
  extends ByAnnotationType {
    private Object key;
    private Object value;

    /**
     * Make a new AnnotationContains that will accept features with an annotation
     * bundle where the value-set assosiated with the property <code>key</code>
     * contains a member equal to <code>value</code>.
     *
     * @param key  the Object used as a key in the annotation
     * @param value the Object associated with key in the annotation
     */
    public AnnotationContains(Object key, Object value) {
      this.key = key;
      this.value = value;

      AnnotationType.Impl type = new AnnotationType.Impl();
      type.setConstraint(
        key,
        new CollectionConstraint.Contains(
            new PropertyConstraint.ExactValue(value),
            CardinalityConstraint.ONE
        )
      );
      setType(type);
    }

    public Object getKey() {
      return key;
    }

    public Object getValue() {
      return value;
    }
  }

  /**
   * Retrieve features that contain a given annotation with any value.
   *
   * @author Matthew Pocock
   * @author Keith James
   * @since 1.1
   */
  public final static class HasAnnotation
  extends ByAnnotationType {
    private Object key;

    /**
     * Make a new ByAnnotation that will accept features with an annotation
     * bundle containing any value associated with 'key'.
     *
     * @param key  the Object used as a key in the annotation
     */
    public HasAnnotation(Object key) {
      this.key = key;

      AnnotationType.Impl type = new AnnotationType.Impl();
      type.setConstraints(
        key,
        PropertyConstraint.ANY,
        CardinalityConstraint.ONE_OR_MORE
      );
      setType(type);
    }

    public Object getKey() {
      return key;
    }
  }

    /**
     * Filter by applying a nested <code>FeatureFilter</code> to the
     * parent feature.  Always <code>false</code> if the parent
     * is not a feature (e.g. top-level features, where the
     * parent is a sequence).
     *
     * @author Thomas Down
     * @since 1.2
     */

    public static class ByParent implements OptimizableFilter, Up {
      static { WalkerFactory.getInstance().addTypeWithParent(ByParent.class); }

        private FeatureFilter filter;

        public ByParent(FeatureFilter ff) {
            filter = ff;
        }

        public FeatureFilter getFilter() {
            return filter;
        }

        public boolean accept(Feature f) {
            FeatureHolder fh = f.getParent();
            if (fh instanceof Feature) {
                return filter.accept((Feature) fh);
            }

            return false;
        }

        public int hashCode() {
            return filter.hashCode() + 173;
        }

        public boolean equals(Object o) {
            if (! (o instanceof FeatureFilter.ByParent)) {
                return false;
            }

            FeatureFilter.ByParent ffbp = (FeatureFilter.ByParent) o;
            return ffbp.getFilter().equals(filter);
        }

        public boolean isProperSubset(FeatureFilter ff) {
            FeatureFilter ancFilter = null;
            if (ff instanceof FeatureFilter.ByParent) {
                ancFilter = ((FeatureFilter.ByParent) ff).getFilter();
            } else if (ff instanceof FeatureFilter.ByAncestor) {
                ancFilter = ((FeatureFilter.ByAncestor) ff).getFilter();
            }

            if (ancFilter != null) {
                return FilterUtils.areProperSubset(ancFilter, filter);
            } else {
                return false;
            }
        }

        public boolean isDisjoint(FeatureFilter ff) {
            // System.err.println("Disjunction test for " + toString());
            // System.err.println("Against " + ff.toString());
            if (ff instanceof IsTopLevel) {
                return true;
            } else if (ff instanceof FeatureFilter.ByParent) {
                return FilterUtils.areDisjoint(
                        ((FeatureFilter.ByParent) ff).getFilter(),
                        getFilter()
                );
            }  else if (ff instanceof FeatureFilter.ByAncestor) {
                return FilterUtils.areDisjoint(
                        ((FeatureFilter.ByAncestor) ff).getFilter(),
                        getFilter()
                );
            } else {
                FeatureFilter childFilter = FilterUtils.getOnlyChildrenFilter(getFilter());
                if (childFilter != null) {
                    return FilterUtils.areDisjoint(
                            childFilter,
                            ff
                    );
                }
            }

            return false;
        }
    }

    /**
     * Filter by applying a nested <code>FeatureFilter</code> to all
     * ancestor features.  Returns <code>true</code> if at least one
     * of them matches the filter.  Always <code>false</code> if the
     * parent is not a feature (e.g. top-level features, where the
     * parent is a sequence).
     *
     * @author Thomas Down
     * @since 1.2
     */

    public static class ByAncestor implements OptimizableFilter, Up {
      static { WalkerFactory.getInstance().addTypeWithParent(ByAncestor.class); }

        private FeatureFilter filter;

        public ByAncestor(FeatureFilter ff) {
            filter = ff;
        }

        public FeatureFilter getFilter() {
            return filter;
        }

        public boolean accept(Feature f) {
            do {
                FeatureHolder fh = f.getParent();
                if (fh instanceof Feature) {
                    f = (Feature) fh;
                    if (filter.accept(f)) {
                        return true;
                    }
                } else {
                    return false;
                }
            } while (true);
        }

        public int hashCode() {
            return filter.hashCode() + 186;
        }

        public boolean equals(Object o) {
            if (! (o instanceof FeatureFilter.ByAncestor)) {
                return false;
            }

            FeatureFilter.ByAncestor ffba = (FeatureFilter.ByAncestor) o;
            return ffba.getFilter().equals(filter);
        }

        public boolean isProperSubset(FeatureFilter ff) {
            FeatureFilter ancFilter = null;
            if (ff instanceof FeatureFilter.ByAncestor) {
                ancFilter = ((FeatureFilter.ByAncestor) ff).getFilter();
            }

            if (ancFilter != null) {
                return FilterUtils.areProperSubset(ancFilter, filter);
            } else {
                return false;
            }
        }

        public boolean isDisjoint(FeatureFilter ff) {
            // System.err.println("Disjunction test for " + toString());
            // System.err.println("Against " + ff.toString());
            if (ff instanceof IsTopLevel) {
                return true;
            }

            if (ff instanceof FeatureFilter.ByParent) {
                return FilterUtils.areDisjoint(
                        ((FeatureFilter.ByParent) ff).getFilter(),
                        getFilter()
                );
            }  else if (ff instanceof FeatureFilter.ByAncestor) {
                return FilterUtils.areDisjoint(
                        ((FeatureFilter.ByAncestor) ff).getFilter(),
                        getFilter()
                );
            } else {
                FeatureFilter descFilter = FilterUtils.getOnlyDescendantsFilter(getFilter());
                if (descFilter != null) {
                    return FilterUtils.areDisjoint(
                            descFilter,
                            ff
                    );
                }

                FeatureFilter childFilter = FilterUtils.getOnlyChildrenFilter(getFilter());
                if (childFilter != null) {
                    // System.err.println("Child filter is: " + childFilter.toString());
                    if (FilterUtils.areProperSubset(childFilter, leaf)) {
                        // System.err.println("Leaf case");
                        return FilterUtils.areDisjoint(
                                childFilter,
                                ff
                        );
                    } else {
                        // System.err.println("Tree case");
                        return FilterUtils.areDisjoint(
                                new FeatureFilter.Or(
                                        childFilter,
                                        new FeatureFilter.ByAncestor(childFilter)
                                ),
                                ff
                        );
                    }
                }
            }

            return false;
        }

        public String toString() {
            return "ByAncestor(" + getFilter().toString() + ")";
        }
    }

    /**
     * Accepts features where all immediate children meet the supplied filter.  This
     * will be <code>true</code> in the case where no child features exist.  Mainly useful
     * for defining schemas of feature-trees.
     *
     * @author Thomas Down
     * @since 1.3
     */

    public static class OnlyChildren implements OptimizableFilter, ByHierarchy {
      static { WalkerFactory.getInstance().addTypeWithParent(OnlyChildren.class); }

        private FeatureFilter filter;

        public OnlyChildren(FeatureFilter ff) {
            this.filter = ff;
        }

        public FeatureFilter getFilter() {
            return filter;
        }

        public boolean accept(Feature f) {
            for (Iterator i = f.features(); i.hasNext(); ) {
                if (!filter.accept((Feature) i.next())) {
                    return false;
                }
            }
            return true;
        }

        public int hashCode() {
            return filter.hashCode() + 762;
        }

        public boolean equals(Object o) {
            if (! (o instanceof FeatureFilter.OnlyChildren)) {
                return false;
            }

            FeatureFilter.OnlyChildren ffoc = (FeatureFilter.OnlyChildren) o;
            return ffoc.getFilter().equals(filter);
        }

        public boolean isProperSubset(FeatureFilter ff) {
            if (ff == FeatureFilter.all) {
                return true;
            } else if (ff instanceof OnlyChildren) {
                return FilterUtils.areProperSubset(
                    getFilter(),
                    ((OnlyChildren) ff).getFilter()
                ) ;
            } else if (ff instanceof OnlyDescendants) {
                return FilterUtils.areProperSubset(
                    getFilter(),
                    ((OnlyDescendants) ff).getFilter()
                ) ;
            } else {
                return false;
            }
        }

        public boolean isDisjoint(FeatureFilter ff) {
            if (ff instanceof ByChild) {
                return FilterUtils.areDisjoint(
                    getFilter(),
                    ((ByChild) ff).getFilter()
                );
            } else {
                return false;
            }
        }

        public String toString() {
            return "OnlyChildren(" + filter.toString() + ")";
        }
    }

    /**
     * Accepts features where all descendants meet the supplied filter.  This
     * will be <code>true</code> in the case where no child features exist.  Mainly useful
     * for defining schemas of feature-trees.
     *
     * @author Thomas Down
     * @since 1.3
     */

    public static class OnlyDescendants implements OptimizableFilter, ByHierarchy {
      static { WalkerFactory.getInstance().addTypeWithParent(OnlyDescendants.class); }

        private FeatureFilter filter;

        public OnlyDescendants(FeatureFilter ff) {
            this.filter = ff;
        }

        public FeatureFilter getFilter() {
            return filter;
        }

        public boolean accept(Feature f) {
            return f.filter(FeatureFilter.all).countFeatures() == f.filter(filter).countFeatures();
        }

        public int hashCode() {
            return filter.hashCode() + 763;
        }

        public boolean equals(Object o) {
            if (! (o instanceof FeatureFilter.OnlyDescendants)) {
                return false;
            }

            FeatureFilter.OnlyDescendants ffoc = (FeatureFilter.OnlyDescendants) o;
            return ffoc.getFilter().equals(filter);
        }


        public boolean isProperSubset(FeatureFilter ff) {
            if (ff == FeatureFilter.all) {
                return true;
            } else if (ff instanceof OnlyDescendants) {
                return FilterUtils.areProperSubset(
                    getFilter(),
                    ((OnlyDescendants) ff).getFilter()
                ) ;
            } else {
                return false;
            }
        }

        public boolean isDisjoint(FeatureFilter ff) {
            if (ff instanceof ByChild) {
                return FilterUtils.areDisjoint(
                    getFilter(),
                    ((ByChild) ff).getFilter()
                );
            } else if (ff instanceof ByDescendant) {
                return FilterUtils.areDisjoint(
                    getFilter(),
                    ((ByDescendant) ff).getFilter()
                );
            } else {
                return false;
            }
        }
    }

    /**
     * Filter by applying a nested <code>FeatureFilter</code> to the
     * child features.  Always <code>false</code> if there are no children.
     *
     * @author Matthew Pocock
     * @author Thomas Down
     * @since 1.3
     */

    public static class ByChild implements OptimizableFilter, Down {
      static { WalkerFactory.getInstance().addTypeWithParent(ByChild.class); }

        private FeatureFilter filter;

        public ByChild(FeatureFilter ff) {
            filter = ff;
        }

        public FeatureFilter getFilter() {
            return filter;
        }

        public boolean accept(Feature f) {
          for(Iterator i = f.features(); i.hasNext(); ) {
            if(filter.accept((Feature) i.next())) {
              return true;
            }
          }

          return false;
        }

        public int hashCode() {
            return filter.hashCode() + 173;
        }

        public boolean equals(Object o) {
            if (! (o instanceof FeatureFilter.ByChild)) {
                return false;
            }

            FeatureFilter.ByChild ffbc = (FeatureFilter.ByChild) o;
            return ffbc.getFilter().equals(filter);
        }

        public boolean isProperSubset(FeatureFilter ff) {
            FeatureFilter descFilter = null;
            if (ff instanceof FeatureFilter.ByChild) {
                descFilter = ((FeatureFilter.ByChild) ff).getFilter();
            } else if (ff instanceof FeatureFilter.ByDescendant) {
                descFilter = ((FeatureFilter.ByDescendant) ff).getFilter();
            }

            if (descFilter != null) {
                return FilterUtils.areProperSubset(descFilter, filter);
            } else {
                return false;
            }
        }

        public boolean isDisjoint(FeatureFilter ff) {
            if (ff instanceof OnlyChildren) {
                return FilterUtils.areDisjoint(
                    getFilter(),
                    ((OnlyChildren) ff).getFilter()
                );
            } else if (ff instanceof OnlyDescendants) {
                return FilterUtils.areDisjoint(
                    getFilter(),
                    ((OnlyDescendants) ff).getFilter()
                );
            } else {
                return false;
            }
        }
    }


    /**
     * Filter by applying a nested <code>FeatureFilter</code> to all
     * descendant features.  Returns <code>true</code> if at least one
     * of them matches the filter.  Always <code>false</code> if the
     * feature has no children.
     *
     * @author Matthew Pocock
     * @author Thomas Down
     * @since 1.2
     */

    public static class ByDescendant implements OptimizableFilter, Down {
      static { WalkerFactory.getInstance().addTypeWithParent(ByDescendant.class); }

        private FeatureFilter filter;

        public ByDescendant(FeatureFilter ff) {
            filter = ff;
        }

        public FeatureFilter getFilter() {
            return filter;
        }

        public boolean accept(Feature f) {
            do {
                FeatureHolder fh = f.getParent();
                if (fh instanceof Feature) {
                    f = (Feature) fh;
                    if (filter.accept(f)) {
                        return true;
                    }
                } else {
                    return false;
                }
            } while (true);
        }

        public int hashCode() {
            return filter.hashCode() + 186;
        }

        public boolean equals(Object o) {
            if (! (o instanceof FeatureFilter.ByDescendant)) {
                return false;
            }

            FeatureFilter.ByDescendant ffba = (FeatureFilter.ByDescendant) o;
            return ffba.getFilter().equals(filter);
        }

        public boolean isProperSubset(FeatureFilter ff) {
            FeatureFilter ancFilter = null;
            if (ff instanceof FeatureFilter.ByDescendant) {
                ancFilter = ((FeatureFilter.ByDescendant) ff).getFilter();
            }

            if (ancFilter != null) {
                return FilterUtils.areProperSubset(ancFilter, filter);
            } else {
                return false;
            }
        }

        public boolean isDisjoint(FeatureFilter ff) {
           if (ff instanceof OnlyDescendants) {
                return FilterUtils.areDisjoint(
                    getFilter(),
                    ((OnlyDescendants) ff).getFilter()
                );
            } else {
                return false;
            }
        }
    }

  /**
   * Accept features with a given reading frame.
   *
   * @author Mark Schreiber
   * @since 1.2
   */
  public final static class FrameFilter implements OptimizableFilter {
    private FramedFeature.ReadingFrame frame;

    /**
     * Build a new filter that matches all features of a reading frame.
     *
     * @param frame the ReadingFrame to match
     */
    public FrameFilter(FramedFeature.ReadingFrame frame) {
      this.frame = frame;
    }

    /**
     * Retrieve the reading frame this filter matches.
     */
     public FramedFeature.ReadingFrame getFrame(){
       return frame;
     }

    /**
     * Accept the Feature if it is an instance of FramedFeature and matches
     * the value of getFrame().
     *
     * @param f the Feature to check
     * @return true if the frame matches, or false otherwise
     */
    public boolean accept(Feature f) {
      if(f instanceof FramedFeature) {
        FramedFeature ff = (FramedFeature) f;
        return ff.getReadingFrame() == frame;
      } else {
        return false;
      }
    }

    public int hashCode() {
        return frame.getFrame() + 99;
    }

    public boolean equals(Object o) {
      return (o instanceof FrameFilter && ((FrameFilter) o).getFrame() == getFrame());
    }

    public boolean isProperSubset(FeatureFilter sup) {
      return this.equals(sup);
    }

    public boolean isDisjoint(FeatureFilter filt) {
      return (filt instanceof AcceptNoneFilter) || (
        (filt instanceof FrameFilter) &&
        ((FrameFilter) filt).getFrame() == getFrame()
      );
    }
  }

    /**
     * <code>ByPairwiseScore</code> is used to filter
     * <code>SimilarityPairFeature</code>s by their score. Features
     * are accepted if their score falls between the filter's minimum
     * and maximum values, inclusive. Features are rejected if they
     * are not <code>SimilarityPairFeature</code>s. The minimum value
     * accepted must be less than the maximum value.
     *
     * @author Keith James
     * @since 1.3
     */
    public static final class ByPairwiseScore implements OptimizableFilter {
        private double minScore;
        private double maxScore;
        private double score;
        private int    hashCode;

        /**
         * Creates a new <code>ByPairwiseScore</code>.
         *
         * @param minScore a <code>double</code>.
         * @param maxScore a <code>double</code>.
         */
        public ByPairwiseScore(double minScore, double maxScore) {
            if (minScore > maxScore)
                throw new IllegalArgumentException("Filter minimum score must be less than maximum score");

            this.minScore = minScore;
            this.maxScore = maxScore;

            hashCode += (minScore == 0.0 ? 0L : Double.doubleToLongBits(minScore));
            hashCode += (maxScore == 0.0 ? 0L : Double.doubleToLongBits(maxScore));
        }

        /**
         * Accept a Feature if it is an instance of
         * SimilarityPairFeature and its score is <= filter's minimum
         * score and >= filter's maximum score.
         *
         * @param f a <code>Feature</code>.
         * @return a <code>boolean</code>.
         */
        public boolean accept(Feature f) {
            if (! (f instanceof SimilarityPairFeature)) {
                return false;
            }

            score = ((SimilarityPairFeature) f).getScore();
            return (score >= minScore &&
                    score <= maxScore);
        }

        /**
         * <code>getMinScore</code> returns the minimum score
         * accepted.
         *
         * @return a <code>double</code>.
         */
        public double getMinScore() {
            return minScore;
        }

        /**
         * <code>getMaxScore</code> returns the maximum score
         * accepted.
         *
         * @return a <code>double</code>.
         */
        public double getMaxScore() {
            return maxScore;
        }

        public boolean equals(Object o) {
            if (o instanceof ByPairwiseScore) {
                ByPairwiseScore psf = (ByPairwiseScore) o;
                if (psf.getMinScore() == minScore &&
                    psf.getMaxScore() == maxScore) {
                    return true;
                }
            }
            return false;
        }

        public int hashCode() {
            return hashCode;
        }

        public boolean isProperSubset(FeatureFilter sup) {
            if (sup instanceof ByPairwiseScore) {
                ByPairwiseScore psf = (ByPairwiseScore) sup;
                return (psf.getMinScore() <= minScore &&
                        psf.getMaxScore() >= maxScore);
            }
            return false;
        }

        public boolean isDisjoint(FeatureFilter filt) {
            if (filt instanceof AcceptNoneFilter)
                return true;

            if (filt instanceof ByPairwiseScore) {
                ByPairwiseScore psf = (ByPairwiseScore) filt;
                return (psf.getMaxScore() < minScore ||
                        psf.getMinScore() > maxScore);
            }
            return false;
        }

        public String toString() {
            return minScore + " >= score <= " + maxScore;
        }
    }

    /**
     * Accepts features which are ComponentFeatures and have a <code>componentSequenceName</code>
     * property of the specified value.
     *
     * @author Thomas Down
     * @since 1.3
     */

    public final static class ByComponentName implements OptimizableFilter {
        private String cname;

        public ByComponentName(String cname) {
            this.cname = cname;
        }

        public boolean accept(Feature f) {
            if (f instanceof ComponentFeature) {
                return cname.equals(((ComponentFeature) f).getComponentSequenceName());
            } else {
                return false;
            }
        }

        public String getComponentName() {
            return cname;
        }

        public boolean equals(Object o) {
            return (o instanceof ByComponentName) && ((ByComponentName) o).getComponentName().equals(cname);
        }

        public int hashCode() {
            return getComponentName().hashCode();
        }

        public boolean isProperSubset(FeatureFilter sup) {
            if (sup instanceof ByComponentName) {
                return equals(sup);
            } else if (sup instanceof ByClass) {
                return ((ByClass) sup).getTestClass().isAssignableFrom(ComponentFeature.class);
            } else {
                return (sup instanceof AcceptAllFilter);
            }
        }

        public boolean isDisjoint(FeatureFilter feat) {
            if (feat instanceof ByComponentName) {
                return !equals(feat);
            } else if (feat instanceof ByClass) {
                Class featC = ((ByClass) feat).getTestClass();
                return ! (featC.isAssignableFrom(ComponentFeature.class));
            } else {
                return (feat instanceof AcceptNoneFilter);
            }
        }

        public String toString() {
            return "ByComponentName(" + cname + ")";
        }
    }

    /**
     * A filter which accepts only top-level Features.  This is true
     * is <code>getParent()</code> returns a <code>Sequence</code> instance.
     *
     * @since 1.3
     */

    public static final FeatureFilter top_level = new IsTopLevel();

    /**
     * A filter which accepts features with no children
     *
     * @since 1.3
     */

    // public static final FeatureFilter leaf = new IsLeaf();
    public static final FeatureFilter leaf = new FeatureFilter.OnlyChildren(FeatureFilter.none);

    // Note: this implements OptimizableFilter, but cheats :-).  Consequently,
    // other optimizablefilters don't know anything about it.  The convenience
    // methods on FilterUtils give ByFeature a higher precedence to make
    // sure this works out.

    /**
     * Accept only features which are equal to the specified feature
     *
     * @author Thomas Down
     * @since 1.3
     */

    public static final class ByFeature implements OptimizableFilter {
        private final Feature feature;

        public ByFeature(Feature f) {
            this.feature = f;
        }

        public Feature getFeature() {
            return feature;
        }

        public boolean accept(Feature f) {
            return f.equals(feature);
        }

        public boolean isProperSubset(FeatureFilter ff) {
            return ff.accept(feature);
        }

        public boolean isDisjoint(FeatureFilter ff) {
            return !ff.accept(feature);
        }

        public int hashCode() {
            return feature.hashCode() + 65;
        }

        public boolean equals(Object o) {
            if (o instanceof FeatureFilter.ByFeature) {
                return ((FeatureFilter.ByFeature) o).getFeature().equals(feature);
            } else {
                return false;
            }
        }
    }
}

