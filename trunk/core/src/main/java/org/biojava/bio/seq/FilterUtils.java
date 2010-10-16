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
import java.util.HashSet;
import java.util.List;

import org.biojava.bio.AnnotationTools;
import org.biojava.bio.AnnotationType;
import org.biojava.bio.BioException;
import org.biojava.bio.CardinalityConstraint;
import org.biojava.bio.PropertyConstraint;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.utils.walker.Visitor;
import org.biojava.utils.walker.Walker;
import org.biojava.utils.walker.WalkerFactory;

/**
 * A set of FeatureFilter algebraic operations.
 *
 * @since 1.2
 * @author Matthew Pocock
 * @author Thomas Down
 */

public class FilterUtils {
  
    //that block avoids problems in WalkerFactory.generateWalker(Class visitorClass),
    //if none of these filters are instanciated before invoking that method.
    static {
      new FeatureFilter.And(null,null);
      new FeatureFilter.Or(null,null);
      new FeatureFilter.Not(null);
    } 
  
    private FilterUtils() {
    }
    
    /**
     * Determines if the set of features matched by sub can be <code>proven</code> to be a
     * proper subset of the features matched by sup.
     * <p>
     * If the filter sub matches only features that are matched by sup, then it is
     * a proper subset. It is still a proper subset if it does not match every
     * feature in sup, as long as no feature matches sub that is rejected by sup.
     * </p>
     *
     * @param sub the subset filter
     * @param sup the superset filter
     * @return <code>true</code> if <code>sub</code> is a proper subset of <code>sup</code>
     */

    public static boolean areProperSubset(FeatureFilter sub, FeatureFilter sup) {
      // Preconditions

      if (sub == null) {
        throw new NullPointerException("Null FeatureFilter: sub");
      }
      if (sup == null) {
        throw new NullPointerException("Null FeatureFilter: sup");
      }

      // Body

      if(sub.equals(sup)) {
        return true;
      } else if (sup == all()) {
        return true;
      } else if (sub == none()) {
        return true;
      } else if (sup instanceof FeatureFilter.And) {
        FeatureFilter.And and_sup = (FeatureFilter.And) sup;
        return areProperSubset(sub, and_sup.getChild1()) && areProperSubset(sub, and_sup.getChild2());
      } else if (sub instanceof FeatureFilter.And) {
        FeatureFilter.And and_sub = (FeatureFilter.And) sub;
        return areProperSubset(and_sub.getChild1(), sup) || areProperSubset(and_sub.getChild2(), sup);
      } else if (sub instanceof FeatureFilter.Or) {
        FeatureFilter.Or or_sub = (FeatureFilter.Or) sub;
        return areProperSubset(or_sub.getChild1(), sup) && areProperSubset(or_sub.getChild2(), sup);
      } else if (sup instanceof FeatureFilter.Or) {
        FeatureFilter.Or or_sup = (FeatureFilter.Or) sup;
        return areProperSubset(sub, or_sup.getChild1()) || areProperSubset(sub, or_sup.getChild2());
      } else if (sup instanceof FeatureFilter.Not) {
        FeatureFilter not_sup = ((FeatureFilter.Not) sup).getChild();
        return areDisjoint(sub, not_sup);
      } else if (sub instanceof FeatureFilter.Not) {
        // How do we prove this one?
        return false;  // return false for now
      } else if (sub instanceof OptimizableFilter) {
          // The works okay for ByFeature, too.
        return ((OptimizableFilter) sub).isProperSubset(sup);
      } else if(sup.equals(sub)) {
        return true;
      } else {
        return false;
      }
    }

    /**
     * Determines if two queries can be proven to be disjoint.
     * <p>
     * They are disjoint if there is no element that is matched by both filters
     * - that is, they have an empty intersection.  Order of arguments to this
     * method is not significant.
     * </p>
     *
     * @param a   the first FeatureFilter
     * @param b   the second FeatureFilter
     * @return <code>true</code> if they are proved to be disjoint, <code>false</code> otherwise
     */

    public static boolean areDisjoint(FeatureFilter a, FeatureFilter b) {
        // System.err.println("Disjunction test of " + a + " | " + b);

      // Preconditions

      if (a == null) {
        throw new NullPointerException("Null FeatureFilter: a");
      }
      if (b == null) {
        throw new NullPointerException("Null FeatureFilter: b");
      }

      // Body

      if(a.equals(b)) {
        return false;
      }

      if(a == none() || b == none()) {
        return true;
      }

      if (a == all()) {
        return areProperSubset(b, FeatureFilter.none);
      } else if(b == all()) {
        return areProperSubset(a, FeatureFilter.none);
      } else if (a == none() || b == none()) {
        return true;
      } if (a instanceof FeatureFilter.And) {
        FeatureFilter.And and_a = (FeatureFilter.And) a;
        return areDisjoint(and_a.getChild1(), b) || areDisjoint(and_a.getChild2(), b);
      } else if (b instanceof FeatureFilter.And) {
        FeatureFilter.And and_b = (FeatureFilter.And) b;
        return areDisjoint(a, and_b.getChild1()) || areDisjoint(a, and_b.getChild2());
      } else if (a instanceof FeatureFilter.Or) {
        FeatureFilter.Or or_a = (FeatureFilter.Or) a;
        return areDisjoint(or_a.getChild1(), b) && areDisjoint(or_a.getChild2(), b);
      } else if (b instanceof FeatureFilter.Or) {
        FeatureFilter.Or or_b = (FeatureFilter.Or) b;
        return areDisjoint(a, or_b.getChild1()) && areDisjoint(a, or_b.getChild2());
      } else if (a instanceof FeatureFilter.Not) {
        FeatureFilter not_a = ((FeatureFilter.Not) a).getChild();
        return areProperSubset(b, not_a);
      } else if (b instanceof FeatureFilter.Not) {
        FeatureFilter not_b = ((FeatureFilter.Not) b).getChild();
        return areProperSubset(a, not_b);
      } else if (a instanceof FeatureFilter.ByFeature) {
          return ((FeatureFilter.ByFeature) a).isDisjoint(b);
      } else if (b instanceof FeatureFilter.ByFeature) {
          return ((FeatureFilter.ByFeature) b).isDisjoint(a);
      } else if (a instanceof FeatureFilter.ByAncestor) {
        return ((OptimizableFilter) a).isDisjoint(b);
      } else if (b instanceof FeatureFilter.ByAncestor) {
        return ((OptimizableFilter) b).isDisjoint(a);
      } else if (a instanceof FeatureFilter.ByParent) {
        return ((OptimizableFilter) a).isDisjoint(b);
      } else if (b instanceof FeatureFilter.ByParent) {
        return ((OptimizableFilter) b).isDisjoint(a);
      } else if (a instanceof OptimizableFilter) {
        return ((OptimizableFilter) a).isDisjoint(b);
      } else if (b instanceof OptimizableFilter) {
        return ((OptimizableFilter) b).isDisjoint(a);
      }

      // *SIGH* we don't have a proof here...

      return false;
    }

    /**
     * Try to determine the minimal location which all features matching a given
     * filter must overlap.
     *
     * @param ff A feature filter
     * @return the minimal location which any features matching <code>ff</code>
     *          must overlap, or <code>null</code> if no proof is possible
     *          (normally indicates that the filter has nothing to do with
     *          location).
     * @since 1.2
     */

    public static Location extractOverlappingLocation(FeatureFilter ff) {
	if (ff instanceof FeatureFilter.OverlapsLocation) {
	    return ((FeatureFilter.OverlapsLocation) ff).getLocation();
	} else if (ff instanceof FeatureFilter.ContainedByLocation) {
	    return ((FeatureFilter.ContainedByLocation) ff).getLocation();
	} else if (ff instanceof FeatureFilter.And) {
	    FeatureFilter.And ffa = (FeatureFilter.And) ff;
	    Location l1 = extractOverlappingLocation(ffa.getChild1());
	    Location l2 = extractOverlappingLocation(ffa.getChild2());

	    if (l1 != null) {
		if (l2 != null) {
		    return l1.intersection(l2);
		} else {
		    return l1;
		}
	    } else {
		if (l2 != null) {
		    return l2;
		} else {
		    return null;
		}
	    }
	} else if (ff instanceof FeatureFilter.Or) {
	    FeatureFilter.Or ffo = (FeatureFilter.Or) ff;
	    Location l1 = extractOverlappingLocation(ffo.getChild1());
	    Location l2 = extractOverlappingLocation(ffo.getChild2());

	    if (l1 != null && l2 != null) {
		return LocationTools.union(l1, l2);
	    }
	}

	return null;
    }

  /**
   * Decide if two feature filters accept exactly the same set of features.
   *
   * <p>
   * Two feature filters are equal if it can be proven that
   * <code>f1.accept(feature) == f2.accept(feature)</code> for all values of
   * <code>feature</code>. If areEqual returns false, this may indicate that
   * they accept clearly different sets of features. It may also, howerver,
   * indicate that the method was unable to prove that they were equal.
   * </p>
   *
   * <p>
   * Note that given a finite set of features, f1 and f2 may match the same
   * sub-set of those features even if they are not equal.
   * </p>
   *
   * @param f1  the first filter
   * @param f2  the second filter
   * @return true if they can be proven to be equivalent
   */
  public final static boolean areEqual(FeatureFilter f1, FeatureFilter f2) {
    if(f1 instanceof FeatureFilter.And && f2 instanceof FeatureFilter.And) {
      List f1f = new ArrayList();
      List f2f = new ArrayList();

      expandAnd(f1, f1f);
      expandAnd(f2, f2f);

      return new HashSet(f1f).equals(new HashSet(f2f));
    } else if(f1 instanceof FeatureFilter.And || f2 instanceof FeatureFilter.And) {
      return false;
    } else if(f1 instanceof FeatureFilter.Or && f2 instanceof FeatureFilter.Or) {
      List f1f = new ArrayList();
      List f2f = new ArrayList();

      expandOr(f1, f1f);
      expandOr(f2, f2f);

      return new HashSet(f1f).equals(new HashSet(f2f));
    } else if(f1 instanceof FeatureFilter.Or || f2 instanceof FeatureFilter.Or) {
      return false;
    } else {
      return f1.equals(f2);
    }
  }

  /**
   * Construct a filter which matches features with a specific <code>type</code>
   * value.
   */
  
  public final static FeatureFilter byType(String type) {
    return new FeatureFilter.ByType(type);
  }

  /**
   * Construct a filter which matches features with a specific <code>source</code>
   * value.
   */
  
  public final static FeatureFilter bySource(String source) {
    return new FeatureFilter.BySource(source);
  }

  /**
   * Construct a filter which matches features which are assignable to the
   * specified class or interface.
   */
  
  public final static FeatureFilter byClass(Class clazz)
  throws ClassCastException {
    return new FeatureFilter.ByClass(clazz);
  }

  /** 
   * Construct a filter which matches features with locations wholly contained
   * by the specified <code>Location</code>.
   */
  
  public final static FeatureFilter containedByLocation(Location loc) {
    return new FeatureFilter.ContainedByLocation(loc);
  }

  /** 
   * Construct a filter which matches features with locations contained by or
   * overlapping the specified <code>Location</code>.
   */
  
  public final static FeatureFilter overlapsLocation(Location loc) {
    return new FeatureFilter.OverlapsLocation(loc);
  }

  /** 
   * Construct a filter which matches features with locations where the interval
   * between the <code>min</code> and <code>max</code> positions are contained by or
   * overlap the specified <code>Location</code>.
   */
  
  public final static FeatureFilter shadowOverlapsLocation(Location loc) {
      return new FeatureFilter.ShadowOverlapsLocation(loc);
  }
  
  /** 
   * Construct a filter which matches features with locations where the interval
   * between the <code>min</code> and <code>max</code> positions are contained by 
   * the specified <code>Location</code>.
   */
  
  public final static FeatureFilter shadowContainedByLocation(Location loc) {
      return new FeatureFilter.ShadowContainedByLocation(loc);
  }

  /**
   * Match features attached to sequences with a specified name.
   */
  
  public final static FeatureFilter bySequenceName(String name) {
    return new FeatureFilter.BySequenceName(name);
  }

  /**
   * Construct a new filter which is the negation of <code>filter</code>.
   */
  
  public final static FeatureFilter not(FeatureFilter filter) {
    return new FeatureFilter.Not(filter);
  }

  /**
   * Construct a new filter which matches the intersection of two other
   * filters.
   */
  
  public final static FeatureFilter and(FeatureFilter c1, FeatureFilter c2) {
    return new FeatureFilter.And(c1, c2);
  }

  /**
   * Constructs a new filter which matches the intersection of a set of filters.
   */
  
  public final static FeatureFilter and(FeatureFilter[] filters) {
    if(filters.length == 0) {
      return all();
    } else if(filters.length == 1) {
      return filters[0];
    } else {
      FeatureFilter f = and(filters[0], filters[1]);
      for(int i = 2; i < filters.length; i++) {
        f = and(f, filters[i]);
      }
      return f;
    }
  }

  /**
   * Construct a new filter which matches the union of two filters.
   */
  
  public final static FeatureFilter or(FeatureFilter c1, FeatureFilter c2) {
    return new FeatureFilter.Or(c1, c2);
  }

  /**
   * Construct a new filter which matches the intersection of two filters.
   */
  
  public final static FeatureFilter or(FeatureFilter[] filters) {
    if(filters.length == 0) {
      return none();
    } else if(filters.length == 1) {
      return filters[0];
    } else {
      FeatureFilter f = or(filters[0], filters[1]);
      for(int i = 2; i < filters.length; i++) {
        f = or(f, filters[i]);
      }
      return f;
    }
  }

  /**
   * Match features with annotations matching the specified <code>AnnotationType</code>
   */
  
  public final static FeatureFilter byAnnotationType(AnnotationType type) {
    return new FeatureFilter.ByAnnotationType(type);
  }

  /**
   * Match features where the annotation property named <code>key</code> is
   * equal to <code>value</code>.
   */
  
  public final static FeatureFilter byAnnotation(Object key, Object value) {
    return new FeatureFilter.ByAnnotation(key, value);
  }

  /**
   * Match features where the annotation property named <code>key</code> is
   * an instance of <code>valClass</code>.
   */
  
  public final static FeatureFilter byAnnotationType(Object key, Class valClass) {
    AnnotationType type = new AnnotationType.Impl();
    type.setConstraints(key, new PropertyConstraint.ByClass(valClass), CardinalityConstraint.ANY);
    return byAnnotationType(type);
  }

  /**
   * Match features where the property <code>key</code> has been defined as having
   * some value, regardless of the exact value.
   */
  
  public final static FeatureFilter hasAnnotation(Object key) {
    return new FeatureFilter.HasAnnotation(key);
  }

  /**
   * Match StrandedFeatures on the specified strand.
   */
  
  public final static FeatureFilter byStrand(StrandedFeature.Strand strand) {
    return new FeatureFilter.StrandFilter(strand);
  }

  /**
   * Match features where the parent feature matches the specified filter.
   * This cannot match top-level features.
   */
  
  public final static FeatureFilter byParent(FeatureFilter parentFilter) {
    return new FeatureFilter.ByParent(parentFilter);
  }

  /**
   * Match features where at least one of the ancestors matches the specified
   * filter.  This cannot match top-level features.
   */
  
  public final static FeatureFilter byAncestor(FeatureFilter ancestorFilter) {
    return new FeatureFilter.ByAncestor(ancestorFilter);
  }

  /**
   * Match features where at least one child feature matches the supplied filter.
   * This does not match leafFeatures.
   */
  
  public final static FeatureFilter byChild(FeatureFilter childFilter) {
    return new FeatureFilter.ByChild(childFilter);
  }

  /**
   * Match features where at least one decendant feature -- possibly but not necessarily an
   * immediate child -- matches the specified filter.
   */
  
  public final static FeatureFilter byDescendant(FeatureFilter descFilter) {
    return new FeatureFilter.ByDescendant(descFilter);
  }

  /**
   * Construct a filter which matches features whose children all match the
   * specified filter.  This filter always matches leaf features.
   */
  
  public final static FeatureFilter onlyChildren(FeatureFilter child) {
    return new FeatureFilter.OnlyChildren(child);
  }

  /**
   * Construct a filter which matches features whose decendants all match the
   * specified filter.  This filter always matches leaf features.
   */
  
  public final static FeatureFilter onlyDescendants(FeatureFilter desc) {
    return new FeatureFilter.OnlyDescendants(desc);
  }

  /**
   * Construct a filter which matches FramedFeatures with the specified reading
   * frame.
   */
  
  public final static FeatureFilter byFrame(FramedFeature.ReadingFrame frame) {
    return new FeatureFilter.FrameFilter(frame);
  }

  /**
   * Match SeqSimilaritiy features with scores in the specified range.
   */
  
  public final static FeatureFilter byPairwiseScore(double minScore, double maxScore) {
    return new FeatureFilter.ByPairwiseScore(minScore, maxScore);
  }

  /**
   * Construct a filter which matches all features which implement the
   * <code>ComponentFeature</code> interface and have a <code>componentName</code>
   * property equal to the specified value
   */
  
  public final static FeatureFilter byComponentName(String compName) {
    return new FeatureFilter.ByComponentName(compName);
  }

  /**
   * Return a filter which matches all top-level features.  These are features
   * which are direct children of a <code>Sequence</code> rather than another
   * <code>Feature</code>.
   */
  
  public final static FeatureFilter topLevel() {
    return FeatureFilter.top_level;
  }

  /**
   * Return a filter which matches features with zero children.
   */
  
  public final static FeatureFilter leaf() {
    return FeatureFilter.leaf;
  }

  /**
   * Return a filter which matches all features.
   */
  
  public final static FeatureFilter all() {
    return FeatureFilter.all;
  }

  /**
   * Return a filter which matches no features.
   */
  
  public final static FeatureFilter none() {
    return FeatureFilter.none;
  }

  /**
   * Attempts to reduce a FeatureFilter to an equivalent FeatureFilter with
   * fewer terms.
   *
   * <p>
   * This will attempt to push all leaf constraints as far from the root of the
   * filter expression as possible, in an attept to prove an empty or universal
   * set. It will then propogate these through the logical operators in an
   * attempt to reduce the entire expression to the empty or universal set. If
   * filters can be combined (for example, overlapping constraints), then this
   * will happen on the way.
   * </p>
   *
   * <p>
   * The resulting filter is guaranteed to accept exactly the same set of\
   * features that is accepted by the argument. In particular,
   * <code>areEqual(filter, optimize(filter))</code> is always true.
   * </p>
   *
   * @param filter  the FeatureFilter to optimize
   * @return an optimized version
   */
  public final static FeatureFilter optimize(FeatureFilter filter) {
    //depth++;
    //System.out.println(depth + ":" + "Optimizing " + filter);
    try {
    if(filter instanceof FeatureFilter.And) {
      //System.out.println(depth + ":" + "is AND");

      List filters = new ArrayList();
      expandAnd(filter, filters);
      //System.out.println(depth + ":" + "as list: " + filters);

      // get all children of this AND, and of all AND children
      for(int i = 0; i < filters.size(); i++) {
        filters.set(i, optimize((FeatureFilter) filters.get(i)));
      }

      // now scan this list for all OR children
      List ors = new ArrayList();
      for(int i = 0; i < filters.size(); i++) {
        FeatureFilter f = (FeatureFilter) filters.get(i);
        if(f instanceof FeatureFilter.Or) {
          filters.remove(i);
          ors.add(f);
          i--;
        }
      }

      // optimize all simple filters
      for(int i = 1; i < filters.size(); i++) {
        //System.out.println(depth + ":" + "i: " + i);
        FeatureFilter a = (FeatureFilter) filters.get(i);
        for(int j = 0; j < i; j++) {
          //System.out.println(depth + ":" + "j: " + j);
          FeatureFilter b = (FeatureFilter) filters.get(j);

          //System.out.println(depth + ":" + "Comparing " + a + ", " + b + " of " + filters);

          if(areDisjoint(a, b)) {
            // a n b = E
            //System.out.println(depth + ":" + "Disjoint. Returning none()");
            return none();
          } else if(areProperSubset(a, b)) {
            //System.out.println(depth + ":" + "a < b. Removing b");
            // if a < b then a n b = a
            filters.remove(j);
            j--; i--;
          } else if(areProperSubset(b, a)) {
            //System.out.println(depth + ":" + "a > b. Removing a");
            // if a > b then a n b = b
            filters.remove(i);
            i--;
            break;
          } else {
            //System.out.println(depth + ":" + "Attempting to calculate intersection");
            FeatureFilter intersect = intersection(a, b);
            if(intersect != null) {
              //System.out.println(depth + ":" + "got intersection: " + intersect);
              filters.set(i, intersect);
              filters.remove(j);
              j--; i-=2;
            } else {
              //System.out.println(depth + ":" + "no luck - moving on");
            }
          }
        }
      }
      //System.out.println(depth + ":" + "Reduced to: " + filters);

      if(filters.isEmpty()) {
        if(ors.isEmpty()) {
          //System.out.println("All empty. Returning none()");
          return none();
        } else {
          FeatureFilter andedOrs = FilterUtils.and((FeatureFilter[]) filters.toArray(new FeatureFilter[] {}));
          //System.out.println("No filters, some ors, returning: " + andedOrs);
          return andedOrs;
        }
      } else {
        FeatureFilter ands = and((FeatureFilter[]) filters.toArray(new FeatureFilter[] {}));
        if(ors.isEmpty()) {
          //System.out.println(depth + ":" + "No ors, just returning anded values: " + ands);
          return ands;
        } else {
          //System.out.println(depth + ":" + "Mixing in ors: " + ors);
          //System.out.println(depth + ":" + "to: " + ands);
          List combs = new ArrayList();
          combs.add(ands);
          for(int i = 0; i < ors.size(); i++) {
            List newCombs = new ArrayList();
            FeatureFilter.Or or = (FeatureFilter.Or) ors.get(i);
            for(int j = 0; j < combs.size(); j++) {
              FeatureFilter f = (FeatureFilter) combs.get(j);
              newCombs.add(and(f, or.getChild1()));
              newCombs.add(and(f, or.getChild2()));
            }
            combs = newCombs;
            //System.out.println("Expanded To: " + combs);
          }
          FeatureFilter res = optimize(or((FeatureFilter[]) combs.toArray(new FeatureFilter[] {})));
          //System.out.println(depth + ":" + "Returning optimized or: " + res);
          return res;
        }
      }
    } else if(filter instanceof FeatureFilter.Or) {
      //System.out.println(depth + ":" + "is OR");

      List filters = new ArrayList();
      expandOr(filter, filters);

      //System.out.println(depth + ":" + "as list: " + filters);

      for(int i = 0; i < filters.size(); i++) {
        filters.set(i, optimize((FeatureFilter) filters.get(i)));
      }

      // now scan this list for all OR children
      List ands = new ArrayList();
      for(int i = 0; i < filters.size(); i++) {
        FeatureFilter f = (FeatureFilter) filters.get(i);
        if(f instanceof FeatureFilter.And) {
          filters.remove(i);
          ands.add(f);
          i--;
        }
      }
      //System.out.println(depth + ":" + "filters: " + filters);
      //System.out.println(depth + ":" + "ands: " + ands);

      for(int i = 1; i < filters.size(); i++) {
        //System.out.println(depth + ":" + "i: " + i);
        FeatureFilter a = (FeatureFilter) filters.get(i);
        for(int j = 0; j < i; j++) {
          //System.out.println(depth + ":" + "j: " + j);
          FeatureFilter b = (FeatureFilter) filters.get(j);

          //System.out.println(depth + ":" + "Comparing " + a + ", " + b + " of " + filters);

          if(a == all() || b == all()) {
            //System.out.println(depth + ":" + "Found an all. Returning all()");
            return all();
          } else if(areProperSubset(a, b)) {
            //System.out.println(depth + ":" + "a < b. Removing a");
            filters.remove(i);
            i--;
            break;
          } else if(areProperSubset(b, a)) {
            //System.out.println(depth + ":" + "a > b. Removing b");
            filters.remove(j);
            j--; i-=2;
          } else {
            //System.out.println(depth + ":" + "Trying to calculate union");
            FeatureFilter union = union(a, b);
            if(union != null) {
              //System.out.println(depth + ":" + "Got union: " + union);
              filters.set(i, union);
              filters.remove(j);
              j--; i--;
            } else {
              //System.out.println(depth + ":" + "no luck - moving on");
            }
          }
        }
      }
      //System.out.println(depth + ":" + "Reduced to: " + filters);

      if(filters.isEmpty()) {
        if(ands.isEmpty()) {
          //System.out.println("All empty. Returning none()");
          return none();
        } else {
          FeatureFilter oredAnds = or((FeatureFilter[]) ands.toArray(new FeatureFilter[] {}));
          //System.out.println("No filters, some ands. returning: " + oredAnds);
          return oredAnds;
        }
      } else {
        FeatureFilter ors = or((FeatureFilter[]) filters.toArray(new FeatureFilter[] {}));
        if(ands.isEmpty()) {
          //System.out.println(depth + ":" + "no ands, returning ors: " + ors);
          return ors;
        } else {
          //System.out.println(depth + ":" + "Mixing in ands: " + ands);
          //System.out.println(depth + ":" + "to: " + ors);
          List combs = new ArrayList();
          combs.add(ors);
          for(int i = 0; i < ands.size(); i++) {
            List newCombs = new ArrayList();
            FeatureFilter.And and = (FeatureFilter.And) ands.get(i);
            for(int j = 0; j < combs.size(); j++) {
              FeatureFilter f = (FeatureFilter) combs.get(j);
              newCombs.add(or(f, and.getChild1()));
              newCombs.add(or(f, and.getChild2()));
            }
            combs = newCombs;
          }
          FeatureFilter val = optimize(and((FeatureFilter[]) combs.toArray(new FeatureFilter[] {})));
          //System.out.println(depth + ":" + "returning anded values: " + val);
          return val;
        }
      }
    } else if(filter instanceof FeatureFilter.Not) {
      FeatureFilter.Not not = (FeatureFilter.Not) filter;
      return not(optimize(not.getChild()));
    } else {
      return filter;
    }
    } finally {
      // depth--;
    }
  }

  private static FeatureFilter intersection(FeatureFilter f1, FeatureFilter f2) {
    if(
      f1 instanceof FeatureFilter.ContainedByLocation &&
      f2 instanceof FeatureFilter.ContainedByLocation
    ) {
      Location loc = LocationTools.intersection(
        ((FeatureFilter.ContainedByLocation) f1).getLocation(),
        ((FeatureFilter.ContainedByLocation) f2).getLocation()
      );
      if(loc == Location.empty) {
        return none();
      } else {
        return containedByLocation(loc);
      }
    } else if(
      f1 instanceof FeatureFilter.OverlapsLocation &&
      f2 instanceof FeatureFilter.OverlapsLocation
    ) {
      // can't do much here
    } else if(
      f1 instanceof FeatureFilter.ByAnnotationType &&
      f2 instanceof FeatureFilter.ByAnnotationType
    ) {
      FeatureFilter.ByAnnotationType f1t = (FeatureFilter.ByAnnotationType) f1;
      FeatureFilter.ByAnnotationType f2t = (FeatureFilter.ByAnnotationType) f2;

      AnnotationType intersect = AnnotationTools.intersection(
        f1t.getType(),
        f2t.getType()
      );

      if(intersect == AnnotationType.NONE) {
        return none();
      } else {
        return byAnnotationType(intersect);
      }
    } else if(
      f1 instanceof ByHierarchy &&
      f2 instanceof ByHierarchy
    ) {
      ByHierarchy f1h = (ByHierarchy) f1;
      ByHierarchy f2h = (ByHierarchy) f2;

      if(f1 instanceof Up && f2 instanceof Up) {
        FeatureFilter filt = optimize(or(f1h.getFilter(), f2h.getFilter()));
        if(filt == none()) {
          return none();
        }

        if(
          f1h instanceof FeatureFilter.ByParent ||
          f2h instanceof FeatureFilter.ByParent
        ) {
          return byParent(filt);
        } else {
          return byAncestor(filt);
        }
      } else if(f1 instanceof Down && f2 instanceof Down) {
        FeatureFilter filt = optimize(or(f1h.getFilter(), f2h.getFilter()));
        if(filt == none()) {
          return none();
        }

        if(
          f1h instanceof FeatureFilter.ByChild ||
          f2h instanceof FeatureFilter.ByChild
        ) {
          return byChild(filt);
        } else {
          return byDescendant(filt);
        }
      } else {
        return none();
      }
    }

    return null;
  }

  private static FeatureFilter union(FeatureFilter f1, FeatureFilter f2) {
    if(
      f1 instanceof FeatureFilter.ContainedByLocation &&
      f2 instanceof FeatureFilter.ContainedByLocation
    ) {
      return containedByLocation(LocationTools.union(
        ((FeatureFilter.ContainedByLocation) f1).getLocation(),
        ((FeatureFilter.ContainedByLocation) f2).getLocation()
      ));
    } else if(
      f1 instanceof FeatureFilter.OverlapsLocation &&
      f2 instanceof FeatureFilter.OverlapsLocation
    ) {
      return overlapsLocation(LocationTools.intersection(
        ((FeatureFilter.OverlapsLocation) f1).getLocation(),
        ((FeatureFilter.OverlapsLocation) f2).getLocation()
      ));
    } else if(
      f1 instanceof FeatureFilter.ByAnnotationType &&
      f2 instanceof FeatureFilter.ByAnnotationType
    ) {
      FeatureFilter.ByAnnotationType f1t = (FeatureFilter.ByAnnotationType) f1;
      FeatureFilter.ByAnnotationType f2t = (FeatureFilter.ByAnnotationType) f2;

      AnnotationType union = AnnotationTools.union(
        f1t.getType(),
        f2t.getType()
      );

      return byAnnotationType(union);
    } else if(
      f1 instanceof ByHierarchy &&
      f2 instanceof ByHierarchy
    ) {
      ByHierarchy f1h = (ByHierarchy) f1;
      ByHierarchy f2h = (ByHierarchy) f2;

      if(f1 instanceof Up && f2 instanceof Up) {
        FeatureFilter filt = optimize(or(f1h.getFilter(), f2h.getFilter()));
        if(filt == none()) {
          return none();
        }

        if(
          f1h instanceof FeatureFilter.ByAncestor ||
          f2h instanceof FeatureFilter.ByAncestor
        ) {
          return byAncestor(filt);
        } else {
          return byParent(filt);
        }
      } else if(f1 instanceof Down && f2 instanceof Down) {
        FeatureFilter filt = optimize(or(f1h.getFilter(), f2h.getFilter()));
        if(filt == none()) {
          return none();
        }

        if(
          f1h instanceof FeatureFilter.ByDescendant ||
          f2h instanceof FeatureFilter.ByDescendant
        ) {
          return byDescendant(filt);
        } else {
          return byChild(filt);
        }
      } else {
        return none();
      }
    }

    return null;
  }

  private static void expandAnd(FeatureFilter filt, List filters) {
    if(filt instanceof FeatureFilter.And) {
      FeatureFilter.And and = (FeatureFilter.And) filt;
      expandAnd(and.getChild1(), filters);
      expandAnd(and.getChild2(), filters);
    } else {
      filters.add(filt);
    }
  }

  private static void expandOr(FeatureFilter filt, List filters) {
    if(filt instanceof FeatureFilter.Or) {
      FeatureFilter.Or or = (FeatureFilter.Or) filt;
      expandOr(or.getChild1(), filters);
      expandOr(or.getChild2(), filters);
    } else {
      filters.add(filt);
    }
  }

  /**
   * <p>This is a general framework method for transforming one filter into another.
   * This method will handle the logical elements of a query (and, or, not) and delegate
   * all the domain-specific munging to a FilterTransformer object.</p>
   *
   * <p>The transformer could flip strands and locations of elements of a filter, add or
   * remove attributes required in annotations, or systematically alter feature types
   * or sources.</p>
   *
   * @param ff  the FeatureFilter to transform
   * @param trans  a FilterTransformer encapsulating rules about how to transform filters
   */
  public static FeatureFilter transformFilter(FeatureFilter ff, FilterTransformer trans) {
    if(ff == null) {
      throw new NullPointerException("Can't transform null filters");
    }

    if(ff instanceof FeatureFilter.And) {
      FeatureFilter.And and = (FeatureFilter.And) ff;
      return and(
        transformFilter(and.getChild1(), trans),
        transformFilter(and.getChild2(), trans)
      );
    } else if(ff instanceof FeatureFilter.Or) {
      FeatureFilter.Or or = (FeatureFilter.Or) ff;
      return or(
        transformFilter(or.getChild1(), trans),
        transformFilter(or.getChild2(), trans)
      );
    } else if(ff instanceof FeatureFilter.Not) {
      return not(
        transformFilter(((FeatureFilter.Not) ff).getChild(), trans)
      );
    } else {
      FeatureFilter tf = trans.transform(ff);
      if(tf != null) {
        return tf;
      } else {
        return ff;
      }
    }
  }

  /**
   * An object able to transform some FeatureFilter instances sytematically into others.
   *
   * @author Matthew Pocock
   */
  public interface FilterTransformer {
    /**
     * Transform a filter, or return null if it can not be transformed.
     *
     * @param filter  the FeatureFilter to attempt to transform
     * @return a transformed filter, or null
     */
    public FeatureFilter transform(FeatureFilter filter);
  }

  /**
   * An implementation of FilterTransformer that attempts to transform by one transformer,
   * and if that fails, by another.
   *
   * @author Matthew Pocock
   */
  public class DelegatingTransformer
  implements FilterTransformer {
    FilterTransformer t1;
    FilterTransformer t2;

    /**
     * Create a new DelegatingTransformer that will apply t1 and then t2 if t1 fails.
     *
     * @param t1 the first FilterTransformer to try
     * @param t2 the seccond FilterTransformer to try
     */
    public DelegatingTransformer(FilterTransformer t1, FilterTransformer t2) {
      this.t1 = t1;
      this.t2 = t2;
    }

    public FeatureFilter transform(FeatureFilter ff) {
      FeatureFilter res = t1.transform(ff);
      if(res == null) {
        res = t2.transform(ff);
      }
      return res;
    }
  }

    static FeatureFilter getOnlyDescendantsFilter(FeatureFilter ff) {
        if (ff instanceof FeatureFilter.OnlyDescendants) {
            return ((FeatureFilter.OnlyDescendants) ff).getFilter();
        } else if (ff instanceof FeatureFilter.And) {
            FeatureFilter.And ffa = (FeatureFilter.And) ff;
            FeatureFilter ocf1 = getOnlyDescendantsFilter(ffa.getChild1());
            FeatureFilter ocf2 = getOnlyDescendantsFilter(ffa.getChild2());
            if (ocf1 == null) {
                return ocf2;
            } else if (ocf2 == null) {
                return ocf1;
            } else {
                return new FeatureFilter.And(ocf1, ocf2);
            }
        } else if (ff instanceof FeatureFilter.Or) {
            FeatureFilter.Or ffa = (FeatureFilter.Or) ff;
            FeatureFilter ocf1 = getOnlyDescendantsFilter(ffa.getChild1());
            FeatureFilter ocf2 = getOnlyDescendantsFilter(ffa.getChild2());
            if (ocf1 == null) {
                return ocf2;
            } else if (ocf2 == null) {
                return ocf1;
            } else {
                return new FeatureFilter.Or(ocf1, ocf2);
            }
        } else {
            return null;
        }
    }

    static FeatureFilter getOnlyChildrenFilter(FeatureFilter ff) {
        // System.err.println("In getOnlyChildrenFilter: " + ff.toString());

        if (ff instanceof FeatureFilter.OnlyChildren) {
            return ((FeatureFilter.OnlyChildren) ff).getFilter();
        } else if (ff instanceof FeatureFilter.OnlyDescendants) {
            return ((FeatureFilter.OnlyDescendants) ff).getFilter();
        } else if (ff instanceof FeatureFilter.And) {
            FeatureFilter.And ffa = (FeatureFilter.And) ff;
            FeatureFilter ocf1 = getOnlyChildrenFilter(ffa.getChild1());
            FeatureFilter ocf2 = getOnlyChildrenFilter(ffa.getChild2());
            if (ocf1 == null) {
                return ocf2;
            } else if (ocf2 == null) {
                return ocf1;
            } else {
                return new FeatureFilter.And(ocf1, ocf2);
            }
        } else if (ff instanceof FeatureFilter.Or) {
            FeatureFilter.Or ffa = (FeatureFilter.Or) ff;
            FeatureFilter ocf1 = getOnlyChildrenFilter(ffa.getChild1());
            FeatureFilter ocf2 = getOnlyChildrenFilter(ffa.getChild2());
            if (ocf1 == null) {
                return ocf2;
            } else if (ocf2 == null) {
                return ocf1;
            } else {
                return new FeatureFilter.Or(ocf1, ocf2);
            }
        } else {
            return null;
        }
    }

  /**
   * Applies a visitor to a filter, and returns the visitor's result or null.
   *
   * @param filter  the filter to scan
   * @param visitor the visitor to scan with
   * @return  the result of the visitor or null
   * @throws BioException   if the required walker could not be created
   */
  public static Object visitFilter(FeatureFilter filter, Visitor visitor)
  throws BioException {
    Walker walker = WalkerFactory.getInstance().getWalker(visitor);
    walker.walk(filter, visitor);
    return walker.getValue();
  }
}
