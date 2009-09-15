package org.biojava.bio.seq.filter;

import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.utils.walker.Visitor;

/**
 * Base-class for visitors that re-write a filter tree.
 *
 * <p>
 * This filter transformer will just duplicate a tree, using the same leaf
 * instances, and re-creating all logical filters, like And and ByDescendant.
 * </p>
 *
 * @author Matthew Pocock
 */
public class FilterTransformer
implements Visitor {
  public FeatureFilter featureFilter(FeatureFilter filter) {
    return filter;
  }

  public FeatureFilter and(FeatureFilter.And and,
                           FeatureFilter c1,
                           FeatureFilter c2)
  {
    return FilterUtils.and(c1, c2);
  }

  public FeatureFilter or(FeatureFilter.Or or,
                          FeatureFilter c1,
                          FeatureFilter c2)
  {
    return FilterUtils.or(c1, c2);
  }

  public FeatureFilter not(FeatureFilter.Not not, FeatureFilter c) {
    return FilterUtils.not(c);
  }

  public FeatureFilter byParent(FeatureFilter.ByParent parent, FeatureFilter c) {
    return FilterUtils.byParent(c);
  }

  public FeatureFilter byAncestor(FeatureFilter.ByAncestor ancestor, FeatureFilter c) {
    return FilterUtils.byAncestor(c);
  }

  public FeatureFilter onlyChildren(FeatureFilter.OnlyChildren child, FeatureFilter c) {
    return FilterUtils.onlyChildren(c);
  }

  public FeatureFilter onlyDescendants(FeatureFilter.OnlyDescendants desc, FeatureFilter c) {
    return FilterUtils.onlyDescendants(c);
  }

  public FeatureFilter byChild(FeatureFilter.ByChild child, FeatureFilter c) {
    return FilterUtils.byChild(c);
  }

  public FeatureFilter byDescendant(FeatureFilter.ByDescendant desc, FeatureFilter c) {
    return FilterUtils.byDescendant(c);
  }
}
