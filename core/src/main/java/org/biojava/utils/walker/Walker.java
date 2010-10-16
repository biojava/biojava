package org.biojava.utils.walker;


/**
 * Objects that can walk over a filter expression, showing each element to a
 * visitor.
 *
 * <p>
 * Walker implementations are not guaranteed to be thread-safe. In particular,
 * it is not possible to use the same Walker instance with more than one
 * thread if the visitor has return values. Walker implementations can be
 * re-used once the previous walk has been completed.
 * </p>
 *
 * 
 * You should use FilterUtils.visitFilter to apply a visitor to a feature
 * filter.
 *
 * 
 * You can use WalkerFactory.getInstance().getWalker(visitor) to get a walker
 * that is suitable for your visitor implementation. This will take care of
 * all the magic needed to hook up visitor call-back methods to the walkers
 * traversal of the features.
 *
 * If you don't like the walkers that WalkerFactory produces, you can implement
 * this directly. This will work fine for simple visitors, e.g., which only have
 * a single method for visting all filters, regardless of type.
 *
 * @author Matthew Pocock
 */
public interface Walker {
  /**
   * This walks the feature tree, showing the visitor each filter in
   * the expression.
   *
   * @param filter
   * @param visitor
   */
  public void walk(Object filter, Visitor visitor);

  /**
   * If the visitor has a return value, then the result of applying the visitor
   * to the tree can be obtained using this method, otherwise the result will
   * be null.
   *
   * @return the visitor's return value, or null
   */
  public Object getValue();
}
