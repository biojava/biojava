package org.biojava.bio.program.tagvalue;

/**
 * Interface for objects that change tag names or properties systematically.
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public interface PropertyChanger {
  /**
   * <p>
   * <code>getNewTag</code> returns the tag which substitutes the
   * specified value.
   * </p>
   *
   * <p>
   * If there is no mapping associated with this tag, it is returned
   * unchanged.
   * </p>
   *
   * @param oldTag an <code>Object</code> to substitute.
   *
   * @return an <code>Object</code>.
   */
  Object getNewTag(Object oldTag);
}
