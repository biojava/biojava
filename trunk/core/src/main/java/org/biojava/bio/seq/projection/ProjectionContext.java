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

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Interface that defines the projection between original features and
 * projected features.
 *
 * <p>
 * <em>Note:</em> Implementing this interface is not trivial. It is assumed that
 * if you are implementing this, you are a power user or developer. All method
 * documentation is aimed at those groups, not users.
 * </p>
 *
 * <h2>Core projection methods</h2>
 *
 * <p>
 * This interface directly specifies a core set of methods that a
 * ProjectionContext implementation must provide. These are needed to allow
 * ProjectedFeatureHolder to connect to the context and for projected features
 * to be implemented on top of the context.
 * </p>
 *
 * <p>
 * All of these methods are implemented by @link ReparentContext, and it would
 * be unwise to try to implement them again independently, but it's a free
 * world and this is, after all, an interface. Some of the methods have
 * contracts that are a little tricky to implement without some thought.
 * </p>
 *
 * <p>
 * The methods projectFeature() and revertFeature() should be used as your
 * exclusive means for converting features between the two views. If you have
 * some utility method that produces projections, it /must/ be exclusively
 * invoked from projectFeature(). Likewise, although it is tempting to cast
 * the projected feature to Projection, and then to call getProjectedFeature(),
 * this will break the encapsulation of the context. The only method that should
 * do this cast is revertFeature().
 *
 * <h2>Extra projection methods</h2>
 *
 * <p>Projection methods are used by the projection engine to map projected
 * and unprojected feature properties. They are not declared directly in this
 * interface, but should be supplied by implementations. They will be picked
 * up by the projection engine using introspection, and connected to projected
 * feature instances by run-time code generation.
 * </p>
 *
 * <p>So that the projection methods can be found, they should all follow the
 * same template. If a feature interface has a property <code>foo</code> of
 * type <code>Bar</code>, then the feature interface will have the two methods:
 * <code><pre>
 * public void setFoo(Bar bar)
 * public Bar getFoo()
 * </pre></code>
 * The projection engine will be looking for a pair of methods in the context
 * named:
 * <code>
 * public Bar projectFoo(Bar bar);
 * public Bar revertFoo(Bar bar);
 * </pre></code>
 * If these methods are found, then the projected feature will have get/set
 * methods that resemble:
 * <code><pre>
 * public void setFoo(Bar bar) {
 *   getViewedFeature().setFoo(getContext().revertFoo(bar));
 * }
 *
 * public Bar getFoo() {
 *   return getContext().projectFoo(getViewedFeature().getFoo());
 * }
 * </pre></code>
 * If these methods are not found, then the projected feature will have get/set
 * methods that resembl:
 * <code><pre>
 * public void setFoo(Bar bar) {
 *   getViewedFeature().setFoo(bar);
 * }
 *
 * public Bar getFoo() {
 *   return getContext().getViewedFeature().getFoo();
 * }
 * </pre></code>
 * </p>
 *
 * <p>
 * Only those methods defined by the interface of the unprojected feature
 * will be mapped accross. So, if the context provides projectors for the
 * strand property but the feature is not a stranded feature, it's projection
 * will not have a strand property.
 * </p>
 *
 * You should probably not be implementing these yourself. Use the standard
 * factory methods in SequenceTools to create new sequences that are altered
 * views of other sequences.
 *
 * You will probably want to instantiate @link ReparentContext or @link
 * TranslateFlipContext to achieve the most commonly needed transformations.
 * If you do have to implement this, extend one of these two classes.
 *
 * Consider extending @link ReparentContext or @link TranslateFlipContex. They
 * do a lot of complex work for you.
 *
 * 
 * When projecting an original location <code>origLoc</code> to your new
 * location <code>yourLoc</code>, be sure to make sure the decorations match.
 * The easiest way to do this is return
 * <code>origLoc.newInstance(yourLoc)</code>. 
 *
 * 
 * Every ProjectionContext implementation must be a public or package-private
 * class. This is
 * because ProjectionEngine wires method calls from the projected features back
 * into the context. These methods are not necessarily defined by the context
 * interface (if they are project/revert pairs), so the class itself must be
 * directly accessible. The ProjectionEngine creates accessory classes in the
 * same package as the context it is using.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.4
 */

public interface ProjectionContext {
  /**
   * Get the features before projection.
   *
   * <p>
   * If you are projecting all the features in some feature holder, that is what
   * this method should return.
   * </p>
   *
   * @return the features before projection
   */
  public FeatureHolder getUnprojectedFeatures();

  public FeatureHolder getParent(Feature projFeat);

  /**
   * Get the sequence for a feature.
   *
   * <p>
   * This will be the return value of <code>projFeat.getParent()</code>.
   * </p>
   *
   * @param projFeat the projected Feature
   * @return the Sequence of the Feature
   */
  public Sequence getSequence(Feature projFeat);

  /**
   * Project a single feature.
   *
   * @param feat the Feature to project
   * @return a Feature representing feat after being transformed by this context
   */
  public Feature projectFeature(Feature feat);

  /**
   * Project all of the features in a FeatureHolder.
   *
   * <p>
   * <em>Warning:</em> The results of calling this method for features that are
   * not in getUnprojectedFeatures() is not specified by this API, but it is
   * reasonable to assume that bad things will happen.
   * </p>
   *
   * @param  features  the FeatureHolder containing the features to project
   * @return a FeatureHolder containing all the features projected
   */
  public FeatureHolder projectFeatures(FeatureHolder features);

  /**
   * Unproject a feature.
   *
   * <p>
   * This is the inverse opperation to @link projectFeature().
   * </p>
   *
   * <p>
   * <em>Note: </em> The result of calling this method for a feature that is
   * not projected through this context is not specified by this API, but it
   * is reasonable to assume that bad things will happen.
   * </p>
   *
   * @param  projFeat  the Feature to un-project
   * @return the unprojected feature
   */
  public Feature revertFeature(Feature projFeat);

  /**
   * Transform a filter on unprojected features so that it applies to projected
   * features.
   *
   * @param filt  the FeatureFilter to transform
   * @return the transformed FeatureFilter
   */
  public FeatureFilter projectFilter(FeatureFilter filt);

  /**
   * Transform a filter on projected features so that it applies to unprojected
   * features.
   *
   * @param filt  the FeatureFilter to transform
   * @return the transformed FeatureFilter
   */
  public FeatureFilter revertFilter(FeatureFilter filt);

  /**
   * Project all features that are children of feature so that they become
   * children of parent.
   *
   * @param feature  the Feature to project all children of
   * @param parent   the new parent feature holder
   * @return  a FeatureHolder containing projections of all children of feature
   *   so that f.getParent() is equal to parent
   */
  public FeatureHolder projectChildFeatures(Feature feature, FeatureHolder parent);

  /**
   * Create a projected feature with properties matching the template.
   *
   * <p>
   * You will probably implement this by delegating to the unprojected feature
   * holder. It is imperative that the template properties are unprojected first
   * so that when the newly created feature is projected, the properties match
   * up.
   * <p>
   * Not every projection context has fully reversible semantics. Use your
   * discression and come up with a reasonable plan that causes least supprise
   * to the user.
   * </p>
   *
   * @param projTempl the Feature.Template to instantiate
   * @return a new projected Feature matching the template as closely as possible
   * @throws BioException  if there was a problem instantiating the template
   * @throws ChangeVetoException  if the feature creation was vetoed
   */
  public Feature createFeature(Feature.Template projTempl)
      throws BioException, ChangeVetoException;

  /**
   * Remove the dying child.
   *
   * @param dyingChild  a projected feature to remove
   * @throws BioException if there is an error removing the feature
   * @throws ChangeVetoException if the removal of the feature was vetoed
   */
  public void removeFeature(Feature dyingChild)
      throws BioException, ChangeVetoException;

  /**
   * Create a new projected feature.
   *
   * <p>See the notes for @link createFeature(Feature.Template) for
   * implementation advice.
   * </p>
   *
   * @param projParent  the parent for the newly created feature
   * @param projTempl  the Feature.Template specifying the new feature
   * @return a new ProjectedFeature that is a child of projParent
   * @throws BioException if there was a problem creating the feature
   * @throws ChangeVetoException if the creation of the feature was vetoed
   */
  public Feature createFeature(Feature projParent, Feature.Template projTempl)
      throws BioException, ChangeVetoException;

  /**
   * Remove the dying child.
   *
   * @param projParent  the projected parent Feature
   * @param dyingChild  the child Feature to remove
   * @throws ChangeVetoException if the removal of the feature was vetoed
   */
  public void removeFeature(Feature projParent, Feature dyingChild)
      throws ChangeVetoException, BioException; // fixme: should we be throwing BioException?

  /**
   * Add a ChangeListener to a projected feature.
   *
   * @param projFeat the projected Feature to add the listener for
   * @param cl  the ChangeListener to add
   * @param ct  the ChangeType to register it for
   */
  public void addChangeListener(Feature projFeat, ChangeListener cl, ChangeType ct);

  /**
   * Remove a ChangeListener from a projected feature.
   *
   * @param projFeat  the projected Feature to remove the listener for
   * @param cl  the ChangeListener to remove
   * @param ct  the ChangeType it is registered for
   */
  public void removeChangeListener(Feature projFeat, ChangeListener cl, ChangeType ct);
}

