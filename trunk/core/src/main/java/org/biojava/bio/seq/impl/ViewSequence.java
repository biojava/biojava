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

package org.biojava.bio.seq.impl;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.OverlayAnnotation;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FeatureRealizer;
import org.biojava.bio.seq.MergeFeatureHolder;
import org.biojava.bio.seq.RealizingFeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.seq.projection.ProjectedFeatureHolder;
import org.biojava.bio.seq.projection.ReparentContext;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;

/**
 * A view onto another Sequence object.  This class allows new
 * features and annotations to be overlaid onto an existing
 * Sequence without modifying it.
 *
 * You will almost certainly want to be calling
 * {@link org.biojava.bio.seq.SequenceTools#view(Sequence seq)} instead of instantiating this
 * class directly.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

public class ViewSequence
  extends
        Unchangeable
  implements
        Sequence,
        RealizingFeatureHolder,
        Serializable
{
    private static final long serialVersionUID = 9866447;
    /**
     * Delegate Sequence.
     */
    private Sequence seqDelegate;

    /**
     * FeatureHolder support
     */
    private MergeFeatureHolder exposedFeatures;
    private SimpleFeatureHolder addedFeatures;

    /**
     * IDs
     */
    private String name;
    private String urn;

    /**
     * Our annotation.
     */
    private Annotation anno;

    /**
     * The FeatureRealizer we use.
     */
    private transient FeatureRealizer featureRealizer;

    private void readObject(ObjectInputStream s)throws IOException, ClassNotFoundException{
        s.defaultReadObject();
        this.featureRealizer = FeatureImpl.DEFAULT;
    }

  /**
   * Construct a view onto an existing sequence and give it a new
   * name.
   * <p>The prefered method is SequenceTools.view(Sequence seq, String name)
   */
  public ViewSequence(Sequence seq, String name) {
    this.name = name;

    seqDelegate = seq;
    addedFeatures = new SimpleFeatureHolder();
    exposedFeatures = new MergeFeatureHolder();
    try {
      exposedFeatures.addFeatureHolder(new ProjectedFeatureHolder(
              new ReparentContext(this, seqDelegate)));
      exposedFeatures.addFeatureHolder(addedFeatures);
    } catch (ChangeVetoException cve) {
      throw new BioError("Modification of hidden featureholder vetoed!", cve);
    }

    urn = seqDelegate.getURN();
    if (urn.indexOf('?') >= 0)
      urn = urn + "&view=" + hashCode();
    else
      urn = urn + "?view=" + hashCode();

    anno = new OverlayAnnotation(seqDelegate.getAnnotation());

    featureRealizer = FeatureImpl.DEFAULT;
  }

    /**
     * Construct a view onto an existing sequence which takes on that
     * sequence's name.
     * <p>The prefered method is SequenceTools.view(Sequence seq)
     */
    public ViewSequence(Sequence seq) {
        this(seq, seq.getName());
    }

    /**
     * Construct a view onto a sequence, using a specific FeatureRealizer.
     *
     * <p>The prefered method is SequenceTools.view(Sequence seq, FeatureRealizer fr)
     */
    public ViewSequence(Sequence seq, FeatureRealizer fr) {
        this(seq);
        this.featureRealizer = fr;
    }

    //
    // We implement SymbolList by delegation
    //

    public Alphabet getAlphabet() {
        return seqDelegate.getAlphabet();
    }

    public Iterator iterator() {
        return seqDelegate.iterator();
    }

    public int length() {
        return seqDelegate.length();
    }

    public String seqString() {
        return seqDelegate.seqString();
    }

    public String subStr(int start, int end) {
        return seqDelegate.subStr(start, end);
    }

    public SymbolList subList(int start, int end) {
        return seqDelegate.subList(start, end);
    }

    public Symbol symbolAt(int indx) {
        return seqDelegate.symbolAt(indx);
    }

    public List toList() {
        return seqDelegate.toList();
    }

    //
    // ID methods -- we have our own.
    //

    public String getURN() {
        return urn;
    }

    public String getName() {
        return name;
    }

    //
    // Basic FeatureHolder methods -- delegate to exposedFeatures
    //

    public int countFeatures() {
        return exposedFeatures.countFeatures();
    }

    public Iterator features() {
        return exposedFeatures.features();
    }

    public FeatureHolder filter(FeatureFilter fc, boolean recurse) {
        return exposedFeatures.filter(fc, recurse);
    }

    public FeatureHolder filter(FeatureFilter fc) {
        return exposedFeatures.filter(fc);
    }

    public FeatureFilter getSchema() {
        return exposedFeatures.getSchema();
    }

    //
    // MutableFeatureHolder methods -- delegate to addedFeatures
    //

    /**
     * Remove a feature from this sequence.  <strong>NOTE:</strong> This
     * method will only succeed for features which were added to this
     * ViewSequence.  Trying to remove a Feature from the underlying
     * sequence will cause an IllegalArgumentException.  I think this
     * is the correct behaviour.
     */

    public void removeFeature(Feature f)
        throws ChangeVetoException
    {
      addedFeatures.removeFeature(f);
    }

    public boolean containsFeature(Feature f) {
      return exposedFeatures.containsFeature(f);
    }

    //
    // Get our annotation
    //

    public Annotation getAnnotation() {
        return anno;
    }

    //
    // Feature realization stuff
    //

    public Feature realizeFeature(FeatureHolder parent, Feature.Template template)
        throws BioException
    {
        return featureRealizer.realizeFeature(this, parent, template);
    }

    public Feature createFeature(Feature.Template template)
        throws BioException, ChangeVetoException
    {
      Location loc = template.location;
      if(loc.getMin() < 1 || loc.getMax() > this.length()) {
          throw new BioException("Failed to create a feature with a location "
                                 + loc
                                 + " outside the sequence: name '"
                                 + getName()
                                 + "', URN '"
                                 + getURN()
                                 + "' length "
                                 + length());
      }
      Feature f = realizeFeature(this, template);
      addedFeatures.addFeature(f);
      return f;
    }

    public FeatureHolder getAddedFeatures() {
        return addedFeatures;
    }

  public void edit(Edit edit) throws ChangeVetoException {
    throw new ChangeVetoException("ViewSequence is immutable");
  }
}
