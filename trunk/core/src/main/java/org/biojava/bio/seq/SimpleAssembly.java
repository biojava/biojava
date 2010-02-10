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

import java.util.Iterator;
import java.util.List;

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.impl.AssembledSymbolList;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * A Sequence which is assembled from other sequences contained
 * in a set of ComponentFeature objects.
 *
 * <p>
 * There is still some potential for optimising SymbolList
 * operations on this class.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.1
 */

public class SimpleAssembly
  extends
    AbstractChangeable
  implements
    Sequence,
    RealizingFeatureHolder
{
    private final String name;
    private final String uri;
    private final Annotation annotation = new SimpleAnnotation();
    private final SimpleFeatureHolder features;
    private final AssembledSymbolList assembly;

    private final FeatureRealizer featureRealizer = org.biojava.bio.seq.impl.FeatureImpl.DEFAULT;

    protected transient ChangeForwarder annotationForwarder;

    {
	assembly = new AssembledSymbolList();
	features = new SimpleFeatureHolder();
    }

    /**
     * Construct a new SimpleAssembly using the DNA alphabet.
     * Initially, the sequence will just contain a list of `N's.
     * Sequence data can be added by adding one or more
     * ComponentFeatures.
     *
     * @param length The length of the sequence
     * @param name The name of the sequence (returned by getName())
     * @param uri The identifier of the sequence (returned by getURN());
     */

    public SimpleAssembly(int length, String name, String uri) {
	this.name = name;
	this.uri = uri;
	assembly.setLength(length);
    }

    /**
     * Construct a new SimpleAssembly using the DNA alphabet.
     * Initially, the sequence will just contain a list of `N's.
     * Sequence data can be added by adding one or more
     * ComponentFeatures.
     *
     * @param name The name of the sequence (returned by getName())
     * @param uri The identifier of the sequence (returned by getURN());
     */

    public SimpleAssembly(String name, String uri) {
	this.name = name;
	this.uri = uri;
    }

    //
    // SymbolList
    //

    public Alphabet getAlphabet() {
	return assembly.getAlphabet();
    }

    public int length() {
	return assembly.length();
    }

    public Symbol symbolAt(int pos) {
	return assembly.symbolAt(pos);
    }

    public SymbolList subList(int start, int end) {
	return assembly.subList(start, end);
    }

    public String seqString() {
	return assembly.seqString();
    }

    public String subStr(int start, int end) {
	return assembly.subStr(start, end);
    }

    public Iterator iterator() {
	return assembly.iterator();
    }

    public List toList() {
	return assembly.toList();
    }

    public void edit(Edit e)
        throws IllegalAlphabetException, ChangeVetoException
    {
	assembly.edit(e);
    }

    //
    // Sequence identification
    //

    public String getName() {
	return name;
    }

    public String getURN() {
	return uri;
    }

    //
    // Annotatable
    //

    public Annotation getAnnotation() {
	return annotation;
    }

    //
    // FeatureHolder
    //

    public Iterator features() {
	return features.features();
    }

    public int countFeatures() {
	return features.countFeatures();
    }

    public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
	    return features.filter(ff, recurse);
    }

    public FeatureHolder filter(FeatureFilter ff) {
	    return features.filter(ff);
    }

    public boolean containsFeature(Feature f) {
      return features.containsFeature(f);
    }

    public Feature createFeature(Feature.Template temp)
        throws BioException, ChangeVetoException
    {
	if (temp.location.getMin() < 1)
	    throw new BioException("Coordinates out of range");

	if (temp instanceof ComponentFeature.Template) {
	    for (Iterator i = assembly.getComponentLocationSet().iterator(); i.hasNext(); ) {
		Location l = (Location) i.next();
		if (l.overlaps(temp.location))
		    throw new BioError("Can't create overlapping ComponentFeature");
	    }
	}

	Feature f = realizeFeature(this, temp);
	features.addFeature(f);
	if (f instanceof ComponentFeature) {
	    ComponentFeature cf = (ComponentFeature) f;
	    Location loc = cf.getLocation();
	    if (loc.getMax() > assembly.length()) {
		assembly.setLength(loc.getMax());
	    }
	    assembly.putComponent(cf);
	}
	return f;
    }

    public void removeFeature(Feature f)
    throws ChangeVetoException {
      if (f instanceof ComponentFeature) {
        assembly.removeComponent(f.getLocation());
      }
      features.removeFeature(f);
    }

    //
    // Feature realization
    //

    public Feature realizeFeature(FeatureHolder fh, Feature.Template temp)
            throws BioException {
      if (temp instanceof ComponentFeature.Template) {
        if (fh != this) {
          throw new BioException("ComponentFeatures can only be attached directly to SimpleAssembly objects");
        }
        ComponentFeature.Template cft = (ComponentFeature.Template) temp;
        return new SimpleComponentFeature(this, cft);
      } else {
        FeatureHolder gopher = fh;
        while (gopher instanceof Feature) {
          if (gopher instanceof ComponentFeature) {
            // fixme: this should delegate onto the ComponentFeature, which
            // should in turn delegate to its ProjectedFeatureHolder and let
            // the projection magic sort this out
            throw new BioException("Cannot [currently] realize features on components of SimpleAssemblies");
          }
          gopher = ((Feature) gopher).getParent();
        }
        return featureRealizer.realizeFeature(this, fh, temp);
      }
    }

    protected ChangeSupport getChangeSupport(ChangeType ct){
      ChangeSupport cs = super.getChangeSupport(ct);

      if(annotationForwarder == null &&
        (ct == null || ct == Annotatable.ANNOTATION)){
        annotationForwarder =
                new ChangeForwarder.Retyper(this, cs, Annotation.PROPERTY);
        getAnnotation().addChangeListener(
            annotationForwarder,
            Annotatable.ANNOTATION);
      }
      return cs;
    }

    public FeatureFilter getSchema() {
        return FeatureFilter.top_level;
    }
}
