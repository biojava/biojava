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

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.projection.ProjectedFeatureHolder;
import org.biojava.bio.seq.projection.TranslateFlipContext;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.ontology.OntoTools;
import org.biojava.ontology.Term;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;

/**
 * ComponentFeature implementation used by SimpleAssembly.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.1
 */

class SimpleComponentFeature
        extends Unchangeable
        implements ComponentFeature
{
  private FeatureHolder parent;

  private FeatureHolder projectedFeatures;

  private Location location;
  private String type;
  private String source;
  private Term typeTerm;
  private Term sourceTerm;
  private Annotation annotation;

  private StrandedFeature.Strand strand;

  private Sequence componentSequence;
  private Location componentLocation;
  private int translation;

  public SimpleComponentFeature(FeatureHolder parent,
                                ComponentFeature.Template temp)
          throws BioException {
    if (temp.componentSequence == null) {
      throw new BioException("This implementation of ComponentFeature requires that sequences must already be resolved");
    }

    if (locationContent(temp.location) !=
            locationContent(temp.componentLocation)) {
      throw new BioException("Component and container locations must contain an equal number of symbols.");
    }

    if (!temp.location.isContiguous() || !temp.componentLocation.isContiguous()) {
      throw new BioException("Can only include contiguous segments in an assembly [may change in future]");
    }

    this.parent = parent;

    this.location = temp.location;
    this.type = temp.type;
    this.source = temp.source;
    this.typeTerm = (temp.typeTerm == null) ? OntoTools.ANY : temp.typeTerm;
    this.sourceTerm = (temp.sourceTerm == null) ? OntoTools.ANY : temp.sourceTerm;
    this.annotation = temp.annotation;
    this.strand = temp.strand;

    this.componentSequence = temp.componentSequence;
    this.componentLocation = temp.componentLocation;

    if (temp.strand == StrandedFeature.NEGATIVE) {
      this.translation = temp.location.getMax() + temp.componentLocation.getMin();
      this.projectedFeatures = new ProjectedFeatureHolder(
              new TranslateFlipContext(
                      this, componentSequence, translation, true));
    } else if (temp.strand == StrandedFeature.POSITIVE) {
      this.translation = temp.location.getMin() - temp.componentLocation.getMin();
      this.projectedFeatures = new ProjectedFeatureHolder(
              new TranslateFlipContext(
                      this, componentSequence, translation, false));
    } else {
      throw new BioException("Strand must be specified when creating a ComponentFeature");
    }
  }

  public Feature.Template makeTemplate() {
    ComponentFeature.Template cft = new ComponentFeature.Template();
    cft.location = getLocation();
    cft.type = getType();
    cft.source = getSource();
    cft.annotation = getAnnotation();
    cft.strand = getStrand();
    cft.componentSequence = getComponentSequence();
    cft.componentLocation = getComponentLocation();
    return cft;
  }

  private int locationContent(Location l) {
    if (l.isContiguous())
      return l.getMax() - l.getMin() + 1;
    int content = 0;
    for (Iterator i = l.blockIterator(); i.hasNext();) {
      Location sl = (Location) i.next();
      content += (sl.getMax() - sl.getMin() + 1);
    }
    return content;
  }

  /**
   * SimpleComponentFeatures are always constructed with an explicit
   * Sequence object, and are thus resolvable by definition.
   */

  public boolean isComponentResolvable() {
    return true;
  }

  public String getComponentSequenceName() {
    return getComponentSequence().getName();
  }

  public StrandedFeature.Strand getStrand() {
    return strand;
  }

  public void setStrand(Strand strand)
          throws ChangeVetoException {
    throw new ChangeVetoException(
            new ChangeEvent(this, STRAND, strand, this.strand),
            "Can't change strand as it is immutable"
    );
  }

  public Location getLocation() {
    return location;
  }

  public void setLocation(Location loc)
          throws ChangeVetoException {
    throw new ChangeVetoException(
            new ChangeEvent(this, LOCATION, loc, this.location),
            "Can't change location as it is immutable"
    );
  }

  public FeatureHolder getParent() {
    return parent;
  }

  public Sequence getSequence() {
    FeatureHolder fh = parent;
    while (fh instanceof Feature)
      fh = ((Feature) fh).getParent();
    return (Sequence) fh;
  }

  public String getType() {
    if (typeTerm != OntoTools.ANY) {
      return typeTerm.getName();
    } else {
      return type;
    }
  }

  public Term getTypeTerm() {
    return typeTerm;
  }

  public void setSource(String source)
          throws ChangeVetoException {
    throw new ChangeVetoException(
            new ChangeEvent(this, SOURCE, source, this.source),
            "Can't change source as it is immutable"
    );
  }

  public void setSourceTerm(Term source)
          throws ChangeVetoException {
    throw new ChangeVetoException(
            new ChangeEvent(this, SOURCE, source, this.sourceTerm),
            "Can't change source as it is immutable"
    );
  }

  public String getSource() {
    if (sourceTerm != OntoTools.ANY) {
      return sourceTerm.getName();
    } else {
      return source;
    }
  }

  public Term getSourceTerm() {
    return sourceTerm;
  }


  public void setType(String type)
          throws ChangeVetoException {
    throw new ChangeVetoException(
            new ChangeEvent(this, TYPE, type, this.type),
            "Can't change type as it is immutable"
    );
  }

  public void setTypeTerm(Term type)
          throws ChangeVetoException {
    throw new ChangeVetoException(
            new ChangeEvent(this, SOURCE, type, this.typeTerm),
            "Can't change type as it is immutable"
    );
  }

  public Annotation getAnnotation() {
    return annotation;
  }

  public SymbolList getSymbols() {
    SymbolList syms = componentLocation.symbols(componentSequence);
//  	if (strand == StrandedFeature.NEGATIVE) {
//  	    try {
//  		syms = DNATools.reverseComplement(syms);
//  	    } catch (IllegalAlphabetException ex) {
//  		throw new BioError(ex);
//  	    }
//  	}
    return syms;
  }

  public Sequence getComponentSequence() {
    return componentSequence;
  }

  public Location getComponentLocation() {
    return componentLocation;
  }

  protected FeatureHolder getProjectedFeatures() {
//  if (projectedFeatures == null) {
//  	    projectedFeatures = new ProjectedFeatureHolder(componentSequence,
//  							   this, translation);
//  	}
    return projectedFeatures;
  }

  public int countFeatures() {
    return componentSequence.countFeatures();
  }

  public Iterator features() {
    return getProjectedFeatures().features();
  }

  public boolean containsFeature(Feature f) {
    return getProjectedFeatures().containsFeature(f);
  }

  public FeatureHolder filter(FeatureFilter ff) {
    FeatureFilter childFilter = new FeatureFilter.And(new FeatureFilter.ByParent(new FeatureFilter.ByClass(ComponentFeature.class)),
                                                      new FeatureFilter.Not(FeatureFilter.top_level));

    if (FilterUtils.areDisjoint(ff, childFilter)) {
      return FeatureHolder.EMPTY_FEATURE_HOLDER;
    } else {
      return getProjectedFeatures().filter(ff);
    }
  }

  public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
    FeatureFilter childFilter = new FeatureFilter.And(new FeatureFilter.ByParent(new FeatureFilter.ByClass(ComponentFeature.class)),
                                                      new FeatureFilter.Not(FeatureFilter.top_level));

    if (FilterUtils.areDisjoint(ff, childFilter)) {
      return FeatureHolder.EMPTY_FEATURE_HOLDER;
    } else {
      return getProjectedFeatures().filter(ff, recurse);
    }
  }

  public Feature createFeature(Feature.Template temp)
          throws BioException {
    throw new BioException("Can't create features in a ComponentFeature (yet?)");
  }

  public void removeFeature(Feature f) {
    throw new UnsupportedOperationException("Can't remove features from a ComponentFeature.");
  }

  public FeatureFilter getSchema() {
    return new FeatureFilter.ByParent(new FeatureFilter.ByFeature(this));
  }
}
