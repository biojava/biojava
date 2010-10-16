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
import java.lang.reflect.Field;
import java.util.Comparator;
import java.util.Iterator;

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.ontology.InvalidTermException;
import org.biojava.ontology.Term;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * A feature within a sequence, or nested within another feature.
 *
 * <h2>Common operations</h2>
 *
 * <pre>
 * // loop over all features in a sequence
 * for(Iterator fi = mySeq.features(); fi.hasNext(); ) {
 *   Feature f = (Feature) fi.next();
 *   System.out.println(f.getType() + "\t" + f.getLocation());
 * }
 *
 * // loop over all features that are children of this one, such as exons in a
 * // gene
 * for(Iterator cfi = feature.features(); cfi.hasNext(); ) {
 *   ...
 *
 * // extract all stranded features that are directly on a sequence
 * FeatureHolder strandedFeatures = mySeq.filter(
 *   new FeatureFilter.ByClass(StrandedFeature.class),
 *   false
 *  );
 * for(fi = strandedFeatures.features(); ...
 *
 * // find all features with the type property set to "EXON" no matter how
 * // far down the feature hierachy they are
 * FeatureHolder repeats = mySeq.filter(
 *   new FeatureFilter.ByType("EXON"),
 *   true;
 * );
 * </pre>
 *
 * <h2>Description</h2>
 *
 *<p> Features contain annotation and a location. The type of the
 * feature is something like 'Repeat' or 'BetaStrand'. Where the
 * feature has been read from an EMBL or Genbank source the type will
 * be the same as the feature key in the feature table e.g. 'gene',
 * 'CDS', 'repeat_unit', 'misc_feature'. The source of the feature is
 * something like 'genscan', 'repeatmasker' or 'made-up'. </p>
 *
 * <p>
 * Features are <em>always</em> contained by a parent <code>FeatureHolder</code>,
 * which may either be a <code>Sequence</code> or another <code>Feature</code>. 
 * Feature instances should never be constructed directly by client
 * code, and the BioJava core does not contain any publicly accessible
 * implementations of the <code>Feature</code> interface.  Instead, you
 * should create a suitable <code>Feature.Template</code>, then pass this
 * to the <code>createFeature</code> method of a <code>Sequence</code>
 * or <code>Feature</code>.
 * </p>
 *
 * <p>
 * We may need some standardisation for what the fields mean. In particular, we
 * should be compliant where sensible with GFF.
 * </p>
 * @see org.biojavax.bio.seq.RichFeature
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Keith James
 */

public interface Feature extends FeatureHolder, Annotatable {

    /**
     * This is used as a key in the <code>Annotation</code> where it
     * identifies internal data. This is not printed when the
     * <code>Feature</code> is written to a flatfile. E.g. the
     * original feature's EMBL location string (if it has one) is
     * stored here.
     */
    public static final String PROPERTY_DATA_KEY = "internal_data";

    /**
     * The location of this feature is being altered.
     */
    public static final ChangeType LOCATION = new ChangeType(
      "Location has altered",
      Feature.class,
      "LOCATION"
    );
    
    /**
     * The type of this feature has altered.
     */
    public static final ChangeType TYPE = new ChangeType(
      "Type has altered",
      Feature.class,
      "TYPE"
    );
    
    /**
     * The source of this feature has altered
     */
    public static final ChangeType SOURCE = new ChangeType(
      "Source has altered",
      Feature.class,
      "SOURCE"
    );
    
    /**
     * The ontological type of this feature has altered.
     */
    public static final ChangeType TYPETERM = new ChangeType(
      "TypeTerm has altered",
      Feature.class,
      "TYPETERM"
    );
    
    /**
     * The ontological source of this feature has altered
     */
    public static final ChangeType SOURCETERM = new ChangeType(
      "SourceTerm has altered",
      Feature.class,
      "SOURCETERM"
    );
    
    /**
     * The location of this feature.
     * <p>
     * The location may be complicated, or simply a range.
     * The annotation is assumed to apply to all the region contained
     * within the location.
     *
     * @return a Location anchoring this feature
     */
    Location getLocation();
    
    /**
     * The new location for this feature.
     * <p>
     * The location may be complicated or simply a range. The annotation is
     * assumed to apply to the entire region contained within the location.
     * Any values returned from methods that rely on the old location must
     * not be affected.
     *
     * @param loc the new Location for this feature
     * @throws ChangeVetoException if the location can't be altered
     */
    void setLocation(Location loc)
        throws ChangeVetoException;
  
    /**
     * The type of the feature.
     *
     * @return the type of this sequence
     */
    String getType();
    
    /**
     * Change the type of this feature.
     *
     * @param type  new type String
     * @throws ChangeVetoException if the type can't be altered
     */
    void setType(String type)
        throws ChangeVetoException;
  
    /**
     * An ontology term defining the type of feature.  This is
     * optional, and will default to <code>OntoTools.ANY</code>
     * in implementations which aren't ontology aware.
     *
     * @since 1.4
     */
    
    public Term getTypeTerm();
    
    /**
     * Set the type ontology-term for this feature.  If this succeeds,
     * it will generally also change the source name.
     *
     * @since 1.4
     * @throws ChangeVetoException if changes are not allowed
     * @throws InvalidTermException if the specified term is not an acceptable type
     *                           for this feature.
     */
    
    public void setTypeTerm(Term t) throws ChangeVetoException, InvalidTermException;
    
    /**
     * An ontology term defining the source of this feature.  This is
     * optional, and will default to <code>OntoTools.ANY</code>
     * in implementations which aren't ontology aware.
     *
     * @since 1.4
     */
    
    public Term getSourceTerm();
    
    /**
     * Set the source ontology-term for this feature.  If this succeeds,
     * it will generally also change the source name.
     *
     * @since 1.4
     * @throws ChangeVetoException if changes are not allowed
     * @throws InvalidTermException if the specified term is not an acceptable type
     *                           for this feature.
     */
    
    public void setSourceTerm(Term t) throws ChangeVetoException, InvalidTermException;
  
    /**
     * The source of the feature. This may be a program or process.
     *
     * @return the source, or generator
     */
    String getSource();
    
    /**
     * Change the source of the Feature.
     *
     * @param source the new source String
     * @throws ChangeVetoException if the source can't be altered
     */
    void setSource(String source)
        throws ChangeVetoException;
    
    /**
     * Return a list of symbols that are contained in this feature.
     * <p>
     * The symbols may not be contiguous in the original sequence, but they
     * will be concatenated together in the resulting SymbolList.
     * <p>
     * The order of the Symbols within the resulting symbol list will be 
     * according to the concept of ordering within the location object.
     * <p>
     * If the feature location is modified then this does not modify any
     * SymbolList produced by earlier invocations of this method.
     *
     * @return  a SymbolList containing each symbol of the parent sequence contained
     *          within this feature in the order they appear in the parent
     */
    SymbolList getSymbols();
  
    /**
     * Return the <code>FeatureHolder</code> to which this feature has been
     * attached.  This will be a <code>Sequence</code> object for top level
     * features, and a <code>Feature</code> object for features further
     * down the tree.
     */

    public FeatureHolder getParent();

    /**
     * Return the <code>Sequence</code> object to which this feature
     * is (ultimately) attached. For top level features, this will be
     * equal to the <code>FeatureHolder</code> returned by
     * <code>getParent</code>.
     *
     * @return the ultimate parent Sequence
     */
    public Sequence getSequence();

    /**
     * Iterate over any child features which are held by this
     * feature.  The order of iteration <em>MAY</em> be significant
     * for some types of Feature.
     */

    public Iterator features();

    /**
     * Create a new Template that could be used to generate a feature identical
     * to this one. The fields of the template can be edited without changing
     * the feature.
     *
     * @return a new Template that would make a feature like this one
     */
    public Template makeTemplate();
    
    /**
     * Template class for a plain feature.
     * <p>
     * This just has fields for representing the properties of a basic Feature. Each
     * sub-interface should provide a template class that inherits off this, and
     * the constructor or factory methods should make a particular feature
     * implementation from the template.
     *
     * <p>
     * The equals(), hashcode(), toString() and populate() methods are defined
     * such that two templates are equal if all their fields are equal.  These
     * are implemented by reflection, and automatically pick up any extra fields
     * added in subclasses.
     * </p>
     *
     * @author Thomas Down
     * @author Matthew Pocock
     */

    public static class Template implements Serializable, Cloneable {
	public Location location;
	public String type;
	public String source;
    public Term typeTerm;
    public Term sourceTerm;
	public Annotation annotation;

      public Object clone() throws CloneNotSupportedException {
        return super.clone();
      }

	public int hashCode() {
	    Class templClazz = getClass();
	    Field[] fields = templClazz.getFields();
	    int hc = 0;
	    for (int i = 0; i < fields.length; ++i) {
		try {
		    Object o = fields[i].get(this);
		    if (o != null) {
			hc += o.hashCode();
		    }
		} catch (Exception ex) {
		    throw new BioError("Can't access template fields", ex);
		}
	    }
	    
	    return hc;
	}

	public boolean equals(Object b) {
	    Class aClazz = getClass();
	    Class bClazz = b.getClass();
	    if (! aClazz.equals(bClazz)) {
		return false;
	    }

	    Field[] fields = aClazz.getFields();
	    for (int i = 0; i < fields.length; ++i) {
		try {
		    Object ao = fields[i].get(this);
		    Object bo = fields[i].get(b);
		    if (ao != bo) {
			if (ao == null) {
			    return false;
			} else {
			    if (! (ao.equals(bo))) {
				return false;
			    }
			}
		    }
		} catch (Exception ex) {
		    throw new BioError("Can't access template fields", ex);
		}
	    }
	    
	    return true;
	}
        

        public String toString() {
          StringBuffer sbuf = new StringBuffer();
          
          Class us = getClass();
          sbuf.append(us);
          sbuf.append(":");
          
          Field[] fields = us.getFields();
          for(int i = 0; i < fields.length; i++) {
            try {
              sbuf.append(" ");
              sbuf.append(fields[i].getName());
              sbuf.append("=");
              sbuf.append(fields[i].get(this));
            } catch (Exception e) {
              throw new BioError("Couldn't access template fields", e);
            }
          }
          
          return sbuf.toString();
        }
    }

    /**
     * <code>byLocationOrder</code> contains a <code>Feature</code>
     * comparator which compares by the minimum base position of their
     * <code>Location</code>.
     */
    public static final ByLocationComparator byLocationOrder =
        new ByLocationComparator();

    /**
     * <code>ByLocationComparator</code> compares
     * <code>Feature</code>s by the minimum base position of their
     * <code>Location</code>.
     *
     * @author Keith James
     * @since 1.2
     */
    public static final class ByLocationComparator implements Comparator
    {
        public int compare(Object o1, Object o2)
        {
            Feature f1 = (Feature) o1;
            Feature f2 = (Feature) o2;

            // We don't subtract one coordinate from another as one or
            // both may be set to Integer.MAX_VALUE/Integer.MIN_VALUE
            // and the result could wrap around. Convert to Long if
            // necessary.
            if (f1.getLocation().getMin() > f2.getLocation().getMin())
                return 1;
            else if (f1.getLocation().getMin() < f2.getLocation().getMin())
                return -1;
            else
                return 0;
        }
    }
}
