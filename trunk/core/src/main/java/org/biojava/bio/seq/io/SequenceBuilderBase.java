/**
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

package org.biojava.bio.seq.io;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;

/**
 * Basic SequenceBuilder implementation which accumulates all
 * notified information.  Subclass this to implement specific
 * Sequence implementations.
 *
 * More functionality is offered by the 
 * {@link org.biojavax.bio.seq.io.SimpleRichSequenceBuilder SimpleRichSequenceBuilder}.
 *
 * @author Thomas Down
 * @author David Huen (modified SimpleSequence to make this)
 * @version 1.2 [newio proposal]
 * 
 */

public abstract class SequenceBuilderBase implements SequenceBuilder {
    public static Object ERROR_FEATURES_PROPERTY
      = SequenceBuilderBase.class + "ERROR_FEATURES_PROPERTY";

    //
    // State
    //

    protected String name;
    protected String uri;

    // annotation on the sequence itself
    protected Annotation annotation;

    // features directly attached to sequence
    private Set rootFeatures;

    private List featureStack;

    protected Sequence seq;

    {
        annotation = new SimpleAnnotation();
        rootFeatures = new HashSet();
        featureStack = new ArrayList();
//	slBuilder = new ChunkedSymbolListBuilder();
    }

    //
    // SeqIOListener
    //

    public void startSequence() {
    }

    public void endSequence() {
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setURI(String uri) {
        this.uri = uri;
    }

    public abstract void addSymbols(Alphabet alpha, Symbol[] syms, int pos, int len)
        throws IllegalAlphabetException;

    /**
     * Add an annotation-bundle entry to the sequence.  If the annotation key
     * isn't currently defined, the value is added directly.  Otherwise:
     *
     * <ul>
     * <li> If the current value implements the Collection interface,
     *      the new value is added to that collection. </li>
     * <li> Otherwise, the current value is replaced by a List object
     *      containing the old value then the new value in that order. </li>
     * </ul>
     */
    public void addSequenceProperty(Object key, Object value) {
      addProperty(annotation, key, value);
    }

    public void startFeature(Feature.Template templ) {
        TemplateWithChildren t2 = new TemplateWithChildren();
        t2.template = templ;
    if(templ.annotation == Annotation.EMPTY_ANNOTATION) {
        templ.annotation = new SmallAnnotation();
    }
        int stackSize = featureStack.size();
        if (stackSize == 0) {
            rootFeatures.add(t2);
        } else {
            TemplateWithChildren parent = (TemplateWithChildren) featureStack.get(stackSize - 1);
            if (parent.children == null)
                parent.children = new HashSet();
            parent.children.add(t2);
        }
        featureStack.add(t2);
    }

    /**
     * Add an annotation-bundle entry to the feature. If the annotation key
     * isn't currently defined, the value is added directly. Otherwise:
     *
     * <ul>
     * <li> If the current value implements the Collection interface,
     *      the new value is added to that collection. </li>
     * <li> Otherwise, the current value is replaced by a List object
     *      containing the old value then the new value in that order. </li>
     * </ul>
     */
    public void addFeatureProperty(Object key, Object value)
    throws ParseException {
      try {
        int stackSize = featureStack.size();

        TemplateWithChildren top =
        (TemplateWithChildren) featureStack.get(stackSize - 1);

        addProperty(top.template.annotation, key, value);
      } catch (IndexOutOfBoundsException ioobe) {
        throw new ParseException(
          ioobe,
          "Attempted to add annotation to a feature when no startFeature " +
          "had been invoked"
        );
      }
    }

    public void endFeature() {
        if (featureStack.size() == 0)
            throw new BioError("Assertion failed: Not within a feature");
        featureStack.remove(featureStack.size() - 1);
    }

    public Sequence makeSequence()
            throws BioException
    {
      //	SymbolList symbols = slBuilder.makeSymbolList();
      //	Sequence seq = new SimpleSequence(symbols, uri, name, annotation);
      try {
        for (Iterator i = rootFeatures.iterator(); i.hasNext(); ) {
          TemplateWithChildren twc = (TemplateWithChildren) i.next();
          try {
            Feature f = seq.createFeature(twc.template);
            if (twc.children != null) {
              makeChildFeatures(f, twc.children);
            }
          } catch (Exception e) {
            // fixme: we should do something more sensible with this error
            e.printStackTrace();
            Set errFeatures;
            Annotation ann = seq.getAnnotation();
            if(ann.containsProperty(ERROR_FEATURES_PROPERTY)) {
              errFeatures = (Set) ann.getProperty(ERROR_FEATURES_PROPERTY);
            } else {
              ann.setProperty(
                ERROR_FEATURES_PROPERTY,
                errFeatures = new HashSet()
              );
            }
            errFeatures.add(twc);
          }
        }
      } catch (Exception ex) {
        throw new BioError("Couldn't create feature",ex);
      }
      return seq;
    }

    private void makeChildFeatures(Feature parent, Set children)
        throws Exception
    {
        for (Iterator i = children.iterator(); i.hasNext(); ) {
            TemplateWithChildren twc = (TemplateWithChildren) i.next();
            Feature f = parent.createFeature(twc.template);
            if (twc.children != null) {
                makeChildFeatures(f, twc.children);
            }
        }
    }

    protected void addProperty(Annotation ann, Object key, Object value) {
        if (value == null)
            return;

        Object oldValue = null;
        Object newValue = value;

        if(ann.containsProperty(key)) {
            oldValue = ann.getProperty(key);
        }

        if (oldValue != null) {
            if (oldValue instanceof Collection) {
                ((Collection) oldValue).add(newValue);
                newValue = oldValue;
            } else {
                List nvList = new ArrayList();
                nvList.add(oldValue);
                nvList.add(newValue);
                newValue = nvList;
            }
        }

        try {
            ann.setProperty(key, newValue);
        } catch (ChangeVetoException ex) {
            throw new BioError("Annotation should be modifiable",ex);
        }
    }

    private static class TemplateWithChildren {
        Feature.Template template;
        Set children;
    }
}
