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

package org.biojava.bio.seq.io.filterxml;

import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.AnnotationType;
import org.biojava.bio.BioError;
import org.biojava.bio.CardinalityConstraint;
import org.biojava.bio.CollectionConstraint;
import org.biojava.bio.PropertyConstraint;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.xml.XMLWriter;

/**
 * Main class for writing AnnotationTypes as XML.  Knows about all the builtin
 * classes of AnnotationType.  It's possible to plug new ones in by calling
 * one of the addXMLPropertyConstraintWriter methods.
 *
 * @author Thomas Down
 * @since 1.3
 */

public class XMLAnnotationTypeWriter {
    /**
     * XML namespace string used to the AnnotationType representation
     */

    public static final String XML_ANNOTATIONTYPE_NS = "http://www.biojava.org/AnnotationType";

    private Map constraintWritersByObject = new HashMap();
    private Map constraintWritersByClass = new HashMap();
    private Map colConstraintWritersByObject = new HashMap();
    private Map colConstraintWritersByClass = new HashMap();
    private boolean strict = false;

    /**
     * Writer for types of CollectionConstraint.  Implement this to add support
     * for a new type of CollectionConstraint.
     *
     * @author Thomas Down
     * @since 1.3
     */

    public interface XMLCollectionConstraintWriter {
        public void writeCollectionConstraint(CollectionConstraint pc,
                                              XMLWriter xw,
                                              XMLAnnotationTypeWriter config)
            throws ClassCastException, IOException, IllegalArgumentException;
    }

    /**
     * Writer for types of PropertyConstraint.  Implement this to add support
     * for a new type of PropertyConstraint.
     *
     * @author Thomas Down
     * @since 1.3
     */

    public interface XMLPropertyConstraintWriter {
        public void writePropertyConstraint(PropertyConstraint pc,
                                            XMLWriter xw,
                                            XMLAnnotationTypeWriter config)
            throws ClassCastException, IOException, IllegalArgumentException;
    }

    /**
     * Construct a new AnnotationTypeWriter which knows about the builtin types of PropertyConstraint
     */

    public XMLAnnotationTypeWriter() {
        try {
            addXMLPropertyConstraintWriter(
                    PropertyConstraint.ANY,
                    new BlankConstraintWriter(XML_ANNOTATIONTYPE_NS, "any")
            );
            addXMLPropertyConstraintWriter(
                    PropertyConstraint.NONE,
                    new BlankConstraintWriter(XML_ANNOTATIONTYPE_NS, "none")
            );
            addXMLPropertyConstraintWriter(
                    PropertyConstraint.ExactValue.class,
                    new ExactValueConstraintWriter()
            );
            addXMLPropertyConstraintWriter(
                    PropertyConstraint.ByClass.class,
                    new ByClassConstraintWriter()
            );
            addXMLPropertyConstraintWriter(
                    PropertyConstraint.And.class,
                    new AndConstraintWriter()
            );
            addXMLPropertyConstraintWriter(
                    PropertyConstraint.Or.class,
                    new OrConstraintWriter()
            );
            addXMLPropertyConstraintWriter(
                    PropertyConstraint.Enumeration.class,
                    new EnumConstraintWriter()
            );
            addXMLPropertyConstraintWriter(
                    PropertyConstraint.ByAnnotationType.class,
                    new AnnotationTypeConstraintWriter()
            );

            addXMLCollectionConstraintWriter(
                    CollectionConstraint.ANY,
                    new BlankCollectionConstraintWriter(XML_ANNOTATIONTYPE_NS, "any")
            );
            addXMLCollectionConstraintWriter(
                    CollectionConstraint.NONE,
                    new BlankCollectionConstraintWriter(XML_ANNOTATIONTYPE_NS, "none")
            );
            addXMLCollectionConstraintWriter(
                    CollectionConstraint.EMPTY,
                    new BlankCollectionConstraintWriter(XML_ANNOTATIONTYPE_NS, "empty")
            );
            addXMLCollectionConstraintWriter(
                    CollectionConstraint.And.class,
                    new AndCollectionConstraintWriter()
            );
            addXMLCollectionConstraintWriter(
                    CollectionConstraint.Or.class,
                    new OrCollectionConstraintWriter()
            );
            addXMLCollectionConstraintWriter(
                    CollectionConstraint.AllValuesIn.class,
                    new AllValuesInCollectionConstraintWriter()
            );
            addXMLCollectionConstraintWriter(
                    CollectionConstraint.Contains.class,
                    new ContainsCollectionConstraintWriter()
            );
        } catch (Exception ex) {
            throw new BioError("Assertion failed: couldn't initialize XMLFilterWriters", ex);
        }
    }

    /**
     * Register a writer for the specified class of property constraint
     */

    public void addXMLPropertyConstraintWriter(Class clazz, XMLPropertyConstraintWriter xfw) {
        constraintWritersByClass.put(clazz, xfw);
    }

    /**
     * Register a writer for a singleton property constraint.
     */

    public void addXMLPropertyConstraintWriter(PropertyConstraint pc, XMLPropertyConstraintWriter xfw) {
        constraintWritersByObject.put(pc, xfw);
    }

    /**
     * Register a writer for the specified class of collection constraint
     */

    public void addXMLCollectionConstraintWriter(Class clazz, XMLCollectionConstraintWriter xfw) {
        colConstraintWritersByClass.put(clazz, xfw);
    }

    /**
     * Register a writer for a singleton property constraint.
     */

    public void addXMLCollectionConstraintWriter(CollectionConstraint pc, XMLCollectionConstraintWriter xfw) {
        colConstraintWritersByObject.put(pc, xfw);
    }

    /**
     * Determine if this writer is in strict mode.
     */

    public boolean isStrict() {
        return strict;
    }

    /**
     * Selects strict mode.  In strict mode, the writer will throw an <code>IllegalArgumentException</code>
     * if it encounters a type of <code>PropertyConstraint</code> it doesn't recognize.  When not
     * in strict model, unrecognized constraints are silently replaced by <code>PropertyConstraint.ANY</code>.
     * Default is <code>false</code>.
     */

    public void setIsStrict(boolean b) {
        this.strict = b;
    }

    private void writeCardinality(Location card,
                                  XMLWriter xw)
        throws IllegalArgumentException, IOException
    {
        if (card == CardinalityConstraint.ANY) {
            xw.openTag(XML_ANNOTATIONTYPE_NS, "cardinalityAny");
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "cardinalityAny");
        } else if (card == CardinalityConstraint.ZERO) {
            xw.openTag(XML_ANNOTATIONTYPE_NS, "cardinalityZero");
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "cardinalityZero");
        } else if (card == CardinalityConstraint.ONE) {
            xw.openTag(XML_ANNOTATIONTYPE_NS, "cardinalityOne");
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "cardinalityOne");
        } else if (card == CardinalityConstraint.NONE) {
            xw.openTag(XML_ANNOTATIONTYPE_NS, "cardinalityNone");
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "cardinalityNone");
        } else {
            xw.openTag(XML_ANNOTATIONTYPE_NS, "cardinality");
            for (Iterator bi = card.blockIterator(); bi.hasNext(); ) {
                Location bloc = (Location) bi.next();
                xw.openTag("span");
                xw.attribute("start", intToString(bloc.getMin()));
                xw.attribute("stop",  intToString(bloc.getMax()));
                xw.closeTag("span");
            }
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "cardinality");
        }
    }

    private String intToString(int i) {
        if (i == Integer.MAX_VALUE) {
            return "infinity";
        } else {
            return "" + i;
        }
    }

    /**
     * Writes a single property constraint.
     */

    void writePropertyConstraint(PropertyConstraint pc,
                                 XMLWriter xw)
        throws IllegalArgumentException, IOException
    {
        XMLPropertyConstraintWriter xpcw = (XMLPropertyConstraintWriter) constraintWritersByObject.get(pc);
        if (xpcw == null) {
            xpcw = (XMLPropertyConstraintWriter) constraintWritersByClass.get(pc.getClass());
        }
        if (xpcw == null) {
            if (strict) {
                throw new IllegalArgumentException("Couldn't find a writer for constraint of type " + pc.getClass().getName());
            } else {
                xpcw = (XMLPropertyConstraintWriter) constraintWritersByObject.get(PropertyConstraint.ANY);
            }
        }
        xpcw.writePropertyConstraint(pc, xw, this);
    }

    /**
     * Writes a single collection constraint
     */

    void writeCollectionConstraint(CollectionConstraint pc,
                                   XMLWriter xw)
        throws IllegalArgumentException, IOException
    {
        XMLCollectionConstraintWriter xpcw = (XMLCollectionConstraintWriter) colConstraintWritersByObject.get(pc);
        if (xpcw == null) {
            xpcw = (XMLCollectionConstraintWriter) colConstraintWritersByClass.get(pc.getClass());
        }
        if (xpcw == null) {
            if (strict) {
                throw new IllegalArgumentException("Couldn't find a writer for constraint of type " + pc.getClass().getName());
            } else {
                xpcw = (XMLCollectionConstraintWriter) constraintWritersByObject.get(CollectionConstraint.ANY);
            }
        }
        xpcw.writeCollectionConstraint(pc, xw, this);
    }

    /**
     * Write an <code>AnnotationType</code> to the specified XMLWriter.
     *
     * @throws IllegalArgumentException if the AnnotationType contains unrecognized
     *                                  constraints, and the writer is in strict mode.
     * @throws IOException if an error occurs while outputting XML.
     */

    public void writeAnnotationType(AnnotationType at, XMLWriter xw)
        throws IllegalArgumentException, IOException
    {
        xw.openTag(XML_ANNOTATIONTYPE_NS, "annotationType");
        Set propKeys = at.getProperties();

        xw.openTag(XML_ANNOTATIONTYPE_NS, "propertyDefault");
        CollectionConstraint.AllValuesIn defaultcc = (CollectionConstraint.AllValuesIn) at.getDefaultConstraint();
        // writeCardinality(defaultcc.getCardinalityConstraint(), xw);
        // writePropertyConstraint(defaultcc.getPropertyConstraint(), xw);
        writeCollectionConstraint(defaultcc, xw);
        xw.closeTag(XML_ANNOTATIONTYPE_NS, "propertyDefault");

        for (Iterator pi = propKeys.iterator(); pi.hasNext(); ) {
            Object propKey = pi.next();
            CollectionConstraint cc = at.getConstraint(propKey);

            xw.openTag(XML_ANNOTATIONTYPE_NS, "property");
            xw.attribute(XML_ANNOTATIONTYPE_NS, "name", propKey.toString());
            // CollectionConstraint.AllValuesIn cc = (CollectionConstraint.AllValuesIn) at.getConstraint(propKey);
            // writeCardinality(cc.getCardinalityConstraint(), xw);
            // writePropertyConstraint(cc.getPropertyConstraint(), xw);
            writeCollectionConstraint(cc, xw);
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "property");
        }
        xw.closeTag(XML_ANNOTATIONTYPE_NS, "annotationType");
    }

    private static class BlankConstraintWriter implements XMLPropertyConstraintWriter {
        private String nsURI;
        private String localName;

        BlankConstraintWriter(String nsURI, String localName) {
            this.nsURI = nsURI;
            this.localName =  localName;
        }

        public void writePropertyConstraint(PropertyConstraint pc,
                                            XMLWriter xw,
                                            XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            xw.openTag(nsURI, localName);
            xw.closeTag(nsURI, localName);
        }
    }

    private static class ExactValueConstraintWriter implements XMLPropertyConstraintWriter {
        public void writePropertyConstraint(PropertyConstraint pc,
                                            XMLWriter xw,
                                            XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            PropertyConstraint.ExactValue pcev = (PropertyConstraint.ExactValue) pc;
            xw.openTag(XML_ANNOTATIONTYPE_NS, "value");
            xw.print(pcev.getValue().toString());
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "value");
        }
    }

    private static class ByClassConstraintWriter implements XMLPropertyConstraintWriter {
        public void writePropertyConstraint(PropertyConstraint pc,
                                            XMLWriter xw,
                                            XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            PropertyConstraint.ByClass pcev = (PropertyConstraint.ByClass) pc;
            xw.openTag(XML_ANNOTATIONTYPE_NS, "byClass");
            xw.print(pcev.getPropertyClass().getName());
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "byClass");
        }
    }

    private static class AndConstraintWriter implements XMLPropertyConstraintWriter {
        public void writePropertyConstraint(PropertyConstraint pc,
                                            XMLWriter xw,
                                            XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            PropertyConstraint.And pca = (PropertyConstraint.And) pc;
            xw.openTag(XML_ANNOTATIONTYPE_NS, "and");
            writeSubConstraint(pca, xw, config);
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "and");
        }

        private void writeSubConstraint(PropertyConstraint pc,
                                        XMLWriter xw,
                                        XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            if (pc instanceof PropertyConstraint.And) {
                PropertyConstraint.And ffa = (PropertyConstraint.And) pc;
                writeSubConstraint(ffa.getChild1(), xw, config);
                writeSubConstraint(ffa.getChild2(), xw, config);
            } else {
                config.writePropertyConstraint(pc, xw);
            }
        }
    }

    private static class OrConstraintWriter implements XMLPropertyConstraintWriter {
        public void writePropertyConstraint(PropertyConstraint pc,
                                            XMLWriter xw,
                                            XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            PropertyConstraint.Or pca = (PropertyConstraint.Or) pc;
            xw.openTag(XML_ANNOTATIONTYPE_NS, "or");
            writeSubConstraint(pca, xw, config);
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "or");
        }

        private void writeSubConstraint(PropertyConstraint pc,
                                        XMLWriter xw,
                                        XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            if (pc instanceof PropertyConstraint.Or) {
                PropertyConstraint.Or ffa = (PropertyConstraint.Or) pc;
                writeSubConstraint(ffa.getChild1(), xw, config);
                writeSubConstraint(ffa.getChild2(), xw, config);
            } else {
                config.writePropertyConstraint(pc, xw);
            }
        }
    }

    private static class EnumConstraintWriter implements XMLPropertyConstraintWriter {
        public void writePropertyConstraint(PropertyConstraint pc,
                                            XMLWriter xw,
                                            XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            PropertyConstraint.Enumeration pcenum = (PropertyConstraint.Enumeration) pc;
            xw.openTag(XML_ANNOTATIONTYPE_NS, "or");
            for (Iterator vi = pcenum.getValues().iterator(); vi.hasNext(); ) {
                xw.openTag(XML_ANNOTATIONTYPE_NS, "value");
                xw.print(vi.next().toString());
                xw.closeTag(XML_ANNOTATIONTYPE_NS, "value");
            }
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "or");
        }
    }

    private static class AnnotationTypeConstraintWriter implements XMLPropertyConstraintWriter {
        public void writePropertyConstraint(PropertyConstraint pc,
                                            XMLWriter xw,
                                            XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            PropertyConstraint.ByAnnotationType pcbat = (PropertyConstraint.ByAnnotationType) pc;
            config.writeAnnotationType(pcbat.getAnnotationType(), xw);
        }
    }



    private static class BlankCollectionConstraintWriter implements XMLCollectionConstraintWriter {
        private String nsURI;
        private String localName;

        BlankCollectionConstraintWriter(String nsURI, String localName) {
            this.nsURI = nsURI;
            this.localName =  localName;
        }

        public void writeCollectionConstraint(CollectionConstraint pc,
                                              XMLWriter xw,
                                              XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            xw.openTag(nsURI, localName);
            xw.closeTag(nsURI, localName);
        }
    }

    private static class AndCollectionConstraintWriter implements XMLCollectionConstraintWriter {
        public void writeCollectionConstraint(CollectionConstraint pc,
                                              XMLWriter xw,
                                              XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            CollectionConstraint.And pca = (CollectionConstraint.And) pc;
            xw.openTag(XML_ANNOTATIONTYPE_NS, "and");
            writeSubConstraint(pca, xw, config);
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "and");
        }

        private void writeSubConstraint(CollectionConstraint pc,
                                        XMLWriter xw,
                                        XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            if (pc instanceof CollectionConstraint.And) {
                CollectionConstraint.And ffa = (CollectionConstraint.And) pc;
                writeSubConstraint(ffa.getChild1(), xw, config);
                writeSubConstraint(ffa.getChild2(), xw, config);
            } else {
                config.writeCollectionConstraint(pc, xw);
            }
        }
    }

    private static class OrCollectionConstraintWriter implements XMLCollectionConstraintWriter {
        public void writeCollectionConstraint(CollectionConstraint pc,
                                              XMLWriter xw,
                                              XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            CollectionConstraint.Or pca = (CollectionConstraint.Or) pc;
            xw.openTag(XML_ANNOTATIONTYPE_NS, "or");
            writeSubConstraint(pca, xw, config);
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "or");
        }

        private void writeSubConstraint(CollectionConstraint pc,
                                        XMLWriter xw,
                                        XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            if (pc instanceof CollectionConstraint.Or) {
                CollectionConstraint.Or ffa = (CollectionConstraint.Or) pc;
                writeSubConstraint(ffa.getChild1(), xw, config);
                writeSubConstraint(ffa.getChild2(), xw, config);
            } else {
                config.writeCollectionConstraint(pc, xw);
            }
        }
    }

    private static class AllValuesInCollectionConstraintWriter implements XMLCollectionConstraintWriter {
        public void writeCollectionConstraint(CollectionConstraint pc,
                                              XMLWriter xw,
                                              XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            CollectionConstraint.AllValuesIn cc = (CollectionConstraint.AllValuesIn) pc;
            xw.openTag(XML_ANNOTATIONTYPE_NS, "allValuesIn");
            config.writeCardinality(cc.getCardinalityConstraint(), xw);
            config.writePropertyConstraint(cc.getPropertyConstraint(), xw);
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "allValuesIn");
        }

    }

    private static class ContainsCollectionConstraintWriter implements XMLCollectionConstraintWriter {
        public void writeCollectionConstraint(CollectionConstraint pc,
                                              XMLWriter xw,
                                              XMLAnnotationTypeWriter config)
            throws ClassCastException, IllegalArgumentException, IOException
        {
            CollectionConstraint.Contains cc = (CollectionConstraint.Contains) pc;
            xw.openTag(XML_ANNOTATIONTYPE_NS, "contains");
            config.writeCardinality(cc.getCardinalityConstraint(), xw);
            config.writePropertyConstraint(cc.getPropertyConstraint(), xw);
            xw.closeTag(XML_ANNOTATIONTYPE_NS, "contains");
        }

    }
}
