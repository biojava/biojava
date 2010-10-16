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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.AnnotationType;
import org.biojava.bio.BioError;
import org.biojava.bio.CardinalityConstraint;
import org.biojava.bio.CollectionConstraint;
import org.biojava.bio.PropertyConstraint;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ClassTools;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * StAX handler for parsing AnnotationTypes in FilterXML documents.  Mainly
 * used internally by XMLFilterHandler.
 *
 * @author Thomas Down
 * @since 1.3
 */

public class XMLAnnotationTypeHandler extends StAXContentHandlerBase {
    private Map handlerFactories = new HashMap();
    private Map collectionHandlerFactories = new HashMap();
    private AnnotationType annotationType = new AnnotationType.Impl();
    private int depth = 0;

    /**
     * Return the AnnotationType built by this handler
     */

    public AnnotationType getAnnotationType() {
        return annotationType;
    }

    /**
     * Handler for an individual <code>PropertyConstraint</code> in an AnnotationType.
     * Implement this if you want to add support for additional types of
     * PropertyConstraint.
     */

    public static interface PropertyConstraintHandler extends StAXContentHandler {
        public PropertyConstraint getPropertyConstraint() throws SAXException;
    }

    /**
     * Handler Factory for a <code>PropertyConstraint</code> in an AnnotationType.
     * Implement this if you want to add support for additional types of
     * PropertyConstraint.
     */

    public static interface PropertyConstraintHandlerFactory {
        public PropertyConstraintHandler makeHandler(String nsURI,
                                                     String localName)
             throws SAXException;
    }

    /**
     * Handler for an individual <code>CollectionConstraint</code> in an AnnotationType.
     * Implement this if you want to add support for additional types of
     * CollectionConstraint.
     */

    public static interface CollectionConstraintHandler extends StAXContentHandler {
        public CollectionConstraint getCollectionConstraint() throws SAXException;
    }

    /**
     * Handler Factory for a <code>CollectionConstraint</code> in an AnnotationType.
     * Implement this if you want to add support for additional types of
     * CollectionConstraint.
     */

    public static interface CollectionConstraintHandlerFactory {
        public CollectionConstraintHandler makeHandler(String nsURI,
                                                       String localName)
             throws SAXException;
    }

    /**
     * Register a factory used to create handlers for the specified tag in an
     * XML AnnotationType
     */

    public void registerPropertyHandlerFactory(String nsURI, String localName, PropertyConstraintHandlerFactory factory) {
        handlerFactories.put(new QName(nsURI, localName), factory);
    }

    /**
     * Register a factory used to create handlers for the specified tag in an
     * XML AnnotationType
     */

    public void registerCollectionHandlerFactory(String nsURI, String localName, CollectionConstraintHandlerFactory factory) {
        collectionHandlerFactories.put(new QName(nsURI, localName), factory);
    }



    private PropertyConstraintHandler getHandler(String nsURI, String localName)
        throws SAXException
    {
        PropertyConstraintHandlerFactory factory = (PropertyConstraintHandlerFactory) handlerFactories.get(new QName(nsURI, localName));
        if (factory != null) {
            return factory.makeHandler(nsURI, localName);
        } else {
            throw new SAXException("Unrecognized element (property context) " + nsURI + ":" + localName);
        }
    }

    private CollectionConstraintHandler getCollectionHandler(String nsURI, String localName)
        throws SAXException
    {
        CollectionConstraintHandlerFactory factory = (CollectionConstraintHandlerFactory) collectionHandlerFactories.get(new QName(nsURI, localName));
        if (factory != null) {
            return factory.makeHandler(nsURI, localName);
        } else {
            throw new SAXException("Unrecognized element (collection context) " + nsURI + ":" + localName);
        }
    }

    /**
     * Construct a new XMLAnnotationTypeHandler which can parse the builtin PropertyConstraints.
     */

    public XMLAnnotationTypeHandler() {
        registerPropertyHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "any",
            new CDATAHandlerFactory() {
                protected PropertyConstraint stringToConstraint(String s) {
                    return PropertyConstraint.ANY;
                }
            }
        );
        registerPropertyHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "none",
            new CDATAHandlerFactory() {
                protected PropertyConstraint stringToConstraint(String s) {
                    return PropertyConstraint.NONE;
                }
            }
        );
        registerPropertyHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "value",
            new CDATAHandlerFactory() {
                protected PropertyConstraint stringToConstraint(String s) {
                    return new PropertyConstraint.ExactValue(s);
                }
            }
        );
        registerPropertyHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "byClass",
            new CDATAHandlerFactory() {
                protected PropertyConstraint stringToConstraint(String s)
                    throws SAXException
                {
                    try {
                        return new PropertyConstraint.ByClass(ClassTools.getClassLoader(this).loadClass(s));
                    } catch (Exception ex) {
                        throw new SAXException("Couldn't find class " + s);
                    }
                }
            }
        );
        registerPropertyHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "and",
            new ConstraintsHandlerFactory() {
                public PropertyConstraint constraintsToConstraint(List l)
                    throws SAXException
                {
                    PropertyConstraint pc = PropertyConstraint.ANY;
                    Iterator i = l.iterator();
                    pc = (PropertyConstraint) i.next();
                    while (i.hasNext()) {
                        pc = new PropertyConstraint.And(pc, (PropertyConstraint) i.next());
                    }
                    return pc;
                }
            }
        );
        registerPropertyHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "or",
            new ConstraintsHandlerFactory() {
                public PropertyConstraint constraintsToConstraint(List l)
                    throws SAXException
                {
                    boolean notValue = false;
                    for (Iterator i = l.iterator(); i.hasNext(); ) {
                        if (! (i.next() instanceof PropertyConstraint.ExactValue)) {
                            notValue = true;
                        }
                    }
                    if (notValue) {
                        PropertyConstraint pc = PropertyConstraint.ANY;
                        Iterator i = l.iterator();
                        pc = (PropertyConstraint) i.next();
                        while (i.hasNext()) {
                            pc = new PropertyConstraint.Or(pc, (PropertyConstraint) i.next());
                        }
                        return pc;
                    } else {
                        Set values = new HashSet();
                        for (Iterator i = l.iterator(); i.hasNext(); ) {
                            values.add(((PropertyConstraint.ExactValue) i.next()).getValue());
                        }
                        return new PropertyConstraint.Enumeration(values);
                    }
                }
            }
        );
        registerPropertyHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "byAnnotationType",
            new ByAnnotationTypeHandlerFactory()
        );


        registerCollectionHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "any",
            new CDATACollectionHandlerFactory() {
                protected CollectionConstraint stringToConstraint(String s) {
                    return CollectionConstraint.ANY;
                }
            }
        );
        registerCollectionHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "none",
            new CDATACollectionHandlerFactory() {
                protected CollectionConstraint stringToConstraint(String s) {
                    return CollectionConstraint.NONE;
                }
            }
        );
        registerCollectionHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "empty",
            new CDATACollectionHandlerFactory() {
                protected CollectionConstraint stringToConstraint(String s) {
                    return CollectionConstraint.EMPTY;
                }
            }
        );
        registerCollectionHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "and",
            new ConstraintsCollectionHandlerFactory() {
                public CollectionConstraint constraintsToConstraint(List l)
                    throws SAXException
                {
                    CollectionConstraint pc = CollectionConstraint.ANY;
                    Iterator i = l.iterator();
                    pc = (CollectionConstraint) i.next();
                    while (i.hasNext()) {
                        pc = new CollectionConstraint.And(pc, (CollectionConstraint) i.next());
                    }
                    return pc;
                }
            }
        );
        registerCollectionHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "or",
            new ConstraintsCollectionHandlerFactory() {
                public CollectionConstraint constraintsToConstraint(List l)
                    throws SAXException
                {
                    CollectionConstraint pc = CollectionConstraint.ANY;
                    Iterator i = l.iterator();
                    pc = (CollectionConstraint) i.next();
                    while (i.hasNext()) {
                        pc = new CollectionConstraint.Or(pc, (CollectionConstraint) i.next());
                    }
                    return pc;
                }
            }
        );
        registerCollectionHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "allValuesIn",
            new PropCardCollectionHandlerFactory() {
                public CollectionConstraint constraintsToConstraint(Location card, PropertyConstraint pc)
                    throws SAXException
                {
                    return new CollectionConstraint.AllValuesIn(pc, card);
                }
            }
        );
        registerCollectionHandlerFactory(
            XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS,
            "contains",
            new PropCardCollectionHandlerFactory() {
                public CollectionConstraint constraintsToConstraint(Location card, PropertyConstraint pc)
                    throws SAXException
                {
                    return new CollectionConstraint.Contains(pc, card);
                }
            }
        );
    }

    public void startElement(String nsURI,
                                                 String localName,
                                     String qName,
                                     Attributes attrs,
                                     DelegationManager dm)
                throws SAXException
    {
        if (depth == 1) {
            if (localName.equals("propertyDefault")) {
                // System.err.println("Handling default constraints");
                dm.delegate(new PropertyHandler() {
                    protected void setConstraint(CollectionConstraint cc)
                        throws ChangeVetoException
                    {
                        annotationType.setDefaultConstraint(cc);
                    }
                } );
            } else if (localName.equals("property")) {
                final Object propName = attrs.getValue(XMLAnnotationTypeWriter.XML_ANNOTATIONTYPE_NS, "name");
                // System.err.println("Handling constraints on " + propName);
                dm.delegate(new PropertyHandler() {
                    protected void setConstraint(CollectionConstraint cc)
                        throws ChangeVetoException
                    {
                        annotationType.setConstraint(propName, cc);
                    }
                } );
            } else {
                throw new SAXException("Unexpected element " + nsURI + ":" + localName);
            }
        }
        ++depth;
    }

    public void endElement(String nsURI,
                                               String localName,
                                               String qName,
                                               StAXContentHandler delegate)
                 throws SAXException
    {
        --depth;
    }

    private abstract class PropertyHandler extends StAXContentHandlerBase {
        private CollectionConstraint cons;
        private int depth = 0;

        public void startElement(String nsURI,
                                                 String localName,
                                     String qName,
                                     Attributes attrs,
                                     DelegationManager dm)
                throws SAXException
        {
            if (depth == 1) {
                dm.delegate(getCollectionHandler(nsURI, localName));
            }

            ++depth;
        }

        public void endElement(String nsURI,
                                               String localName,
                                               String qName,
                                               StAXContentHandler delegate)
                 throws SAXException
        {
            --depth;

            /*

            if (delegate instanceof CardinalityHandler) {
                cardinality = ((CardinalityHandler) delegate).getCardinality();
            } else if (delegate instanceof PropertyConstraintHandler) {
                cons = ((PropertyConstraintHandler) delegate).getPropertyConstraint();
            }

            */
            if (delegate instanceof CollectionConstraintHandler) {
                cons = ((CollectionConstraintHandler) delegate).getCollectionConstraint();
            }
        }

        public void endTree() {
            try {
                setConstraint(cons);
            } catch (ChangeVetoException ex) {
                throw new BioError("Assertion failed: couldn't modify AnnotationType", ex);
            }
        }

        protected abstract void setConstraint(CollectionConstraint cc) throws ChangeVetoException;
    }

    private class CardinalityHandler extends StAXContentHandlerBase {
        private List spans = new ArrayList();

        public void startElement(String nsURI,
                                                 String localName,
                                     String qName,
                                     Attributes attrs,
                                     DelegationManager dm)
                throws SAXException
        {
            // System.err.println("CardinalityHandler: " + localName);
            if ("span".equals(localName)) {
                int start = Integer.parseInt(attrs.getValue("start"));
                String stops = attrs.getValue("stop");
                int stop;
                if (stops.equals("infinity")) {
                    stop = Integer.MAX_VALUE;
                } else {
                    stop = Integer.parseInt(stops);
                }
                spans.add(new RangeLocation(start, stop));
            }
        }

        public Location getCardinality() {
            // System.err.println("GetCardinality");
            return LocationTools.union(spans);
        }
    }

    private abstract class CDATAHandlerFactory implements PropertyConstraintHandlerFactory {
        private class CDATAHandler extends StringElementHandlerBase implements PropertyConstraintHandler {
            private PropertyConstraint cons;

            public PropertyConstraint getPropertyConstraint() {
                return cons;
            }

            protected void setStringValue(String s)
                throws SAXException
            {
                cons = stringToConstraint(s);
                // System.err.println("stringToConstraint returned " + cons);
            }
        }

        protected abstract PropertyConstraint stringToConstraint(String s) throws SAXException;

        public PropertyConstraintHandler makeHandler(String nsURI, String localName) {
            // System.err.println("Making CDATAHandler for " + localName);
            return new CDATAHandler();
        }
    }

    private abstract class ConstraintsHandlerFactory implements PropertyConstraintHandlerFactory {
        private class ConstraintsHandler extends StAXContentHandlerBase implements PropertyConstraintHandler {
            private List constraintChildren = new ArrayList();
            private int depth = 0;

            public void startElement(String nsURI,
                                                 String localName,
                                     String qName,
                                     Attributes attrs,
                                     DelegationManager dm)
                throws SAXException
            {
                if (depth == 1) {
                    PropertyConstraintHandler childHandler = getHandler(nsURI, localName);
                    dm.delegate(childHandler);
                }
                ++depth;
            }

            public void endElement(String nsURI,
                                               String localName,
                                               String qName,
                                               StAXContentHandler delegate)
                 throws SAXException
            {
                if (delegate instanceof PropertyConstraintHandler) {
                    constraintChildren.add(((PropertyConstraintHandler) delegate).getPropertyConstraint());
                }
                --depth;
            }

            public PropertyConstraint getPropertyConstraint()
                throws SAXException
            {
                return constraintsToConstraint(constraintChildren);
            }
        }

        protected abstract PropertyConstraint constraintsToConstraint(List filters) throws SAXException;

        public PropertyConstraintHandler makeHandler(String nsURI, String localName) {
            return new ConstraintsHandler();
        }
    }

    private class ByAnnotationTypeHandlerFactory implements PropertyConstraintHandlerFactory {
        private class ByAnnotationTypeHandler extends StAXContentHandlerBase implements PropertyConstraintHandler {
            private int depth = 0;
            private AnnotationType annoType;

            public void startElement(String nsURI,
                                                 String localName,
                                     String qName,
                                     Attributes attrs,
                                     DelegationManager dm)
                throws SAXException
            {
                if (depth == 1) {
                    dm.delegate(new XMLAnnotationTypeHandler());
                }
                ++depth;
            }

            public void endElement(String nsURI,
                                               String localName,
                                               String qName,
                                               StAXContentHandler delegate)
                 throws SAXException
            {
                if (delegate instanceof XMLAnnotationTypeHandler) {
                    annoType = ((XMLAnnotationTypeHandler) delegate).getAnnotationType();
                }
                --depth;
            }

            public PropertyConstraint getPropertyConstraint()
                throws SAXException
            {
                return new PropertyConstraint.ByAnnotationType(annoType);
            }
        }

        public PropertyConstraintHandler makeHandler(String nsURI, String localName) {
            return new ByAnnotationTypeHandler();
        }
    }

    private abstract class CDATACollectionHandlerFactory implements CollectionConstraintHandlerFactory {
        private class CDATAHandler extends StringElementHandlerBase implements CollectionConstraintHandler {
            private CollectionConstraint cons;

            public CollectionConstraint getCollectionConstraint() {
                return cons;
            }

            protected void setStringValue(String s)
                throws SAXException
            {
                cons = stringToConstraint(s);
                // System.err.println("stringToConstraint returned " + cons);
            }
        }

        protected abstract CollectionConstraint stringToConstraint(String s) throws SAXException;

        public CollectionConstraintHandler makeHandler(String nsURI, String localName) {
            // System.err.println("Making CDATAHandler for " + localName);
            return new CDATAHandler();
        }
    }

    private abstract class ConstraintsCollectionHandlerFactory implements CollectionConstraintHandlerFactory {
        private class ConstraintsHandler extends StAXContentHandlerBase implements CollectionConstraintHandler {
            private List constraintChildren = new ArrayList();
            private int depth = 0;

            public void startElement(String nsURI,
                                                 String localName,
                                     String qName,
                                     Attributes attrs,
                                     DelegationManager dm)
                throws SAXException
            {
                if (depth == 1) {
                    CollectionConstraintHandler childHandler = getCollectionHandler(nsURI, localName);
                    dm.delegate(childHandler);
                }
                ++depth;
            }

            public void endElement(String nsURI,
                                               String localName,
                                               String qName,
                                               StAXContentHandler delegate)
                 throws SAXException
            {
                if (delegate instanceof CollectionConstraintHandler) {
                    constraintChildren.add(((CollectionConstraintHandler) delegate).getCollectionConstraint());
                }
                --depth;
            }

            public CollectionConstraint getCollectionConstraint()
                throws SAXException
            {
                return constraintsToConstraint(constraintChildren);
            }
        }

        protected abstract CollectionConstraint constraintsToConstraint(List filters) throws SAXException;

        public CollectionConstraintHandler makeHandler(String nsURI, String localName) {
            return new ConstraintsHandler();
        }
    }

    private abstract class PropCardCollectionHandlerFactory implements CollectionConstraintHandlerFactory {
        private class PropCardHandler extends StAXContentHandlerBase implements CollectionConstraintHandler {
            private Location cardinality;
            private PropertyConstraint pc;
            private int depth = 0;

            public void startElement(String nsURI,
                                                 String localName,
                                     String qName,
                                     Attributes attrs,
                                     DelegationManager dm)
                throws SAXException
            {
                if (depth == 1) {
                    if  ("cardinalityAny".equals(localName)) {
                        cardinality = CardinalityConstraint.ANY;
                    } else if ("cardinalityZero".equals(localName)) {
                        cardinality = CardinalityConstraint.ZERO;
                    } else if ("cardinalityOne".equals(localName)) {
                        cardinality = CardinalityConstraint.ONE;
                    } else if ("cardinalityNone".equals(localName)) {
                        cardinality = CardinalityConstraint.NONE;
                    } else if ("cardinality".equals(localName)) {
                        dm.delegate(new CardinalityHandler());
                    } else {
                        dm.delegate(getHandler(nsURI, localName));
                    }
                }
                ++depth;
            }

            public void endElement(String nsURI,
                                               String localName,
                                               String qName,
                                               StAXContentHandler delegate)
                 throws SAXException
            {
                if (delegate instanceof CardinalityHandler) {
                    cardinality = ((CardinalityHandler) delegate).getCardinality();
                } else if (delegate instanceof PropertyConstraintHandler) {
                    pc = ((PropertyConstraintHandler) delegate).getPropertyConstraint();
                }
                --depth;
            }

            public CollectionConstraint getCollectionConstraint()
                throws SAXException
            {
                return constraintsToConstraint(cardinality, pc);
            }
        }

        protected abstract CollectionConstraint constraintsToConstraint(Location card, PropertyConstraint pc) throws SAXException;

        public CollectionConstraintHandler makeHandler(String nsURI, String localName) {
            return new PropCardHandler();
        }
    }
}
