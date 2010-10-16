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
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FramedFeature;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ClassTools;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * Factory producing content handlers for parsing FilterXML elements.
 *
 * <p>
 * An XMLFilterHandler object is a collection of individual StAX handlers for parsing FilterXML
 * documents.  It uses XMLAnnotationTypeHandler to parse <code>byAnnotationType</code> elements.
 * To handle an individual XML filter, you should call the getStAXContentHandler method
 * <p>
 *
 * <p><strong>Example:</strong></p>
 * <pre>
 *       // Setup
 *       XMLFilterHandler filterHandler = new XMLFilterHandler();
 *       Reader xmlFile = new FileReader("featurefilter.xml");
 *
 *       // Create an XML parser
 *       SAXParserFactory spf = SAXParserFactory.newInstance();
 *       spf.setNamespaceAware(true);
 *       XMLReader parser = spf.newSAXParser().getXMLReader();
 *
 *       // Create a new handler for this document
 *       XMLFilterHandler.FilterHandler handler = filterHandler.getStAXContentHandler();
 *       parser.setContentHandler(new SAX2StAXAdaptor(handler));
 *
 *       // Parse the file and retrieve the FeatureFilter
 *       parser.parse(new InputSource(xmlFile));
 *       FeatureFilter filter = handler.getFeatureFilter();
 * </pre>
 *
 *
 * @author Thomas Down
 * @since 1.3
 */
 
public class XMLFilterHandler {
    Map handlerFactories = new HashMap();
    
    /**
     * StAXContentHandler for a particular type of FeatureFilter.  Implement
     * this interface to add parsing support for new types of FeatureFilter.
     *
     * @author Thomas Down
     */
    
    public static interface FilterHandler extends StAXContentHandler {
        public FeatureFilter getFeatureFilter() throws SAXException;
    }            
    
    /**
     * Factory of StAXContentHandlers for a particular type of FeatureFilter.  Implement
     * this interface to add parsing support for new types of FeatureFilter.
     *
     * @author Thomas Down
     */
    
    public static interface FilterHandlerFactory {
        public FilterHandler makeHandler(String nsURI, 
                                   String localName,
                                   XMLFilterHandler config)
             throws SAXException;
    }
    
    /**
     * Register a factory used to create handlers for the specified tag name
     */
    
    public void registerHandlerFactory(String nsURI, String localName, FilterHandlerFactory factory) {
        handlerFactories.put(new QName(nsURI, localName), factory);
    }
    
    /**
     * Retrieve a <code>FilterHandler</code> for the specified tag name.
     */
    
    public FilterHandler getHandler(String nsURI, String localName) 
        throws SAXException
    {
        FilterHandlerFactory factory = (FilterHandlerFactory) handlerFactories.get(new QName(nsURI, localName));
        if (factory != null) {
            return factory.makeHandler(nsURI, localName, this);
        } else {
            throw new SAXException("No handler for " + nsURI + ":" + localName);
        }
    }
    
    /**
     * Construct a new XMLFilterHandler which can parse the builtin types of <code>FeatureFilter</code>.
     */
    
    public XMLFilterHandler() {
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "byType",
            new CDATAHandlerFactory() {
                protected FeatureFilter stringToFilter(String s) {
                    return new FeatureFilter.ByType(s);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "bySource",
            new CDATAHandlerFactory() {
                protected FeatureFilter stringToFilter(String s) {
                    return new FeatureFilter.BySource(s);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "byClass",
            new CDATAHandlerFactory() {
                protected FeatureFilter stringToFilter(String s) 
                    throws SAXException
                {
                    try {
                        return new FeatureFilter.ByClass(ClassTools.getClassLoader(this).loadClass(s));
                    } catch (Exception ex) {
                        throw new SAXException("Couldn't load class " + s);
                    }
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "byStrand",
            new CDATAHandlerFactory() {
                protected FeatureFilter stringToFilter(String s) 
                    throws SAXException
                {
                    s = s.trim();
                    StrandedFeature.Strand strand = StrandedFeature.UNKNOWN;
                    if ("POSITIVE".equalsIgnoreCase(s)) {
                        strand = StrandedFeature.POSITIVE;
                    } else if ("NEGATIVE".equalsIgnoreCase(s)) {
                        strand = StrandedFeature.NEGATIVE;
                    } 
                    return new FeatureFilter.StrandFilter(strand);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "byFrame",
            new CDATAHandlerFactory() {
                protected FeatureFilter stringToFilter(String s) 
                    throws SAXException
                {
                    int frameNum = Integer.parseInt(s.trim());
                    switch (frameNum) {
                        case 0:
                            return new FeatureFilter.FrameFilter(FramedFeature.FRAME_0);
                        case 1:
                            return new FeatureFilter.FrameFilter(FramedFeature.FRAME_1);
                        case 2:
                            return new FeatureFilter.FrameFilter(FramedFeature.FRAME_2);
                    }
                    throw new SAXException("Invalid frame: " + frameNum);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "all",
            new CDATAHandlerFactory() {
                protected FeatureFilter stringToFilter(String s) {
                    return FeatureFilter.all;
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "none",
            new CDATAHandlerFactory() {
                protected FeatureFilter stringToFilter(String s) {
                    return FeatureFilter.none;
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "isTopLevel",
            new CDATAHandlerFactory() {
                protected FeatureFilter stringToFilter(String s) {
                    return FeatureFilter.top_level;
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "isLeaf",
            new CDATAHandlerFactory() {
                protected FeatureFilter stringToFilter(String s) {
                    return FeatureFilter.leaf;
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "bySequenceName",
            new CDATAHandlerFactory() {
                protected FeatureFilter stringToFilter(String s) {
                    return new FeatureFilter.BySequenceName(s);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "byComponentName",
            new CDATAHandlerFactory() {
                protected FeatureFilter stringToFilter(String s) {
                    return new FeatureFilter.ByComponentName(s);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "not",
            new SingleFilterHandlerFactory() {
                protected FeatureFilter filterToFilter(FeatureFilter ff) {
                    return new FeatureFilter.Not(ff);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "byParent",
            new SingleFilterHandlerFactory() {
                protected FeatureFilter filterToFilter(FeatureFilter ff) {
                    return new FeatureFilter.ByParent(ff);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "byAncestor",
            new SingleFilterHandlerFactory() {
                protected FeatureFilter filterToFilter(FeatureFilter ff) {
                    return new FeatureFilter.ByAncestor(ff);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "byChild",
            new SingleFilterHandlerFactory() {
                protected FeatureFilter filterToFilter(FeatureFilter ff) {
                    return new FeatureFilter.ByChild(ff);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "byDescendant",
            new SingleFilterHandlerFactory() {
                protected FeatureFilter filterToFilter(FeatureFilter ff) {
                    return new FeatureFilter.ByDescendant(ff);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "onlyChildren",
            new SingleFilterHandlerFactory() {
                protected FeatureFilter filterToFilter(FeatureFilter ff) {
                    return new FeatureFilter.OnlyChildren(ff);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "onlyDescendants",
            new SingleFilterHandlerFactory() {
                protected FeatureFilter filterToFilter(FeatureFilter ff) {
                    return new FeatureFilter.OnlyDescendants(ff);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "and",
            new FiltersHandlerFactory() {
                protected FeatureFilter filtersToFilter(List l) {
                    FeatureFilter ff;
                    Iterator i = l.iterator();
                    ff = (FeatureFilter) i.next();
                    while (i.hasNext()) {
                        ff = new FeatureFilter.And(ff, (FeatureFilter) i.next());
                    }
                    return ff;
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "or",
            new FiltersHandlerFactory() {
                protected FeatureFilter filtersToFilter(List l) {
                    FeatureFilter ff;
                    Iterator i = l.iterator();
                    ff = (FeatureFilter) i.next();
                    while (i.hasNext()) {
                        ff = new FeatureFilter.Or(ff, (FeatureFilter) i.next());
                    }
                    return ff;
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "overlapsLocation",
            new LocationHandlerFactory() {
                protected FeatureFilter locationToFilter(Location loc) {
                    return new FeatureFilter.OverlapsLocation(loc);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "containedByLocation",
            new LocationHandlerFactory() {
                protected FeatureFilter locationToFilter(Location loc) {
                    return new FeatureFilter.ContainedByLocation(loc);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "shadowOverlapsLocation",
            new LocationHandlerFactory() {
                protected FeatureFilter locationToFilter(Location loc) {
                    return new FeatureFilter.ShadowOverlapsLocation(loc);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "shadowContainedByLocation",
            new LocationHandlerFactory() {
                protected FeatureFilter locationToFilter(Location loc) {
                    return new FeatureFilter.ShadowContainedByLocation(loc);
                }
            }
        );
        registerHandlerFactory(
            XMLFilterWriter.XML_FILTER_NS,
            "byAnnotationType",
            new AnnotationTypeHandlerFactory()
        );
    }
    
    private abstract static class SingleFilterHandlerFactory extends FiltersHandlerFactory {
        protected abstract FeatureFilter filterToFilter(FeatureFilter ff) throws SAXException;
        protected final FeatureFilter filtersToFilter(List filters)
            throws SAXException
        {
            if (filters.size() != 1) {
                throw new SAXException("Expected 1 child filter but got " + filters.size());
            }
            return filterToFilter((FeatureFilter) filters.get(0));
        }
    }
    
    private abstract static class FiltersHandlerFactory implements FilterHandlerFactory {
        private class FiltersHandler extends StAXContentHandlerBase implements FilterHandler {
            private XMLFilterHandler config;
            private List filterChildren = new ArrayList();
            private int depth = 0;
            
            public FiltersHandler(XMLFilterHandler config) {
                super();
                this.config = config;
            }
            
            public void startElement(String nsURI,
			                         String localName,
                                     String qName,
                                     Attributes attrs,
                                     DelegationManager dm)
                throws SAXException
            {
                if (depth == 1) {
                    FilterHandler childHandler = config.getHandler(nsURI, localName);
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
                if (delegate instanceof FilterHandler) {
                    filterChildren.add(((FilterHandler) delegate).getFeatureFilter());
                }
                --depth;
            }
            
            public FeatureFilter getFeatureFilter() 
                throws SAXException
            {
                return filtersToFilter(filterChildren);
            }
        }
        
        protected abstract FeatureFilter filtersToFilter(List filters) throws SAXException;
        
        public FilterHandler makeHandler(String nsURI, String localName, XMLFilterHandler config) {
            return new FiltersHandler(config);
        }
    }
    
    private abstract static class LocationHandlerFactory implements FilterHandlerFactory {
        private class LocationHandler extends StAXContentHandlerBase implements FilterHandler {
            private List spans = new ArrayList();
            
            public void startElement(String nsURI,
			                         String localName,
                                     String qName,
                                     Attributes attrs,
                                     DelegationManager dm)
                throws SAXException
            {
                if ("span".equals(localName)) {
                    int start = Integer.parseInt(attrs.getValue("start"));
                    int stop = Integer.parseInt(attrs.getValue("stop"));
                    spans.add(new RangeLocation(start, stop));
                }
            }
            
            public FeatureFilter getFeatureFilter() 
                throws SAXException
            {
                return locationToFilter(LocationTools.union(spans));
            }
        }
        
        protected abstract FeatureFilter locationToFilter(Location loc) throws SAXException;
        
        public FilterHandler makeHandler(String nsURI, String localName, XMLFilterHandler config) {
            return new LocationHandler();
        }
    }
    
    private abstract static class CDATAHandlerFactory implements FilterHandlerFactory {
        private class CDATAHandler extends StringElementHandlerBase implements FilterHandler {
            private FeatureFilter ourFilter;
            
            public FeatureFilter getFeatureFilter() {
                return ourFilter;
            }
            
            protected void setStringValue(String s) 
                throws SAXException
            {
                ourFilter = stringToFilter(s);
            }
        }
        
        protected abstract FeatureFilter stringToFilter(String s) throws SAXException;
        
        public FilterHandler makeHandler(String nsURI, String localName, XMLFilterHandler config) {
            return new CDATAHandler();
        }
    }
    
    private static class AnnotationTypeHandlerFactory implements FilterHandlerFactory {
        private class AnnotationTypeHandler extends StAXContentHandlerBase implements FilterHandler {
            private int depth = 0;
            private FeatureFilter filter;
            
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
                    filter = new FeatureFilter.ByAnnotationType(((XMLAnnotationTypeHandler) delegate).getAnnotationType());
                }
                --depth;
            }
            
            public FeatureFilter getFeatureFilter() {
                return filter;
            }
        }
        
        public FilterHandler makeHandler(String nsURI, String localName, XMLFilterHandler config) {
            return new AnnotationTypeHandler();
        }
    }
    
    private class AnyFilterHandler extends StAXContentHandlerBase implements FilterHandler {
        private FeatureFilter filter;
        
        public FeatureFilter getFeatureFilter() {
            return filter;
        }
        
        public void startElement(String nsURI,
			     String localName,
			     String qName,
			     Attributes attrs,
			     DelegationManager dm)
            throws SAXException
        {
            FilterHandler h = getHandler(nsURI, localName);
            dm.delegate(h);
        }

        public void endElement(String nsURI,
			   String localName,
			   String qName,
			   StAXContentHandler delegate)
            throws SAXException
        {
            filter = ((FilterHandler) delegate).getFeatureFilter();
        }
    }
    
    /**
     * Return a StAXContentHandler which can deal with any FilterXML construct known to this class.
     */
    
    public FilterHandler getStAXContentHandler() {
        return new AnyFilterHandler();
    }
}
