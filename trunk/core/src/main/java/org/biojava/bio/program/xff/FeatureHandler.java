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

package org.biojava.bio.program.xff;

import org.biojava.bio.Annotation;
import org.biojava.bio.MergeAnnotation;
import org.biojava.bio.OverlayAnnotation;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.AssertionFailure;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * StAX handler for the basic <code>feature</code> type of XFF.
 * This class can also be subclassed to handle other XFF types.
 *
 * <p>
 * In general, to handle a <code>feature</code> subclass, you will
 * with to:
 * </p>
 *
 * <ul>
 * <li>If necessary, override createFeatureTemplate to build the appropriate BioJava
 *     Feature.Template</li>
 * <li>Add your own <code>startElement</code> and <code>endElement</code>
 *     methods which handle extra extra elements in your feature type.  These
 *     should normally pass on all the standard elements to
 *     <code>super.startElement</code> and <code>super.endElement</code>.</li>
 * </ul>
 *
 * <p>
 * Note that, since <code>FeatureHandler</code> does some basic housekeeping,
 * if you `consume' a startElement notification (i.e. don't pass it on to the
 * superclass) you must also consume the matching endElement.  Since FeatureHandler
 * silently ignores all unrecognized elements, it is usually safe to pass on
 * all startElement and endElement notifications -- even those which are specific
 * to your feature type.
 * </p>
 *
 * @author Thomas Down
 * @since 1.2
 */

public class FeatureHandler extends StAXContentHandlerBase {
    public static final XFFPartHandlerFactory FEATURE_HANDLER_FACTORY = new XFFPartHandlerFactory() {
	    public StAXContentHandler getPartHandler(XFFFeatureSetHandler xffenv) {
		return new FeatureHandler(xffenv);
	    }
	} ;

    private XFFFeatureSetHandler xffenv;
    private Feature.Template template = null;
    private boolean startFired = false;
    private boolean endFired = false;
    private int level = 0;

    /**
     * Construct a new Feature handler, passing in an XFF-parsing environment.
     */

    public FeatureHandler(XFFFeatureSetHandler xffenv) {
	this.xffenv = xffenv;
    }

    /**
     * Return the XFF processing environment passed in when this handler was
     * created.
     */

    public XFFFeatureSetHandler getXFFEnvironment() {
	return xffenv;
    }

    /**
     * Get the template for the feature being constructed.  This should
     * usually not be overridden.  Delegates to <code>createFeatureTemplate</code>
     * for template construction.
     */

    protected Feature.Template getFeatureTemplate() {
	if (template == null) {
	    template = createFeatureTemplate();
	}
	return template;
    }

    /**
     * Create a new template of the appropriate type.  Override this method
     * if you wish to use a template type other than Feature.Template.
     */

    protected Feature.Template createFeatureTemplate() {
	return new Feature.Template();
    }

  /**
   * Fire the startFeature event.  You should wrap this method if you want
   * to perform any validation on the Feature.Template before the
   * startFeature is fired.
   *
   * @throws ParseException if the startFeature notification fails, or if
   *                        it has already been made.
   */

  protected void fireStartFeature()
          throws ParseException
  {
    if (startFired) {
      throw new ParseException("startFeature event has already been fired for this feature");
    }

    Feature.Template templ = getFeatureTemplate();
    Annotation ann = getXFFEnvironment().getMergeAnnotation();
    Annotation orig = templ.annotation;
    Annotation res = null;

    if(ann != null && orig != null) {
      try {
        MergeAnnotation ma = new MergeAnnotation();
        ma.addAnnotation(templ.annotation);
        ma.addAnnotation(ann);
        res = ma;
      } catch (ChangeVetoException cve) {
        throw new AssertionFailure(cve);
      }
    } else if(ann != null) {
      res = ann;
    } else if(orig != null) {
      res = orig;
    }

    if(res == null) {
      res = new SmallAnnotation();
    } else {
      res = new OverlayAnnotation(res);
    }

    templ.annotation = res;

    getXFFEnvironment().getFeatureListener().startFeature(templ);
    startFired = true;
  }

    /**
     * Fire the endFeature event.
     */

    protected void fireEndFeature()
        throws ParseException
    {
	if (!startFired) {
	    throw new ParseException("startFeature has not yet been fired for this feature.");
	}

	if (endFired) {
	    throw new ParseException("endFeature event has already been fired for this feature");
	}

	getXFFEnvironment().getFeatureListener().endFeature();
	endFired = true;
    }

    /**
     * Set a property.  If the startFeature notification has not yet been fired,
     * the property is added directly to the annotation bundle in the feature
     * template, otherwise an addFeatureProperty event is generated.
     */

    protected void setFeatureProperty(Object key, Object value)
        throws ChangeVetoException, ParseException
    {
	if (startFired) {
	    getXFFEnvironment().getFeatureListener().addFeatureProperty(key, value);
	} else {
	    Feature.Template ft = getFeatureTemplate();
	    if (ft.annotation == null) {
		ft.annotation = new SmallAnnotation();
	    }
	    ft.annotation.setProperty(key, value);
	}
    }

    /**
     * StAX callback for element starts.  Wrap this method to handle extra
     * elements within your own feature types.
     */

    public void startElement(String nsURI,
			     String localName,
			     String qName,
			     Attributes attrs,
			     DelegationManager dm)
	 throws SAXException
    {
	level++;
	if (level == 1) {
	    // This was our initial startElement.
	    String id = attrs.getValue("id");
	    if (id != null) {
		try {
		    setFeatureProperty(XFFFeatureSetHandler.PROPERTY_XFF_ID, id);
		} catch (Exception ex) {
		    throw new SAXException("Couldn't set property", ex);
		}
	    }
	} else {
	    if (localName.equals("type")) {
		dm.delegate(getTypeHandler());
	    } else if (localName.equals("source")) {
		dm.delegate(getSourceHandler());
	    } else if (localName.equals("location")) {
		dm.delegate(getLocationHandler());
	    } else if (localName.equals("id")) {
		dm.delegate(getOldIDHandler());
	    } else if (localName.equals("details")) {
		if (!startFired) {
		    try {
			fireStartFeature();
		    } catch (ParseException ex) {
			throw new SAXException(ex);
		    }
		}

		dm.delegate(xffenv.getDetailsHandler());

		// Need to handle details properly
	    } else if (localName.equals("featureSet")) {
		if (!startFired) {
		    try {
			fireStartFeature();
		    } catch (ParseException ex) {
			throw new SAXException(ex);
		    }
		}

		dm.delegate(xffenv);
	    } else {
		// throw new SAXException("Couldn't recognize element " + localName + " in namespace " + nsURI);
	    }
	}
    }

    /**
     * StAX callback for element ends.  Wrap this method to handle extra
     * elements within your own feature types.
     */

    public void endElement(String nsURI,
			   String localName,
			   String qName,
			   StAXContentHandler handler)
	throws SAXException
    {
	level--;

	if (level == 0) {
	    // Our tree is done.

	    try {
		if (!startFired) {
		    fireStartFeature();
		}

		if (!endFired) {
		    fireEndFeature();
		}
	    } catch (ParseException ex) {
		throw new SAXException(ex);
	    }
	}
    }

    protected StAXContentHandler getTypeHandler() {
	return new StringElementHandlerBase() {
		protected void setStringValue(String s) {
		    getFeatureTemplate().type = s.trim();
		}
	    } ;
    }

    protected StAXContentHandler getSourceHandler() {
	return new StringElementHandlerBase() {
		protected void setStringValue(String s) {
		    getFeatureTemplate().source = s.trim();
		}
	    } ;
    }

    protected StAXContentHandler getOldIDHandler() {
	return new StringElementHandlerBase() {
		protected void setStringValue(String s)
		    throws SAXException
		{
		    try {
			setFeatureProperty(XFFFeatureSetHandler.PROPERTY_XFF_ID, s.trim());
		    } catch (Exception ex) {
			throw new SAXException("Couldn't set property", ex);
		    }
		}
	    } ;
    }

    protected StAXContentHandler getLocationHandler() {
	return new LocationHandlerBase() {
		protected void setLocationValue(Location l) {
		    getFeatureTemplate().location = l;
		}
	    } ;
    }
}
