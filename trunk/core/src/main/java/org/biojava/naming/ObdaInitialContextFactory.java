package org.biojava.naming;

import java.util.Hashtable;

import javax.naming.Context;
import javax.naming.Name;
import javax.naming.NamingException;
import javax.naming.directory.BasicAttributes;
import javax.naming.spi.InitialContextFactory;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.utils.ClassTools;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.DefaultHandler;

/**
 *
 *
 * @author Matthew Pocock
 */
public class ObdaInitialContextFactory
        implements InitialContextFactory
{
  private static final String CORE = "obda/naming/core.xml";

  public Context getInitialContext(Hashtable environment)
          throws NamingException
  {
    try {
      InputSource iSource = new InputSource(ClassTools.getClassLoader(ObdaInitialContextFactory.class).getResourceAsStream(CORE));
      SAXParserFactory spf = SAXParserFactory.newInstance();
      spf.setValidating(false);
      spf.setNamespaceAware(true);
      XMLReader reader = spf.newSAXParser().getXMLReader();
      ObdaHandler handler = new ObdaHandler(environment);
      reader.setContentHandler(handler);
      reader.parse(iSource);
      return handler.getRoot();
    } catch (Exception e) {
      throw new Error(e);
    }
  }

  private class ObdaHandler
          extends DefaultHandler
  {
    Hashtable env;
    ObdaContext root;
    ObdaContext current;
    StringBuffer description = null;

    ObdaHandler(Hashtable env)
    {
      this.env = env;
    }

    public ObdaContext getRoot()
    {
      return root;
    }

    public void startDocument()
            throws SAXException
    {
      root = new ObdaContext(null, null,
                             new Hashtable(), env, new BasicAttributes());
    }

    public void startElement(String uri, String localName,
                             String qName, Attributes attributes)
            throws SAXException
    {
      if(qName.equals("directory")) {
        // do directory things
      } else if(qName.equals("urn")) {
        // unpack the URN, break it into words, make the contexts (if needed)
        try {
          Name name = ObdaUriParser.getInstance()
                  .parse(attributes.getValue("name"));
          ObdaContext ctxt = root;
          for(int i = 0; i < name.size(); i++) {
            ctxt = resolve(ctxt, name.get(i));
          }
          current = ctxt;
        } catch (NamingException e) {
          throw new SAXException(e);
        }
      } else if(qName.equals("description")) {
        // add the description to the current context - we will first have to
        // collect all the description text together and then set it on
        // end element
        description = new StringBuffer();
      }
    }

    public void endElement(String uri, String localName, String qName)
            throws SAXException
    {
      if(qName.equals("description") && description != null) {
        current.getAttrs().put("description", description.toString());
        description = null;
      }
    }

    public void characters(char ch[], int start, int length)
            throws SAXException
    {
      if(description != null) {
        description.append(ch, start, length);
      }
    }

    private ObdaContext resolve(ObdaContext parent, String name)
    {
      Hashtable bindings = parent.getBindings();
      ObdaContext child = (ObdaContext) bindings.get(name);

      if(child == null) {
        child = new ObdaContext(
                parent, name,
                new Hashtable(), new Hashtable(), new BasicAttributes());
        bindings.put(name, child);
      }

      return child;
    }
  }
}
