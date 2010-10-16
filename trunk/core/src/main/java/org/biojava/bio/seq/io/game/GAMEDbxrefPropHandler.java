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
 
package org.biojava.bio.seq.io.game;

import java.util.List;
import java.util.Vector;

import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.SAXException;

/**
 * Deals with database crossreferences
 *
 * @author David Huen
 */
public class GAMEDbxrefPropHandler 
               extends StAXPropertyHandler
{
  // set up factory method
  public static final StAXHandlerFactory GAME_DBXREF_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new GAMEDbxrefPropHandler(staxenv);
    }
  };

  // the key values in a <db_xref>
  private String XrefDb = null;
  private String DbXrefId = null;

  public class DbXrefElement
  {
    private String XrefDb;
    private String DbXrefId;

    private DbXrefElement(String XrefDb, String DbXrefId)
    {
      this.XrefDb = XrefDb;
      this.DbXrefId = DbXrefId;
    }

    public String getXrefDb()
    {
      return XrefDb;
    }

    public String getDbXrefId()
    {
      return DbXrefId;
    }
  }

  private class XrefDbHandler extends StringElementHandlerBase
  {
    protected void setStringValue(String s)
    {
      XrefDb = s.trim();
    }
  }

  private class DbXrefIdHandler extends StringElementHandlerBase
  {
    protected void setStringValue(String s)
    {
      DbXrefId = s.trim();
    }
  }

  GAMEDbxrefPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);

    setHandlerCharacteristics("dbxref", true); 

    // setup handlers
    super.addHandler(new ElementRecognizer.ByLocalName("xref_db"),
      new StAXHandlerFactory() {
           public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
             return new XrefDbHandler(); }
      }
    );
        
 
    super.addHandler(new ElementRecognizer.ByLocalName("db_xref_id"),
      new StAXHandlerFactory() {
           public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
             return new DbXrefIdHandler(); }
      }
    );
   
  }

/*
  public void setXrefDb(String s)
  {
    XrefDb = s;
  }

  public void setDbXrefId(String s)
  {
    DbXrefId = s;
  }
*/

/*
  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
         throws SAXException
  {
//    System.out.println("GAMEDbxrefPropHandler.startElementHandler entered");
  }
*/

/**
 * when exiting, put the DbXrefElement into the annotation bundle
 */
  public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
                throws SAXException
  {
    if (XrefDb == null || DbXrefId == null) {
      // malformed <db_xref>
      throw new SAXException("Malformed <db_xref> ");
    }

    // create the dbxref List if need be.
    try {

      if (!staxenv.featureTemplate.annotation.containsProperty("dbxref_list")) {
        // create the necessary list element
        staxenv.featureTemplate.annotation.setProperty("dbxref_list", new Vector());
      }

      // stash away the data in it
      List dbxrefList = (List) staxenv.featureTemplate.annotation.getProperty("dbxref_list");
//      System.out.println("DbXrefElement set up with " + XrefDb + " " + DbXrefId);
      dbxrefList.add(new DbXrefElement(XrefDb, DbXrefId));
    }
    catch (ChangeVetoException cve) {
      System.err.println("GAMEDbxrefPropHandler: change vetoed.");
    }
  }
}
