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
package org.biojava.bio.seq.io.agave;
import java.util.ArrayList;
import java.util.List;

import org.biojava.utils.ChangeVetoException;

/**
 * Deals with alternate sequence IDs
 * 
 * @author Hanning Ni     Doubletwist Inc
 */
public class AGAVEAltIdsPropHandler
               extends StAXPropertyHandler implements AGAVEDbIdCallbackItf
{
  // set up factory method
  public static final StAXHandlerFactory AGAVE_ALT_IDS_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEAltIdsPropHandler(staxenv);
    }
  };


  AGAVEAltIdsPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);

    setHandlerCharacteristics("alt_ids", true);

    super.addHandler(new ElementRecognizer.ByLocalName("db_id"),
         AGAVEDbIdPropHandler.AGAVE_DBID_PROP_HANDLER_FACTORY);

  }

   public void addDbId(AGAVEDbId db_id)
   {
      try{
         Object ob = UtilHelper.getProperty( staxenv.featureTemplate.annotation, "alt_ids");
         if( ob != null )
             ((List)ob).add( db_id ) ;
         else
         {
             List kws = new ArrayList(1) ;
             kws.add( db_id ) ;
             staxenv.featureTemplate.annotation.setProperty("alt_ids", kws);
         }
      }catch (ChangeVetoException cve) {
        cve.printStackTrace();
      }
   }
}

