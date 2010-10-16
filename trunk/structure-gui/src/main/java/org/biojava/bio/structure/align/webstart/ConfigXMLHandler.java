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
 *
 *      http://www.biojava.org/

 * @author Andreas Prlic
 *
 */
package org.biojava.bio.structure.align.webstart;


import org.biojava.bio.structure.align.util.UserConfiguration;
import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.Attributes ;



/**
 * XML content handler for serialisation of RegistryConfiguration class
 */
public class ConfigXMLHandler extends DefaultHandler {

   UserConfiguration config ;

   /**
    * 
    */
   public ConfigXMLHandler() {
      super();

      config         = new UserConfiguration();
   }

   public void startElement (String uri, String name, String qName, Attributes atts){
      //System.out.println("new element >" + name + "< >" + qName+"<" + uri);
      if ( qName.equals("PDBFILEPATH")){

         String path = atts.getValue("path");
         // default path is system tmp...
         if ( path != null)
            config.setPdbFilePath(path);

         String isSplit = atts.getValue("split");
         config.setSplit(true);
         if ( isSplit != null)     {   	 
            if ( isSplit.equals("false"))
               config.setSplit(false);
         }

         String autoFetch = atts.getValue("autoFetch");
         config.setAutoFetch(true);
         if ( autoFetch != null){
            if ( autoFetch.equals("false"))
               config.setAutoFetch(false);
         }

         String fileFormat = atts.getValue("fileFormat");
         config.setFileFormat(UserConfiguration.PDB_FORMAT);
         if ( fileFormat != null) {
            if ( fileFormat.equals(UserConfiguration.MMCIF_FORMAT))
               config.setFileFormat(UserConfiguration.MMCIF_FORMAT);
         }

      }
   }







   public UserConfiguration getConfig() {
      return config ;
   }

}
