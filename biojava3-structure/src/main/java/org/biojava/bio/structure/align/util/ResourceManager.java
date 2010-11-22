/*
 *                  BioJava development code
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
 * Created on Aug 4, 2006
 *
 */
package org.biojava.bio.structure.align.util;

import java.util.MissingResourceException;
import java.util.ResourceBundle;



/** A class that manages the Strings that are defined in the spice.properties file. 
 * This will be usefull for internationalisation. 
 * 
 * TODO: provide .properties files for other locales.
 * e.g. jfatcat_de_DE.properties, etc.
 * 
 * @author Andreas Prlic
 * @since 1:43:04 PM
 * @version %I% %G%
 */
public class ResourceManager {

   private String BUNDLE_NAME ;

   private ResourceBundle RESOURCE_BUNDLE ;


   public ResourceManager() {
      BUNDLE_NAME =  "jfatcat"; //$NON-NLS-1$;
      RESOURCE_BUNDLE = ResourceBundle.getBundle(BUNDLE_NAME);
   }

   private ResourceManager(String bundleName){
      try {
         RESOURCE_BUNDLE = ResourceBundle.getBundle(bundleName);
      } catch(Exception e){
         e.printStackTrace();
      }

   }
   public static ResourceManager getResourceManager(String bundleName){
      return new ResourceManager(bundleName);
   }

   public String getString(String key) {

      try {
         return RESOURCE_BUNDLE.getString(key);
      } catch (MissingResourceException e) {
         System.err.println(e.getMessage());
         return '!' + key + '!';
      }
   }
}
