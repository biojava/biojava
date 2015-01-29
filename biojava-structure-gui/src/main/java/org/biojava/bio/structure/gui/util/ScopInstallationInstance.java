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
 * Created on Jun 30, 2010
 * Author: ap3 
 *
 */

package org.biojava.bio.structure.gui.util;

import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.align.webstart.WebStartMain;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopInstallation;

public class ScopInstallationInstance
{
   
   
   static ScopInstallationInstance me = new ScopInstallationInstance();
   ScopDatabase install;
   private ScopInstallationInstance(){
      UserConfiguration config = WebStartMain.getWebStartConfig();
      String cacheLocation = config.getPdbFilePath();

       install = new ScopInstallation(cacheLocation);
   }

   
   public static ScopInstallationInstance getInstance(){
      return me;
   }
   public  ScopDatabase getSCOP(){
      return install;
   }
}
