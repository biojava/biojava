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
 * Created on Apr 6, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.webstart;

import javax.swing.UIManager;
import javax.swing.UIManager.LookAndFeelInfo;


public class AligUIManager
{

   public static void setLookAndFeel(){
      try {
   
         
         String system = UIManager.getSystemLookAndFeelClassName();
         if ( system != null) {
            //System.out.println("setting look and feel to " + system);
            UIManager.setLookAndFeel(system);
            
         }
         
         //System.out.println("Installed Look And Feels:");
         LookAndFeelInfo[] feels = UIManager.getInstalledLookAndFeels();
         
         if ( feels != null){
            //for ( LookAndFeelInfo info: feels){
               //System.out.println(info.getName() + " " + info.getClassName());
           // }
         }
         
         
         //System.out.println("Auxiliary Look And Feels:");
        // LookAndFeel[] looks = UIManager.getAuxiliaryLookAndFeels();
         //printLookAndFeel(looks);
         
         
        

      } catch ( Exception e ) {
         e.printStackTrace();
      }

   }

//   private static void printLookAndFeel(LookAndFeel[] looks)
//   {
//
//      if ( looks != null){
//         System.out.println("got " + looks.length + " lookAndFeels");
//         for (LookAndFeel laf : looks){
//            System.out.println(laf.getDescription());
//         }
//      } else {
//         System.out.println("No other LookAndFeels found.");
//      }
//      
//   }
}
