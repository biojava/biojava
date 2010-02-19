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
 * Created on Feb 8, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.webstart;

import java.lang.reflect.Method;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.model.AFPChain;


/** we can't release jFatCat source code to the public as of yet. As such this provides a wrapper
 * 
 * @author Andreas Prlic
 *
 */
public class JFatCatProxy
{

   public static final String afpTwisterClassName = "org.rcsb.fatcat.calc.AFPTwister"; 

   public static final String fatCatAlignerClassName = "org.rcsb.fatcat.calc.FatCatAligner";

   public static final String fatCatRigidClassName = "org.rcsb.fatcat.FatCatRigid";

   Object fatCatRigid;

   public JFatCatProxy(){
      fatCatRigid = null;
   }

   public void setStructureAlignment(Object fatCatRigid)
   {

      this.fatCatRigid = fatCatRigid;

   }

   public  Group[] twistGroups(AFPChain afpChain, Atom[] ca1, Atom[] ca2) 
   {


      // does the following:
      //Group[] twistedGroups = AFPTwister.twistOptimized(afpChain,ca1,ca2);

      //FatCatAligner aligner =  fatCat.getFatCatAligner();
      //aligner.setTwistedGroups(twistedGroups);
      try {
         Class afpTwister = Class.forName(afpTwisterClassName);

         Method m = afpTwister.getMethod("twistOptimized", new Class[] { AFPChain.class, Atom[].class, Atom[].class});

         Group[] twistedGroups = (Group[]) m.invoke(null,afpChain,ca1,ca2);


         if ( afpChain.getAlgorithmName().startsWith("jFatCat")) {
            Class fatCatRigidC = Class.forName(fatCatRigidClassName);

            Method getFatCatAligner =   fatCatRigidC.getMethod("getFatCatAligner", new Class[]{});

            Class fatCatAligner = Class.forName(fatCatAlignerClassName);

            Object fatCatAlignerInstance = getFatCatAligner.invoke(fatCatRigid, null);

            Method setTwistedGroups = fatCatAligner.getMethod("setTwistedGroups",new Class[]{Group[].class});

            setTwistedGroups.invoke(fatCatAlignerInstance,  (Object) twistedGroups);
         }

         // we only have data for the optimized alignment so far...

         return twistedGroups;

      }  catch (Exception e){
         System.err.println("jfatcat.jar in classpath?");
         e.printStackTrace();

      }
      return null;

   }

}
