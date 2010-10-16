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

package demo;

import org.biojava.bio.structure.align.webstart.WebStartMain;

public class WebStartClientDemo
{

   public static void main(String[] args){
     
      
      //client.main(new String[]{"fatcat", "3BMV.A","2GUY.A", "http://pdb114.rcsb.org:8080/jfatcatserver/align/"});
      //client.main(new String[]{"fatcat", "2GUY.A","3BMV.A"});
      
      //client.main( new String[]{"fatcat", "1EXQ.A","1EX4.B","http://pdb114.rcsb.org:8080/jfatcatserver/align/"} );
      
      //WebStartMain.main( new String[]{"fatcat_flexible", "1cdg.A", "1tim.B"} );
      //WebStartMain.main( new String[]{"ce", "1tim.B", "1cdg.A"} );
      //WebStartMain.main( new String[]{"ce", "1cdg.A", "1tim.B"} );
      //WebStartMain.main( new String[]{"ce_cp", "1vhr.A","2ihb.A"} );
      //WebStartMain.main( new String[]{"fatcat", "2BC3.B","1SWG.D"} );
      WebStartMain.main(new String[]{"fatcat","1P80.D","2IUF.E"});
      //WebStartMain.main(new String[]{"fatcat","1O08.A","1FEZ.A"});
      
      
   }
}
