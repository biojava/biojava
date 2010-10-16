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
 * Created on Feb 2, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.dasobert;


import org.biojava.bio.program.das.dasalignment.Alignment;
import org.biojava.dasobert.das.AlignmentParameters;
import org.biojava.dasobert.das.AlignmentThread;
import org.biojava.dasobert.dasregistry.Das1Source;
import org.biojava.dasobert.eventmodel.AlignmentEvent;
import org.biojava.dasobert.eventmodel.AlignmentListener;

import junit.framework.TestCase;

public class TestSisyphusServer extends TestCase
{

   public void testServer(){
      
      AlignmentParameters params = new AlignmentParameters();
      
      Das1Source dasSource = new Das1Source();
      
      dasSource.setUrl("http://sisyphus.mrc-cpe.cam.ac.uk/sisyphus/das/alignments/");
      
      params.setDasSource(dasSource);
      params.setQuery("AL10057468");
      
      AlignmentThread t = new AlignmentThread(params);
      
      AlignmentListener li = new AlignmentListener()
      {
         
         public void noAlignmentFound(AlignmentEvent e)
         {
            // TODO Auto-generated method stub
            
         }
         
         public void newAlignment(AlignmentEvent e)
         {
           Alignment a = e.getAlignment();
            System.out.println(a);
            
         }
         
         public void clearAlignment()
         {
            // TODO Auto-generated method stub
            
         }
      };
      
      t.addAlignmentListener(li);
      
      t.run();
      
   }
   
}
