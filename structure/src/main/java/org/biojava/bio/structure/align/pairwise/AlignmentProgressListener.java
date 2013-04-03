/*
 *                    PDB web development code
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
 *
 * Created on May 23, 2009
 * Created by Andreas Prlic
 *
 */

package org.biojava.bio.structure.align.pairwise;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.List;

import javax.swing.JFrame;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.helper.JointFragments;
import org.biojava.bio.structure.jama.Matrix;


/** A class to listen to progress of the structure alignmnent calcualtions
 *
 * @author andreas
 *
 */
public class AlignmentProgressListener
{

   String n1;
   String n2;
   int l1;
   int l2;

   Atom[] ca1 ;
   Atom[] ca2;

   public AlignmentProgressListener(){

   }

   public void startingAlignment(String name1, Atom[] ca1, String name2, Atom[] ca2){
      n1 = name1;
      n2 = name2;
      l1 = ca1.length;
      l2 = ca2.length;
      this.ca1 = ca1;
      this.ca2 = ca2;
   }

   public void calculatedFragmentPairs(List<FragmentPair> fragments){
      System.out.println("got: " + fragments.size() + " fragments");

      String title = "Initial FragmentPairs for:" +  n1 + "("+l1+")"+ " vs. " + n2 + " ("+l2+")";
     // ScaleableMatrixPanel panel = new ScaleableMatrixPanel();


      Matrix m = new Matrix(l1,l2,99);

      for (FragmentPair p : fragments){
         for (FragmentPair pair2: fragments){


            //Atom v2 = tmpfidx[j].getCenter2();
            Atom v1 = p.getUnitv();
            Atom v2 = pair2.getUnitv();
            //System.out.println("v1: "+v1);
            //System.out.println("v2: "+v2);
            try{
               double dist = Calc.getDistance(v1,v2);
               for (int i =0 ; i < p.getLength(); i++){
                  int p1 = p.getPos1();
                  int p2 = p.getPos2();
                  m.set(p1+i,p2+i,dist);
               }
               for (int i =0 ; i < pair2.getLength(); i++){
                  int p1 = pair2.getPos1();
                  int p2 = pair2.getPos2();
                  m.set(p1+i,p2+i,dist);
               }

            } catch (StructureException e){
               e.printStackTrace();
            }
         }

      }
    //  panel.setMatrix(m);
      JFrame frame = new JFrame();

      frame.setTitle(title);

      frame.addWindowListener(new WindowAdapter(){
         public void windowClosing(WindowEvent e){
            JFrame f = (JFrame) e.getSource();
            f.setVisible(false);
            f.dispose();
         }



      });
      //frame.getContentPane().add(panel);
      frame.pack();
      frame.setVisible(true);



   }


   public void jointFragments(JointFragments[] fragments){
      System.out.println("numberof Joint fragments: " + fragments.length);

      String title = "JointFragment for:" +  n1 + "("+l1+")"+ " vs. " + n2 + " ("+l2+")";
     // ScaleableMatrixPanel panel = new ScaleableMatrixPanel();

      Matrix m = new Matrix(l1,l2,99);

      for (JointFragments p : fragments){
         for (int[] idx : p.getIdxlist() ){
            m.set(idx[0],idx[1],p.getRms());
         }
      }
     // panel.setMatrix(m);
      JFrame frame = new JFrame();

      frame.setTitle(title);

      frame.addWindowListener(new WindowAdapter(){
         public void windowClosing(WindowEvent e){
            JFrame f = (JFrame) e.getSource();
            f.setVisible(false);
            f.dispose();
         }



      });
     // frame.getContentPane().add(panel);
      frame.pack();
      frame.setVisible(true);




      for (JointFragments f : fragments){
         System.out.println(f);
      }
   }
}



