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
 * Created on Jul 28, 2009
 * Created by ap3
 *
 */

package org.biojava.nbio.structure.align.gui;

import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.gui.ScaleableMatrixPanel;
import org.biojava.nbio.structure.jama.Matrix;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

public class MyDistMaxListener
implements ActionListener{
   AFPChain parent;
   public MyDistMaxListener(AFPChain parent){
      this.parent = parent;
   }
   @Override
public void actionPerformed(ActionEvent arg0)
   {

      System.out.println("show distance matrices");

      if ( parent == null) {
         System.err.println("Not displaying any alignment currently!");
         return;
      }
      if (parent.getDisTable1()!=null) showMatrix(parent.getDisTable1(), "Internal distances for Structure 1");
      if (parent.getDisTable2()!=null) showMatrix(parent.getDisTable2(), "Internal distances for Structure 2");

   }

   private void showMatrix(Matrix m, String title){
      ScaleableMatrixPanel smp = new ScaleableMatrixPanel();
      JFrame frame = new JFrame(title);
      frame.addWindowListener(new WindowAdapter(){
         @Override
		public void windowClosing(WindowEvent e){
            JFrame f = (JFrame) e.getSource();
            f.setVisible(false);
            f.dispose();
         }



      });

      smp.setMatrix(m);

      frame.getContentPane().add(smp);

      frame.pack();
      frame.setVisible(true);
   }


}

