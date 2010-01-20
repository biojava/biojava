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

package org.biojava.bio.structure.align.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;

public class MyExportListener implements ActionListener{

   StructureAlignmentJmol parent;
   MyExportListener(StructureAlignmentJmol parent){
      this.parent = parent;
   }
   public void actionPerformed(ActionEvent arg0)
   {
      final JFileChooser fc = new JFileChooser();

      int returnVal = fc.showSaveDialog(null);

      if (returnVal == JFileChooser.APPROVE_OPTION) {
         File file = fc.getSelectedFile();
         //This is where a real application would open the file.
         System.out.println("Exporting PDB file to: " + file.getName());

         Structure s = parent.getStructure();

         try {
            PrintWriter pw = new PrintWriter(new FileWriter(file));                 
            pw.println(s.toPDB());                   
            pw.close();
         } catch (IOException e){
        	 JOptionPane.showMessageDialog(null,"Could not export file. Excetion: " + e.getMessage());
         }


      } else {
         System.out.println("Export command cancelled by user.");
      }


   }
}
