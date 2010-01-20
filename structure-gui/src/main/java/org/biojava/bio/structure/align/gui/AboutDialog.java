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
 * Created on Jul 21, 2009
 * Created by ap3
 *
 */

package org.biojava.bio.structure.align.gui;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JEditorPane;
import javax.swing.JScrollPane;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.util.ResourceManager;
import org.biojava.bio.structure.align.webstart.BrowserOpener;

public class AboutDialog 
{
   Box vBox;
   public AboutDialog(){
      
   }
   
   public void showDialog(){
      JDialog dialog = new JDialog();

      dialog.setSize(new Dimension(500,650));

      ResourceManager mgr = ResourceManager.getResourceManager("ce");

      String msg = "";
    	  
      msg += mgr.getString("ce.about");
      
      msg += "<b>Currently suported algorithms and version:</b><br>";
      // add the Algorithms and  version nrs.
      
      StructureAlignment[] algorithms = StructureAlignmentFactory.getAllAlgorithms();
      for (StructureAlignment algorithm: algorithms){
    	  msg+="<i>"+algorithm.getAlgorithmName()+"</i> V." +algorithm.getVersion()+"<br>";
      }
      //msg+="<hr>";
      
      JEditorPane txt = new JEditorPane("text/html", msg);
      txt.setEditable(false);

      JScrollPane scroll = new JScrollPane(txt);
      scroll.setSize(new Dimension(300,500));
      vBox= Box.createVerticalBox();
      vBox.add(scroll);
      
      txt.addHyperlinkListener(new HyperlinkListener(){
         
         public void hyperlinkUpdate(HyperlinkEvent e) {
             
             if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
                 String href = e.getDescription();
                 BrowserOpener.showDocument(href);
             }
             if ( e.getEventType() == HyperlinkEvent.EventType.ENTERED) {
                 // change the mouse curor
                 vBox.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
             }
             if (e.getEventType() == HyperlinkEvent.EventType.EXITED) { 
                 vBox.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
             }
         }
     });


      

      JButton close = new JButton("Close");

      close.addActionListener(new ActionListener(){
         public void actionPerformed(ActionEvent event) {
            Object source = event.getSource();

            JButton but = (JButton)source;
            Container parent = but.getParent().getParent().getParent().getParent().getParent().getParent() ;

            JDialog dia = (JDialog) parent;
            dia.dispose();
         }
      });

      Box hBoxb = Box.createHorizontalBox();
      hBoxb.add(Box.createGlue());
      hBoxb.add(close,BorderLayout.EAST);

      vBox.add(hBoxb);

      dialog.getContentPane().add(vBox);
      dialog.setVisible(true);
      
      
   }
}
