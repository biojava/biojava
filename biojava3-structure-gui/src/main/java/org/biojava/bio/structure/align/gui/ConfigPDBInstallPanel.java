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

package org.biojava.bio.structure.align.gui;

import java.awt.BorderLayout;
import java.awt.Dimension;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;

import javax.swing.Action;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;


import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.align.webstart.WebStartMain;
import org.biojava.bio.structure.gui.util.PDBUploadPanel;

public class ConfigPDBInstallPanel extends JPanel
{

   /**
    * 
    */
   private static final long serialVersionUID = -1055193854675583808L;

   /** the system property PDB_DIR can be used to configure the 
    * default location for PDB files.
    */
   public static final String PDB_DIR = "PDB_DIR";
   
   JCheckBox pdbSplit;
   JCheckBox fromFtp;
   JComboBox fileType;
   
   JTextField pdbDir;

   private static ConfigPDBInstallPanel instance = new ConfigPDBInstallPanel();
   
   static JDialog dialog;
   
   private ConfigPDBInstallPanel(){
      
      UserConfiguration   config = WebStartMain.getWebStartConfig();

      fileType = PDBUploadPanel.getFileFormatSelect();

      Box vBox = Box.createVerticalBox();
      Box hBox = Box.createHorizontalBox();

      pdbDir = new JTextField(20);
      
      String conf = System.getProperty(PDB_DIR);
      if ( conf != null){
         pdbDir.setText(conf);
      }

      JLabel l01 = new JLabel("Directory containing PDBs");
      //panel.add(l01);
      //panel.add(f);

      hBox.add(l01);
      hBox.add(pdbDir);
      vBox.add(hBox);

      if ( config != null) {
         pdbDir.setText( config.getPdbFilePath() );
      }

      Action action = new ChooseDirAction(pdbDir, config);

      JButton chooser = new JButton(action);
      hBox.add(chooser);

      Box hBox2 = Box.createHorizontalBox();

      JLabel label = new JLabel("Files in split directories:");

      pdbSplit = new JCheckBox();
      pdbSplit.setMnemonic(KeyEvent.VK_S);

      pdbSplit.setSelected(true);
      if ( config != null){
         pdbSplit.setSelected(config.isSplit());
      }

      hBox2.add(Box.createGlue());
      hBox2.add(label);
      hBox2.add(pdbSplit);

      vBox.add(hBox2);

      Box hBox3 = Box.createHorizontalBox();
      JLabel label2 = new JLabel("Fetch missing PDBs from ftp site:");
      fromFtp = new JCheckBox();
      fromFtp.setMnemonic(KeyEvent.VK_F);
      fromFtp.setSelected(true);
      if ( config != null)
         fromFtp.setSelected(config.getAutoFetch());

      JLabel ftype = new JLabel("File format:");

      hBox3.add(Box.createGlue());
      hBox3.add(ftype);
      hBox3.add(fileType);
      hBox3.add(Box.createGlue());
      hBox3.add(label2);
      hBox3.add(fromFtp);

      vBox.add(hBox3);
      this.add(vBox);

   }
   
   public static void showDialog(){
      if ( dialog != null) {
        
         dialog.setVisible(true);
         return;
      }
      
      dialog = new JDialog();
      dialog.setSize(new Dimension(600,300));
      Box vBox = Box.createVerticalBox();
      
      vBox.add(instance);

      
      
      UIManager.LookAndFeelInfo lookAndFeels[] = UIManager.getInstalledLookAndFeels();
      JPanel panel = new JPanel();
      
      for(int i = 0; i < lookAndFeels.length; i++){
        JButton button = new JButton(lookAndFeels[i].getName());
        button.addActionListener(new MyAction(dialog));
        panel.add(button);
      }
      
      JTabbedPane jTab = new JTabbedPane();
      jTab.addTab("Select User Interface", null, panel, "Select the Look and Feel of the application.");
           
      vBox.add(jTab);
      vBox.add(Box.createGlue());
      
      JButton apply = new JButton("Apply");
      apply.addActionListener(new ActionListener(){
         public void actionPerformed(ActionEvent event) {
            instance.applyValues();
            dialog.dispose();
         }
      });
      
      JButton close = new JButton("Cancel");

      close.addActionListener(new ActionListener(){
         public void actionPerformed(ActionEvent event) {
           dialog.dispose();
         }
      });

      Box hBoxb = Box.createHorizontalBox();
      hBoxb.add(Box.createGlue());
      hBoxb.add(close,BorderLayout.EAST);
      hBoxb.add(apply,BorderLayout.EAST);
      vBox.add(hBoxb);

      dialog.getContentPane().add(vBox);
      dialog.setVisible(true);
      
   }

   protected void applyValues()
   {
      UserConfiguration config = WebStartMain.getWebStartConfig();
      
      String dir = pdbDir.getText();
      config.setPdbFilePath(dir);
      boolean isSplit = pdbSplit.isSelected();
      config.setSplit(isSplit);
      boolean fromFtpF = fromFtp.isSelected();
      config.setAutoFetch(fromFtpF);
      String fileFormat = (String)fileType.getSelectedItem();
      config.setFileFormat(fileFormat);
      
      // and now persist..
      WebStartMain.persistConfig(config);
      
   }

   public JCheckBox getFromFtp()
   {
      return fromFtp;
   }
   public void setFromFtp(JCheckBox fromFtp)
   {
      System.out.println(fromFtp);
      this.fromFtp = fromFtp;
   }
   public JCheckBox getPdbSplit()
   {
      return pdbSplit;
   }
   public void setPdbSplit(JCheckBox pdbSplit)
   {
      this.pdbSplit = pdbSplit;
   }
   public JTextField getPDBDirField(){
      return pdbDir;
   }
   public void setPDBDirField(JTextField dir){
      pdbDir = dir;
   }


   

}

class MyAction implements ActionListener{
   JDialog dialog;
   
   public MyAction(JDialog dialog) {
      this.dialog= dialog;
   }
   
   
   @SuppressWarnings("unused")
public void actionPerformed(ActionEvent ae){
     Object EventSource = ae.getSource();
     String lookAndFeelClassName = null;
     UIManager.LookAndFeelInfo looks[] = UIManager.getInstalledLookAndFeels();
     for(int i = 0; i < looks.length; i++){
       if(ae.getActionCommand().equals(looks[i].getName())){
         lookAndFeelClassName = looks[i].getClassName();
         break;
       }
     }
     try{
       UIManager.setLookAndFeel(lookAndFeelClassName);
      
       SwingUtilities.updateComponentTreeUI(dialog);
       
       SwingUtilities.updateComponentTreeUI( AlignmentGui.getInstanceNoVisibilityChange());
     }
     catch(Exception e){
       JOptionPane.showMessageDialog(dialog, "Can't change look and feel:" 
+ lookAndFeelClassName, "Invalid PLAF", JOptionPane.ERROR_MESSAGE);
     }
   }
 }
