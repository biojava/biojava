package org.biojava.bio.structure.align.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;

import java.io.InputStream;
import java.io.InputStreamReader;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.align.webstart.WebStartMain;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;
import org.biojava3.core.util.InputStreamProvider;


public class MyAlignmentLoadListener implements ActionListener {

   StructureAlignmentJmol jmol;
   public MyAlignmentLoadListener(StructureAlignmentJmol jmol){
      this.jmol = jmol;
   }
   public void actionPerformed(ActionEvent evt) {

      final JFileChooser fc = new JFileChooser();

      //					In response to a button click:
      int returnVal = fc.showOpenDialog(null);

      if ( returnVal == JFileChooser.APPROVE_OPTION) {

         File file = fc.getSelectedFile();

         try {
            InputStreamProvider ip = new InputStreamProvider();
            InputStream stream = ip.getInputStream(file);
            BufferedReader in = new BufferedReader( new InputStreamReader(stream));


            StringBuffer xml = new StringBuffer();
            String str;
            while ((str = in.readLine()) != null) {
               xml.append(str);
            }
            in.close();

            AFPChain[] afps = AFPChainXMLParser.parseMultiXML(xml.toString());
            AFPChain afpChain = afps[0];

            UserConfiguration config = WebStartMain.getWebStartConfig();
            AtomCache cache = new AtomCache(config.getPdbFilePath(),config.isSplit());

            Atom[] ca1 = cache.getAtoms(afpChain.getName1());
            Atom[] ca2 = cache.getAtoms(afpChain.getName2());

            AFPChainXMLParser.rebuildAFPChain(afpChain, ca1, ca2);

            //Chain c1 = ca1[0].getParent().getParent();
            //Chain c2 = ca2[0].getParent().getParent();


            //StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(afpChain.getAlgorithmName());
            StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(afpChain, ca1, ca2); 

            //String result = afpChain.toFatcat(ca1, ca2);

            //String rot = afpChain.toRotMat();
            DisplayAFP.showAlignmentImage(afpChain, ca1,ca2,jmol);



         } catch (Exception e){
            e.printStackTrace();
            JOptionPane.showMessageDialog(null,"Could not load alignment file. Exception: " + e.getMessage());
         }

      }

   }

}
