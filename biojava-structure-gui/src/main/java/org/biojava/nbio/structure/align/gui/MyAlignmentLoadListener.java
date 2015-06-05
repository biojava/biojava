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
 */
package org.biojava.nbio.structure.align.gui;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.align.webstart.WebStartMain;
import org.biojava.nbio.structure.align.xml.AFPChainXMLParser;
import org.biojava.nbio.core.util.InputStreamProvider;

import javax.swing.*;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;


public class MyAlignmentLoadListener implements ActionListener {

   AbstractAlignmentJmol jmol;
   public MyAlignmentLoadListener(AbstractAlignmentJmol jmol){
      this.jmol = jmol;
   }
   @Override
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
            AtomCache cache = new AtomCache(config.getPdbFilePath(),config.getCacheFilePath());

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
