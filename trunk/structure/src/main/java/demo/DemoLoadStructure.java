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
 * Created on Jan 27, 2010
 * Author: Andreas Prlic 
 *
 */

package demo;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.utils.io.InputStreamProvider;

public class DemoLoadStructure
{

   public static void main(String[] args){




      DemoLoadStructure demo  = new DemoLoadStructure();

      demo.basicLoad();

      demo.loadStructureFromCache();
   }

   public void basicLoad(){
      try {

         PDBFileReader reader = new PDBFileReader();

         // the path to the local PDB installation
         reader.setPath("/tmp");

         // are all files in one directory, or are the files split,
         // as on the PDB ftp servers?
         reader.setPdbDirectorySplit(true);

         // should a missing PDB id be fetched automatically from the FTP servers?
         reader.setAutoFetch(true);

         // should the ATOM and SEQRES residues be aligned when creating the internal data model?
         reader.setAlignSeqRes(false);

         // should secondary structure get parsed from the file
         reader.setParseSecStruc(false);

         Structure structure = reader.getStructureById("4hhb");

         System.out.println(structure);
      } catch (Exception e){
         e.printStackTrace();
      }

   }

   public void loadStructureFromCache(){
      String pdbId = "4hhb";
      String chainName = "4hhb.A";
      String entityName = "4hhb:0";

      // split PDB file installation?
      boolean isPdbDirectorySplit = true;

      String pdbFilePath = "/tmp/";

      // we can set a flag if the file should be cached in memory
      // this will enhance IO massively if the same files have to be accessed over and over again.
      // since this is a soft cache, no danger of memory leak
      // this is actually not necessary to provide, since the default is "true" if the AtomCache is being used.
      System.setProperty(InputStreamProvider.CACHE_PROPERTY, "true");

      AtomCache cache = new AtomCache(pdbFilePath,isPdbDirectorySplit);

      try {
         System.out.println("======================");
         Structure s = cache.getStructure(pdbId);

         System.out.println("Full Structure:" + s);

         Atom[] ca = cache.getAtoms(chainName);
         System.out.println("got " + ca.length + " CA atoms");

         Structure firstEntity = cache.getStructure(entityName);
         System.out.println("First entity: " + firstEntity);

      } catch (Exception e){
         e.printStackTrace();
      }

   }
}
