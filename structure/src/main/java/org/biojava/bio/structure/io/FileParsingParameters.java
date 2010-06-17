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
 * Created on Jun 16, 2010
 * Author: ap3 
 *
 */

package org.biojava.bio.structure.io;

import java.io.Serializable;

/** A class that configures parameters that can be sent to the PDB file parsers
 * 
 * @author Andreas Prlic
 *
 */
public class FileParsingParameters implements Serializable
{

   /**
    * 
    */
   private static final long serialVersionUID = 5878292315163939027L;

   /** flag to detect if the secondary structure info should be read
    * 
    */
   boolean parseSecStruc;

   /** Flag to control if SEQRES and ATOM records should be aligned
    * 
    */
   boolean alignSeqRes;


   /** Flag to control if the chemical component info should be downloaded while parsing the files. (files will be cached).
    * 
    */
   boolean loadChemCompInfo;

   /** Set the flag to only read in Ca atoms - this is useful for parsing large structures like 1htq.
    *
    */
   boolean parseCAOnly;

   /** Flag to parse header only
    * 
    */
   boolean headerOnly;


   public FileParsingParameters(){
      setDefault();
   }

   public void setDefault(){

      parseSecStruc = false;

      // by default we now do NOT align Atom and SeqRes records
      alignSeqRes   = false;
      parseCAOnly = false;

      // don't download ChemComp dictionary by default.
      loadChemCompInfo = false;
      headerOnly = false;

   }

   /** is secondary structure assignment being parsed from the file?
    * default is null
    * @return boolean if HELIX STRAND and TURN fields are being parsed
    */
   public boolean isParseSecStruc() {
      return parseSecStruc;
   }

   /** a flag to tell the parser to parse the Author's secondary structure assignment from the file
    * default is set to false, i.e. do NOT parse.
    * @param parseSecStruc if HELIX STRAND and TURN fields are being parsed
    */
   public void setParseSecStruc(boolean parseSecStruc) {
      this.parseSecStruc = parseSecStruc;
   }




   public boolean isLoadChemCompInfo()
   {
      return loadChemCompInfo;
   }

   public void setLoadChemCompInfo(boolean loadChemCompInfo)
   {

      if ( loadChemCompInfo)
         System.setProperty(PDBFileReader.LOAD_CHEM_COMP_PROPERTY, "true");
      this.loadChemCompInfo = loadChemCompInfo;

   }


   public boolean isHeaderOnly()
   {
      return headerOnly;
   }

   public void setHeaderOnly(boolean headerOnly)
   {
      this.headerOnly = headerOnly;
   }

   /** the flag if only the C-alpha atoms of the structure should be parsed.
    *
    * @return the flag
    */
   public boolean isParseCAOnly() {
      return parseCAOnly;
   }
   /** the flag if only the C-alpha atoms of the structure should be parsed.
    *
    * @param parseCAOnly boolean flag to enable or disable C-alpha only parsing
    */
   public void setParseCAOnly(boolean parseCAOnly) {
      this.parseCAOnly = parseCAOnly;
   }



   /** Flag if the SEQRES amino acids should be aligned with the ATOM amino acids.
    *
    * @return flag if SEQRES - ATOM amino acids alignment is enabled
    */
   public boolean isAlignSeqRes() {
      return alignSeqRes;
   }



   /** define if the SEQRES in the structure should be aligned with the ATOM records
    * if yes, the AminoAcids in structure.getSeqRes will have the coordinates set.
    * @param alignSeqRes
    */
   public void setAlignSeqRes(boolean alignSeqRes) {
      this.alignSeqRes = alignSeqRes;
   }





}
