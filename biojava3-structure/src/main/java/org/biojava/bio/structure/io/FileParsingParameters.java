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

import org.biojava.bio.structure.AminoAcid;

/** A class that configures parameters that can be sent to the PDB file parsers
 * 
 * <ul>
 * <li> {@link #setParseCAOnly(boolean)} - parse only the Atom records for C-alpha atoms</li>
 * <li> {@link #setParseSecStruc(boolean)} - a flag if the secondary structure information from the PDB file (author's assignment) should be parsed.
 *      If true the assignment can be accessed through {@link AminoAcid}.getSecStruc(); </li>
 * <li> {@link #setAlignSeqRes(boolean)} - should the AminoAcid sequences from the SEQRES
 *      and ATOM records of a PDB file be aligned? (default:yes)</li>
 * </ul>
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



   /** Should the chemical component information be automatically be downloade from the web?
    * If set to false, a limited set of ChemComps is being used.
    * @return flag if the data should be loaded
    */
   public boolean isLoadChemCompInfo()
   {
      return loadChemCompInfo;
   }

   /**  Sets if chemical component defintions should be loaded from the web
    * 
    * @param loadChemCompInfo flag
    */
   public void setLoadChemCompInfo(boolean loadChemCompInfo)
   {

      if ( loadChemCompInfo)
         System.setProperty(PDBFileReader.LOAD_CHEM_COMP_PROPERTY, "true");
      this.loadChemCompInfo = loadChemCompInfo;

   }

   /** Parse only the PDB file header out of the files
    * 
    * @return flag
    */
   public boolean isHeaderOnly()
   {
      return headerOnly;
   }

   /** Parse only the PDB file header out of the files
    * 
    * @param headerOnly flag
    */
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
