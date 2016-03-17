package org.biojava.nbio.structure.io.mmtf;


import org.biojava.nbio.structure.Structure;
import org.rcsb.mmtf.decoder.DecodeStructure;
import org.rcsb.mmtf.decoder.ParsingParams;

public class ParseUsingBioJava {

   /**
    * Utility function to get a biojava structure from a byte array.
    * @param inputByteArray Must be uncompressed (i.e. with entropy compression methods like gzip)
    * @param parsingParams
    * @return
    */
   public Structure getBiojavaStruct(byte[] inputByteArray, ParsingParams parsingParams) {
     // Make the decoder
     BioJavaStructureDecoder biojavaStructureDecoder = new BioJavaStructureDecoder();
     DecodeStructure ds = new DecodeStructure(inputByteArray);
      ds.getStructFromByteArray(biojavaStructureDecoder, parsingParams);
     // Now return this structure
     return biojavaStructureDecoder.getStructure();
   }
}
