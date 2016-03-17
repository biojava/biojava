package org.biojava.nbio.structure.io.mmtf;

import java.io.FileNotFoundException;
import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.rcsb.mmtf.decoder.ParsingParams;
import org.rcsb.mmtf.examples.HandleIO;

/**
 * An example script to run the decoder.
 * @author Anthony Bradley
 *
 */
public final class ExampleScript {


  /**
   * Make the constructor private.
   */
  private ExampleScript() {

  }

  /**
   * The main method.
   *
   * @param args the arguments
   * @throws FileNotFoundException the file not found exception
   * @throws IOException Signals that an I/O exception has occurred.
   */
  public static void main(final String[] args) throws
  FileNotFoundException, IOException {
    // Set the PDB_CACHE_DIR path
    System.setProperty("PDB_CACHE_DIR", "/Users/anthony/PDB_CACHE");
    HandleIO gbjs = new HandleIO();
    byte[] inputByteArr = gbjs.getFromUrl("1qmz");
    ParsingParams parsingParms = new ParsingParams();
    parsingParms.setParseInternal(false);
    ParseUsingBioJava parseUseBiojava = new ParseUsingBioJava();
    Structure biojavaStruct = parseUseBiojava.getBiojavaStruct(inputByteArr, parsingParms);
    System.out.println(biojavaStruct.getChains());
  }

}
