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
package org.biojava.utils.io;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.RandomAccessFile;

/**
 * @author Matthew Pocock
 */
public class RAF extends RandomAccessFile {
  private File file;
  
  public RAF(File file, String mode)
  throws FileNotFoundException {
    super(file, mode);
    this.file = file;
  }
  
  public RAF(String name, String mode)
  throws FileNotFoundException {
    super(name, mode);
    this.file = new File(name);
  }
  
  public File getFile() {
    return file;
  }
}
