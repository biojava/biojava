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

package org.biojava.bio.program.indexdb;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;

import org.biojava.utils.io.RAF;

/**
 * @author Matthew Pocock
 * @author Keith James
 */
class PrimaryIDList
extends SearchableFileAsList {
  private Comparator INDEX_COMPARATOR = new Comparator() {
    public int compare(Object a, Object b) {
      String as;
      String bs;
      
      if(a instanceof Record) {
        as = ((Record) a).getID();
      } else {
        as = (String) a;
      }
      
      if(b instanceof Record) {
        bs = ((Record) b).getID();
      } else {
        bs = (String) b;
      }
      
      return BioStore.STRING_CASE_SENSITIVE_ORDER.compare(as, bs);
    }
  };
  
  private BioStore store;
  
  public PrimaryIDList(File file, int recordLen, BioStore store)
  throws IOException {
    super(file, recordLen);
    this.store = store;
  }
  
  public PrimaryIDList(File file, BioStore store, boolean mutable)
  throws IOException {
    super(file, mutable);
    this.store = store;
  }
  
  protected Object parseRecord(byte[] buffer) {
    int lastI = 0;
    int newI = 0;
    while(buffer[newI] != '\t') {
      newI++;
    }
    String id = new String(buffer, lastI, newI - lastI);
    
    lastI = ++newI;
    while(buffer[newI] != '\t') {
      newI++;
    }
    RAF file = store.getFileForID(Integer.parseInt(new String(buffer, lastI, newI - lastI).trim()));
    
    lastI = ++newI;
    while(buffer[newI] != '\t') {
      newI++;
    }
    long start = Long.parseLong(new String(buffer, lastI, newI - lastI));
    
    newI++;
    int length = Integer.parseInt(
      new String(buffer, newI, buffer.length - newI).trim()
    );
    
    return new Record.Impl(id, file, start, length);
  }
  
  protected void generateRecord(byte[] buffer, Object item)
  throws IOException {
    String id = null;
    int fileID = -1;
    String start = null;
    String length = null;

    try {
      Record indx = (Record) item;
      
      id = indx.getID();
      if(id == null) {
        throw new NullPointerException("Can't process null ID: " + indx);
      }
      fileID = store.getIDForFile(indx.getFile());
      start = String.valueOf(indx.getOffset());
      length = String.valueOf(indx.getLength());

      int i = 0;
      byte[] str;
      
      str = id.getBytes();
      for(int j = 0; j < str.length; j++) {
        buffer[i++] = str[j];
      }
      
      buffer[i++] = '\t';
      
      str = String.valueOf(fileID).getBytes();
      for(int j = 0; j < str.length; j++) {
        buffer[i++] = str[j];
      }
      
      buffer[i++] = '\t';
      
      str = start.getBytes();
      for(int j = 0; j < str.length; j++) {
        buffer[i++] = str[j];
      }
      
      buffer[i++] = '\t';
      
      str = length.getBytes();
      for(int j = 0; j < str.length; j++) {
        buffer[i++] = str[j];
      }
      
      while(i < buffer.length) {
        buffer[i++] = ' ';
      }
    } catch (IOException ex) {
      String attemptedLine =
        id + "\t" +
        fileID + "\t" +
        start + "\t" +
        length;
      throw (IOException) new IOException(
        "Could not build record. Record length: " + buffer.length +
        " Line length: " + attemptedLine.length() +
        " " + attemptedLine).initCause(ex);
    } catch (ArrayIndexOutOfBoundsException ex) {
      String attemptedLine =
        id + "\t" +
        fileID + "\t" +
        start + "\t" +
        length;
      throw (IOException) new IOException(
        "Could not build record. Record length: " + buffer.length +
        " Line length: " + attemptedLine.length() +
        " " + attemptedLine).initCause(ex);
    }
  }
  
  public Comparator getComparator() {
    return INDEX_COMPARATOR;
  }
}

