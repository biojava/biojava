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

/**
 * @author Matthew Pocock
 */
class SecondaryFileAsList
extends SearchableFileAsList {
  private Comparator KEY_VALUE_COMPARATOR = new Comparator() {
    public int compare(Object a, Object b) {
      String as = a.toString();
      String bs = b.toString();

      return BioStore.STRING_CASE_SENSITIVE_ORDER.compare(as, bs);
    }
  };

  public SecondaryFileAsList(File file, int recordLen)
  throws IOException {
    super(file, recordLen);
  }

  public SecondaryFileAsList(File file, boolean mutable)
  throws IOException {
    super(file, mutable);
  }

  protected Object parseRecord(byte[] buffer) {
    int tabI = 0;
    while(buffer[tabI] != '\t') {
      tabI++;
    }
    String prim = new String(buffer, 0, tabI);
    tabI++;
    String sec = new String(buffer, tabI, buffer.length - tabI).trim();

    return new KeyPair.Impl(prim, sec);
  }

  protected void generateRecord(byte[] buffer, Object item) {
    KeyPair kp = (KeyPair) item;

    int i = 0;
    byte[] str = null;

    try {
      str = kp.getPrimary().getBytes();
      for(int j = 0; j < str.length; j++) {
        buffer[i++] = str[j];
      }
    } catch (ArrayIndexOutOfBoundsException e) {
      throw new ArrayIndexOutOfBoundsException(
        "Over ran buffer with primary ID: " + new String(str) +
        " "+ str.length +
        " : " + buffer.length +" index: "+i
      );
    }

    buffer[i++] = '\t';

    try {
      str = kp.getSecondary().getBytes();
      for(int j = 0; j < str.length; j++) {
        if(i < buffer.length)
          buffer[i++] = str[j];
        else
          break;
      }
    } catch (ArrayIndexOutOfBoundsException e) {
      throw new ArrayIndexOutOfBoundsException(
        "Over ran buffer with secondary ID: " + new String(str) +
        " "+ str.length +
        " : " + buffer.length +
        " index: "+i
      );
    }

    while(i < buffer.length) {
      buffer[i++] = ' ';
    }
  }

  public Comparator getComparator() {
    return KEY_VALUE_COMPARATOR;
  }
}
