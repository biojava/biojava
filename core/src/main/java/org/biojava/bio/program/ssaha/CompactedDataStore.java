package org.biojava.bio.program.ssaha;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.nio.IntBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

import org.biojava.bio.BioError;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Packing;
import org.biojava.bio.symbol.SymbolList;

/**
 * An implementation of DataStore that will map onto a file using the NIO
 * constructs. You should obtain one of these by using the methods in
 * MappedDataStoreFactory.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */

public class CompactedDataStore implements DataStore {
  private final Packing packing;
  private final int wordLength;
  private final IntBuffer hashTable;
  private final MappedByteBuffer hitTable;
  private final IntBuffer nameArray;
  private final MappedByteBuffer nameTable;
  private final int numSequences;
  
  CompactedDataStore(File dataStoreFile)
      throws IOException 
  {
    FileChannel channel = new FileInputStream(dataStoreFile).getChannel();
    
    MappedByteBuffer rootBuffer = channel.map(
      FileChannel.MapMode.READ_ONLY,
      0,
      4 * 6
    );
    rootBuffer.position(0);
    
    final int hashTablePos = rootBuffer.getInt();
    final int hitTablePos = rootBuffer.getInt();
    final int nameArrayPos = rootBuffer.getInt();
    final int nameTablePos = rootBuffer.getInt();
    wordLength = rootBuffer.getInt();
    
    // extend root map to include the serialized packing
    int packingStreamLength = rootBuffer.getInt();
    //System.out.println("hashTablePos:\t" + hashTablePos);
    //System.out.println("hitTablePos:\t" + hitTablePos);
    //System.out.println("nameArrayPos:\t" + nameArrayPos);
    //System.out.println("nameTablePos:\t" + nameTablePos);
    //System.out.println("packingStreamLength:\t" + packingStreamLength);
    rootBuffer = channel.map(
      FileChannel.MapMode.READ_ONLY,
      0,
      4 * 6 + packingStreamLength
    );
    rootBuffer.position(4 * 6);
    byte[] packingBuffer = new byte[packingStreamLength];
    rootBuffer.get(packingBuffer);
    ByteArrayInputStream packingStream = new ByteArrayInputStream(packingBuffer);
    ObjectInputStream packingSerializer = new ObjectInputStream(packingStream);
    
    try {
      this.packing = (Packing) packingSerializer.readObject();
    } catch (ClassNotFoundException cnfe) {
      throw new Error("Can't restore packing", cnfe);
    }
    
    // map in buffer for the hash table
    MappedByteBuffer hashTable_MB = channel.map(
      FileChannel.MapMode.READ_ONLY,
      hashTablePos,
      4
    );
    hashTable_MB.position(0);
    int hashTableSize = hashTable_MB.getInt();
    hashTable = channel.map(
      FileChannel.MapMode.READ_ONLY,
      hashTablePos + 4,
      hashTableSize - 4
    ).asIntBuffer();
    
    // map in buffer for hit table
    MappedByteBuffer hitTable_MB = channel.map(
      FileChannel.MapMode.READ_ONLY,
      hitTablePos,
      4
    );
    hitTable_MB.position(0);
    int hitTableSize = hitTable_MB.getInt();
    hitTable = channel.map(
      FileChannel.MapMode.READ_ONLY,
      hitTablePos + 4,
      hitTableSize - 4
    );
    
    // map in buffer for names array
    MappedByteBuffer nameArray_MB = channel.map(
      FileChannel.MapMode.READ_ONLY,
      nameArrayPos,
      4
    );
    nameArray_MB.position(0);
    int nameArraySize = nameArray_MB.getInt();
    numSequences = nameArraySize / 8;
    // System.err.println("numSequences: " + numSequences);
    nameArray = channel.map(
      FileChannel.MapMode.READ_ONLY,
      nameArrayPos + 4,
      nameArraySize - 4
    ).asIntBuffer();
    
    // map in buffer for names table
    MappedByteBuffer nameTable_MB = channel.map(
      FileChannel.MapMode.READ_ONLY,
      nameTablePos,
      4
    );
    nameTable_MB.position(0);
    int nameTableSize = nameTable_MB.getInt();
    nameTable = channel.map(
      FileChannel.MapMode.READ_ONLY,
      nameTablePos + 4,
      nameTableSize - 4
    );
  }
  
  public FiniteAlphabet getAlphabet() {
    return packing.getAlphabet();
  }
  
  public void search(
    String seqID,
    SymbolList symList,
    SearchListener listener
  ) {
    try {
	int word = 0;
	int lengthFromUnknown = 0;
	listener.startSearch(seqID);
	for(int pos = 1; pos <= symList.length(); pos++) {
	    word = word >> (int) packing.wordSize();
	    int p = packing.pack(symList.symbolAt(pos));
	    if (p < 0) {
		lengthFromUnknown = 0;
	    } else {
		lengthFromUnknown++;
		word |= (int) p << ((int) (wordLength - 1) * packing.wordSize());
	    }
	    
	    if (lengthFromUnknown >= wordLength) {
		fireHits(word, pos - wordLength + 1, listener);
	    }
	}
	listener.endSearch(seqID);
    } catch (IllegalSymbolException ise) {
      throw new BioError("Assertion Failure: Symbol dissapeared");
    }
  }
  
  public String seqNameForID(int id) {
    int offset = nameArray.get(id);
    nameTable.position(offset);
    int length = nameTable.getInt();
    StringBuffer sbuff = new StringBuffer(length);
    for(int i = 0; i < length; i++) {
      sbuff.append(nameTable.getChar());
    }
    return sbuff.toString();
  }
  
    private int seqIDForPos(int pos) {
	if (numSequences == 1) {
	    return 0;
	} else {
	    int maxBound = numSequences - 1;
	    int minBound = 0;

	    while (true) {
		int mid = (minBound + maxBound) / 2;
		int offset = nameArray.get((mid * 2) + 1);
		int endOffset = Integer.MAX_VALUE;
		if (mid < (numSequences - 1)) {
		    endOffset = nameArray.get((mid * 2) + 3);
		}
		if (pos > offset && pos < endOffset) {
		    return mid * 2;
		} else if (pos < offset) {
		    maxBound = mid - 1;
		} else if (pos > endOffset) {
		    minBound = mid + 1;
		} else {
		    throw new Error("Ooops: could not locate seq name for " +
                                    "\tpos: " + pos +
                                    "\tmid: " + mid +
                                    "\toffset: " + offset +
                                    "\tendOffset: " + endOffset +
                                    "\tminBound: " + minBound +
                                    "\tmaxBound: " + maxBound);
		}
	    }
	}
    }

    private int offsetForID(int id) {
	if (numSequences == 1) {
	    return 0;
	} else {
	    return nameArray.get(id + 1);
	}
    }

  public void fireHits(
    int word,
    int offset,
    SearchListener listener
  ) {
    int hitOffset = hashTable.get(word);
    if(hitOffset >= 0) {
      try {
        hitTable.position(hitOffset);
      } catch (IllegalArgumentException e) {
        System.out.println("word:\t" + word);
        System.out.println("offset:\t" + offset);
        System.out.println("hitOffset\t" + hitOffset);
        throw e;
        
      }
      int hits = hitTable.getInt();
      
      for(int i = 0; i < hits; i++) {
	  int pos = hitTable.getInt();
	  int id = seqIDForPos(pos);

	  listener.hit(
		       id,
		       offset,
		       pos - offsetForID(id),
		       wordLength
		      );
      }
    } else if (hitOffset == -2) {
	System.err.println("Hit an elided word!");
    }
  }
}
