package org.biojava.bio.program.ssaha;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.nio.LongBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

import org.biojava.bio.BioError;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Packing;
import org.biojava.bio.symbol.PackingFactory;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.Constants;
import org.biojava.utils.io.LargeBuffer;

/**
 * An implementation of DataStore that will map onto a file using the NIO
 * constructs. You should obtain one of these by using the methods in
 * NIODataStoreFactory.
 *
 * @author Matthew Pocock
 */
class NIODataStore implements DataStore {
  private final Packing packing;
  private final int wordLength;
  private final LongBuffer hashTable;
  private final LargeBuffer hitTable;
  private final MappedByteBuffer nameTable;
  
  NIODataStore(File dataStoreFile)
  throws IOException {
    FileChannel channel = new FileInputStream(dataStoreFile).getChannel();
    
    MappedByteBuffer rootBuffer = channel.map(
      FileChannel.MapMode.READ_ONLY,
      0,
      3 * Constants.BYTES_IN_LONG + // positions
      2 * Constants.BYTES_IN_INT    // word length & packing length
    );
    rootBuffer.position(0);

    final long hashTablePos = rootBuffer.getLong();
    final long hitTablePos = rootBuffer.getLong();
    final long nameTablePos = rootBuffer.getLong();
    wordLength = rootBuffer.getInt();
    
    // extend root map to include the serialized packing
    int packingStreamLength = rootBuffer.getInt();
    System.out.println("hashTablePos:\t" + hashTablePos);
    System.out.println("hitTablePos:\t" + hitTablePos);
    System.out.println("nameTablePos:\t" + nameTablePos);
    System.out.println("packingStreamLength:\t" + packingStreamLength);
    rootBuffer = channel.map(
      FileChannel.MapMode.READ_ONLY,
      0,
      3 * Constants.BYTES_IN_LONG + // positions
      2 * Constants.BYTES_IN_INT +  // word length & packing length
      packingStreamLength           // serialized packing
    );
    rootBuffer.position(3 * Constants.BYTES_IN_LONG + 2 * Constants.BYTES_IN_INT);
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
      Constants.BYTES_IN_LONG
    );
    hashTable_MB.position(0);
    long hashTableSize = hashTable_MB.getLong();
    System.out.println("hashTableSize: " + hashTableSize);
    hashTable = channel.map(
      FileChannel.MapMode.READ_ONLY,
      hashTablePos + Constants.BYTES_IN_LONG,
      hashTableSize - Constants.BYTES_IN_LONG
    ).asLongBuffer();
    
    // map in buffer for names table
    MappedByteBuffer nameTable_MB = channel.map(
      FileChannel.MapMode.READ_ONLY,
      nameTablePos,
      Constants.BYTES_IN_INT
    );
    nameTable_MB.position(0);
    int nameTableSize = nameTable_MB.getInt();
    nameTable = channel.map(
      FileChannel.MapMode.READ_ONLY,
      nameTablePos + Constants.BYTES_IN_INT,
      nameTableSize - Constants.BYTES_IN_INT
    );

    // map in buffer for hit table
    MappedByteBuffer hitTable_MB = channel.map(
      FileChannel.MapMode.READ_ONLY,
      hitTablePos,
      Constants.BYTES_IN_LONG
    );
    hitTable_MB.position(0);
    long hitTableSize = hitTable_MB.getLong();
    hitTable = new LargeBuffer(
      channel,
      FileChannel.MapMode.READ_ONLY,
      hitTablePos + Constants.BYTES_IN_LONG,
      hitTableSize - Constants.BYTES_IN_LONG
    );
  }
  
  public FiniteAlphabet getAlphabet() {
    return packing.getAlphabet();
  }
  
  public void search(
    String seqID,
    SymbolList symList,
    SearchListener listener
  ) throws SearchException {
    try {
      int word = PackingFactory.primeWord(symList, wordLength, packing);
      listener.startSearch(seqID);
      fireHits(word, 1, listener);
      for(int j = wordLength + 2; j <= symList.length(); j++) {
        word = PackingFactory.nextWord(symList, word, j, wordLength, packing);
        fireHits(word, j - wordLength, listener);
      }
      listener.endSearch(seqID);
    } catch (IllegalSymbolException ise) {
      throw new BioError("Assertion Failure: Symbol dissapeared");
    } catch (IOException ioe) {
      throw new SearchException(ioe);
    }
  }
  
  public String seqNameForID(int id)
  throws SearchException {;
//    try {
      nameTable.position(id);
      int length = nameTable.getShort();
      StringBuffer sbuff = new StringBuffer(length);
      for(int i = 0; i < length; i++) {
        sbuff.append(nameTable.getChar());
      }
      return sbuff.toString();
//    } catch (IOException ioe) {
//      throw new NestedException(ioe);
//    }
  }
  
  private void fireHits(
    int word,
    int offset,
    SearchListener listener
  ) throws IOException {
    long hitOffset = hashTable.get(word);
    //System.out.println("hitOffset: " + hitOffset);
    if(hitOffset != -1) {
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
        listener.hit(
          hitTable.getInt(),
          offset,
          hitTable.getInt(),
          wordLength
        );
      }
    }
  }
}

