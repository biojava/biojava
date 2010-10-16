package org.biojava.bio.program.ssaha;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.RandomAccessFile;
import java.nio.BufferOverflowException;
import java.nio.LongBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.Packing;
import org.biojava.bio.symbol.PackingFactory;
import org.biojava.utils.Constants;
import org.biojava.utils.io.LargeBuffer;

/**
 * <p>
 * Builder for a datastore that has no practical file size limit.
 * </p>
 *
 * <p>
 * This implementation of the data store factory uses longs as indecies
 * internaly, so can be used with files exceeding 2 gigs in size.
 * </p>
 *
 * <p>
 * The data store file has the following structure.
 * <pre>
 * file: header, hash table, nameTable, hitTable
 *
 * header:
 *   long hashTablePos, // byte offset in file
 *   long hitTablePos,  // byte offset in file
 *   long nameTablePos, // byte offset in file
 *   int wordLength,
 *   int serializedPackingLength,
 *   byte[] serializedPacking
 *
 * hashTable:
 *   long hashTableLength,
 *   long[hashTableLength] hits // byte offset into hitTable
 *
 * nameTable:
 *   int nameTableSize, // size in bytes
 *   (short nameLength, char[nameLength] name)[nameTableSize] names
 *
 * hitTable:
 *   long hitTableSize, // size in bytes
 *   hitTableRecord[hitTableSize] hits
 *
 * hitTableRecord:
 *   int hitCount,
 *   hitRecord[hitCount] hit
 *
 * hit:
 *   long seqOffset, // byte offset into sequence names table
 *   int pos         // biological position in sequence
 * </pre>
 * </p>
 * @author Matthew Pocock
 */
public class NIODataStoreFactory
implements DataStoreFactory {
  public DataStore getDataStore(File storeFile)
  throws IOException {
    return new NIODataStore(storeFile);
  }
  
  public DataStore buildDataStore(
    File storeFile,
    SequenceDB seqDB,
    Packing packing,
    int wordLength,
    int threshold
  ) throws
    IllegalAlphabetException,
    IOException,
    BioException
  {
    ByteArrayOutputStream packingStream = new ByteArrayOutputStream();
    ObjectOutputStream packingSerializer = new ObjectOutputStream(packingStream);
    packingSerializer.writeObject(packing);
    packingSerializer.flush();
    
    final int structDataSize =
      3 * Constants.BYTES_IN_LONG + // positions
      2 * Constants.BYTES_IN_INT +  // word length & packing length
      packingStream.toByteArray().length;
    
    final long hashTablePos;
    final long hitTablePos;
    final long nameTablePos;
    
    storeFile.createNewFile();
    final RandomAccessFile store = new RandomAccessFile(storeFile, "rw");
    final FileChannel channel = store.getChannel();
    
    // allocate array for k-tuple -> hit list
    System.out.println("Word length:\t" + wordLength);
    int words = 1 << (
      (int) packing.wordSize() *
      (int) wordLength
    );
    System.out.println("Words:\t" + words);
    
    hashTablePos = structDataSize;
    long hashTableSize =
      Constants.BYTES_IN_LONG +        // hash table length
      words * Constants.BYTES_IN_LONG; // hash table entries
    
    System.out.println("Allocated:\t" + hashTableSize);
    final MappedByteBuffer hashTable_MB = channel.map(
      FileChannel.MapMode.READ_WRITE,
      hashTablePos,
      hashTableSize
    );
    final LongBuffer hashTable = hashTable_MB.asLongBuffer();
    hashTable.put(0, hashTableSize); // write length of k-tuple array
    
    // initialize counts to zero
    for(int i = 0; i < words; i++) {
      hashTable.put(i+1, 0);
    }
    hashTable.position(0);
    System.out.println("Initialized hash table");
    
    // 1st pass
    // writes counts as longs for each k-tuple
    // count up the space required for sequence names
    //
    System.out.println("First parse");
    int seqCount = 0;
    int nameChars = 0;
    for(SequenceIterator i = seqDB.sequenceIterator(); i.hasNext(); ) {
      Sequence seq = i.nextSequence();
      System.out.println(seq.getName());
      if(seq.length() > wordLength) {
        seqCount++;
        nameChars += seq.getName().length();
        
        int word = PackingFactory.primeWord(seq, wordLength, packing);
        //PackingFactory.binary(word);
        addCount(hashTable, word);
        for(int j = wordLength+2; j <= seq.length(); j++) {
          word = PackingFactory.nextWord(seq, word, j, wordLength, packing);
          //PackingFactory.binary(word);
          addCount(hashTable, word);
        }
      }
    }
    System.out.println();
    System.out.println("Done");
     
    // map the space for sequence names as short length, char* name
    //
    nameTablePos = hashTablePos + hashTableSize;
    int nameTableSize =
      Constants.BYTES_IN_INT +              // bytes in table
      seqCount * Constants.BYTES_IN_SHORT + // string lengths
      nameChars * Constants.BYTES_IN_CHAR;  // characters
    System.out.println("nameTableSize:\t" + nameTableSize);
    final MappedByteBuffer nameTable = channel.map(
      FileChannel.MapMode.READ_WRITE,
      nameTablePos,
      nameTableSize
    );
    nameTable.putInt(0, nameTableSize);
    nameTable.position(Constants.BYTES_IN_INT);
    
    // add up the number of k-tuples
    //
    long kmersUsed = 0; // number of kmers with valid hits
    long hitCount = 0;  // total number of individual hits
    for(int i = 0; i < words; i++) {
      long counts = hashTable.get(i + 1);
      if(counts > 0 && counts < threshold) {
        hitCount += counts;
        kmersUsed++;
      }
    }
    System.out.println("hitCount: " + hitCount);
    System.out.println("kmersUsed: " + kmersUsed);
    
    hitTablePos = nameTablePos + nameTableSize;
    long hitTableSize =
      (long) Constants.BYTES_IN_LONG +                            // size
      (long) hitCount * (Constants.BYTES_IN_INT + Constants.BYTES_IN_INT) +  // list elements
      (long) kmersUsed * Constants.BYTES_IN_INT;                  // size of lists
    System.out.println("hitTableSize:\t" + hitTableSize);
    System.out.println("hitTablePos:\t" + hitTablePos);
    
    final LargeBuffer hitTable = new LargeBuffer(
      channel,
      FileChannel.MapMode.READ_WRITE,
      hitTablePos,
      hitTableSize
    );
    hitTable.position(0);
    hitTable.putLong(hitTableSize);
    
    // write locations of hit arrays
    System.out.println("Writing pointers from hash to hits");
    long hitOffset = 0;
    for(int i = 0; i < words; i++) {
      long counts = hashTable.get(i+1);
      if(counts > 0 && counts < threshold) {
        try {
          // record location of a block of the form:
          // n,(seqID,offset)1,(seqID,offset)2,...,(seqID,offset)n
          if(hitOffset < 0) {
            throw new IndexOutOfBoundsException("Hit offset negative");
          }
          hashTable.put(i + 1, hitOffset); // wire hash table to hit table
          hitTable.putInt(hitOffset + Constants.BYTES_IN_LONG, 0); // initialy we have no hits
          hitOffset +=
            Constants.BYTES_IN_INT +
            counts * (Constants.BYTES_IN_INT + Constants.BYTES_IN_INT);
        } catch (IndexOutOfBoundsException e) {
          System.out.println("counts:\t" + counts);
          System.out.println("word:\t" + i);
          System.out.println("hitOffset:\t" + hitOffset);
          throw e;
        }
      } else {
        // too many hits - set the number of hits to the flag value -1
        hashTable.put(i + 1, -1);
      }
    }
    
    // 2nd parse
    // write sequence array and names
    // write hitTable
    System.out.println("2nd parse");
    int seqNumber = 0;
    nameTable.position(Constants.BYTES_IN_INT);
    for(SequenceIterator i = seqDB.sequenceIterator(); i.hasNext(); ) {
      Sequence seq = i.nextSequence();
      System.out.println(seq.getName());
      seqNumber++;
      if(seq.length() > wordLength) {
        try {
          // write sequence name length and chars into nameTable
          int namePos = nameTable.position() - Constants.BYTES_IN_INT;
          String name = seq.getName();
          nameTable.putShort((short) name.length());
          for(int j = 0; j < name.length(); j++) {
            nameTable.putChar((char) name.charAt(j));
          }
          
          // write k-mer seq,offset
          int word = PackingFactory.primeWord(seq, wordLength, packing);
          writeRecord(hashTable, hitTable, 1, namePos, word);
          for(int j = wordLength+2; j <= seq.length(); j++) {
            word = PackingFactory.nextWord(seq, word, j, wordLength, packing);
            writeRecord(hashTable, hitTable, j - wordLength, namePos, word);
          }
        } catch (BufferOverflowException e) {
          System.out.println("name:\t" + seq.getName());
          System.out.println("seqNumber:\t" + seqNumber);
          System.out.println("nt pos:\t" + nameTable.position());
          throw e;
        }
      }
    }
    System.out.println();
    System.out.println("Done");

    final MappedByteBuffer rootBuffer = channel.map(
      FileChannel.MapMode.READ_WRITE,
      0,
      structDataSize
    );
    
    rootBuffer.position(0);
    rootBuffer.putLong(hashTablePos);
    rootBuffer.putLong(hitTablePos);
    rootBuffer.putLong(nameTablePos);
    rootBuffer.putInt(wordLength);
    rootBuffer.putInt(packingStream.toByteArray().length);
    rootBuffer.put(packingStream.toByteArray());
    
    rootBuffer.force();
    hashTable_MB.force();
    hitTable.force();
    nameTable.force();
    
    System.out.println("Committed");
    
    return getDataStore(storeFile);
  }
  
  private void addCount(LongBuffer buffer, int word) {
    long count = buffer.get(word+1);
    count++;
    buffer.put(word+1, count);
  }
  
  private void writeRecord(
    LongBuffer hashTable,
    LargeBuffer hitTable,
    int offset,
    int seqNumber,
    int word
  ) throws IOException {
    long kmerPointer = hashTable.get(word+1);
    if(kmerPointer != -1) {
      kmerPointer += Constants.BYTES_IN_LONG;

      int hitCount = hitTable.getInt(kmerPointer);
      long pos = kmerPointer + hitCount * (Constants.BYTES_IN_INT + Constants.BYTES_IN_INT) + Constants.BYTES_IN_INT;
      
      hitTable.position(pos);
      hitTable.putInt(seqNumber);
      hitTable.putInt(offset);
      hitTable.putInt(kmerPointer, hitCount + 1);
    }
  }
  
}
