package org.biojava.bio.program.ssaha;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.RandomAccessFile;
import java.nio.BufferOverflowException;
import java.nio.IntBuffer;
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

/**
 * <p>
 * Builder for a data store that is backed by a java.nio.MappedByteBuffer.
 * This has a limitation that the total size of the mapped buffer and
 * therefore the hash table can not exceed 2 gigs.
 * </p>
 *
 * <p>
 * The data store file has the following structure.
 * <pre>
 * file: header, hash table, nameArray, nameTable, hitTable
 *
 * header:
 *   int hashTablePos, // byte offset in file
 *   int hitTablePos,  // byte offset in file
 *   int nameArrayPos, // byte offset in file
 *   int nameTablePos, // byte offset in file
 *   int wordLength,
 *   int serializedPackingLength,
 *   byte[] serializedPacking
 *
 *   hash table:
 *     int hashTableLength,
 *     int[hashTableLength] hits // index into hitTable
 *
 *  nameArray:
 *    int nameArrayLength,
 *    int[nameArrayLength] nameArray // byte offset into nameTable
 *
 *  nameTable:
 *    int nameTableSize, // size in bytes
 *    (short nameLength, char[nameLength] name)[nameTableSize] names
 *
 *  hitTable:
 *    int hitTableSize, // size in bytes
 *    hitTableRecord[hitTableSize] hits
 *
 *  hitTableRecord:
 *    int hitCount,
 *    hitRecord[hitCount] hit
 *
 *  hit:
 *    int seqIndex, // index into nameArray
 *    int offset    // offset into the sequence
 * </pre>
 * </p>
 *
 * @author Matthew Pocock
 */
public class MappedDataStoreFactory
implements DataStoreFactory {
  public DataStore getDataStore(File storeFile)
  throws IOException {
    return new MappedDataStore(storeFile);
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
      6 * Constants.BYTES_IN_INT +
      packingStream.toByteArray().length;
    
    final int hashTablePos;
    final int hitTablePos;
    final int nameArrayPos;
    final int nameTablePos;
    
    storeFile.createNewFile();
    final RandomAccessFile store = new RandomAccessFile(storeFile, "rw");
    final FileChannel channel = store.getChannel();
    
    // allocate array for k-tuple -> hit list
    //System.out.println("Word length:\t" + wordLength);
    int words = 1 << (
      (int) packing.wordSize() *
      (int) wordLength
    );
    //System.out.println("Words:\t" + words);
    
    hashTablePos = structDataSize;
    int hashTableSize =
      (int) Constants.BYTES_IN_INT + // hash table length
      words * (int) Constants.BYTES_IN_INT; // hash table entries
    
    //System.out.println("Allocated:\t" + hashTableSize);
    final MappedByteBuffer hashTable_MB = channel.map(
      FileChannel.MapMode.READ_WRITE,
      hashTablePos,
      hashTableSize
    );
    final IntBuffer hashTable = hashTable_MB.asIntBuffer();
    hashTable.put(0, hashTableSize); // write length of k-tuple array
    
    // initialize counts to zero
    for(int i = 0; i < words; i++) {
      hashTable.put(i+1, 0);
    }
    hashTable.position(0);
    
    // 1st pass
    // writes counts as ints for each k-tuple
    // count up the space required for sequence names
    //
    int seqCount = 0;
    int nameChars = 0;
    for(SequenceIterator i = seqDB.sequenceIterator(); i.hasNext(); ) {
      Sequence seq = i.nextSequence();
      if(seq.length() > wordLength) {
        seqCount++;
        nameChars += seq.getName().length();
        
        int word = PackingFactory.primeWord(seq, wordLength, packing);
        //PackingFactory.binary(word);
        addCount(hashTable, word);
        for(int j = wordLength + 1; j <= seq.length(); j++) {
          word = PackingFactory.nextWord(seq, word, j, wordLength, packing);
          //PackingFactory.binary(word);
          addCount(hashTable, word);
        }
      }
    }
    
    // map the space for sequence index->name
    //
    nameArrayPos = hashTablePos + hashTableSize;
    int nameArraySize = (seqCount + 1) * Constants.BYTES_IN_INT;
    //System.out.println("seqCount:\t" + seqCount);
    //System.out.println("nameArraySize:\t" + nameArraySize);
    final MappedByteBuffer nameArray_MB = channel.map(
      FileChannel.MapMode.READ_WRITE,
      nameArrayPos,
      nameArraySize
    );
    final IntBuffer nameArray = nameArray_MB.asIntBuffer();
    nameArray.put(0, nameArraySize);
    
    // map the space for sequence names as short length, char* name
    //
    nameTablePos = nameArrayPos + nameArraySize;
    int nameTableSize =
      Constants.BYTES_IN_INT +
      seqCount * Constants.BYTES_IN_INT +
      nameChars * Constants.BYTES_IN_CHAR;
    //System.out.println("nameTableSize:\t" + nameTableSize);
    final MappedByteBuffer nameTable = channel.map(
      FileChannel.MapMode.READ_WRITE,
      nameTablePos,
      nameTableSize
    );
    nameTable.putInt(0, nameTableSize);
    nameTable.position(Constants.BYTES_IN_INT);
    
    // add up the number of k-tuples
    //
    int kmersUsed = 0;
    int hitCount = 0;
    for(int i = 0; i < words; i++) {
      int counts = hashTable.get(i + 1);
      if(counts > 0 && counts < threshold) {
        hitCount++;
        kmersUsed += counts;
      }
    }
    
    // map the space for hits
    hitTablePos = nameTablePos + nameTableSize;
    long hitTableSize =
      (long) Constants.BYTES_IN_INT +                            // size
      (long) kmersUsed * (Constants.BYTES_IN_INT + Constants.BYTES_IN_INT) +  // list elements
      (long) hitCount * Constants.BYTES_IN_INT;                  // size of lists
    //System.out.println("hitTableSize:\t" + hitTableSize);
    //System.out.println("hitTableSize:\t" + (int) hitTableSize);
    //System.out.println("hitTablePos:\t" + hitTablePos);
    final MappedByteBuffer hitTable = channel.map(
      FileChannel.MapMode.READ_WRITE,
      hitTablePos,
      (int) hitTableSize
    );
    hitTable.putInt(0, (int) hitTableSize);
    hitTable.position(Constants.BYTES_IN_INT);
    
    // write locations of hit arrays
    int hitOffset = 0;
    for(int i = 0; i < words; i++) {
      int counts = hashTable.get(i+1);
      if(counts > 0 && counts < threshold) {
        try {
        // record location of a block of the form:
        // n,(seqID,offset)1,(seqID,offset)2,...,(seqID,offset)n
        if(hitOffset < 0) {
          throw new IndexOutOfBoundsException("Hit offset negative");
        }
        hashTable.put(i + 1, hitOffset); // wire hash table to hit table
        hitTable.putInt(hitOffset + Constants.BYTES_IN_INT, 0); // initialy we have no hits
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
    int seqNumber = 0;
    nameTable.position(Constants.BYTES_IN_INT);
    for(SequenceIterator i = seqDB.sequenceIterator(); i.hasNext(); ) {
      Sequence seq = i.nextSequence();
      
      if(seq.length() > wordLength) {
        try {
          
          // write sequence name reference into nameArray
          nameArray.put(seqNumber + 1, nameTable.position()-Constants.BYTES_IN_INT);
          
          // write sequence name length and chars into nameTable
          String name = seq.getName();
          nameTable.putInt(name.length());
          for(int j = 0; j < name.length(); j++) {
            nameTable.putChar((char) name.charAt(j));
          }
          
          // write k-mer seq,offset
          int word = PackingFactory.primeWord(seq, wordLength, packing);
          writeRecord(hashTable, hitTable, 1, seqNumber, word);
          for(int j = wordLength+1; j <= seq.length(); j++) {
            word = PackingFactory.nextWord(seq, word, j, wordLength, packing);
            writeRecord(hashTable, hitTable, j - wordLength + 1, seqNumber, word);
          }
        } catch (BufferOverflowException e) {
          System.out.println("name:\t" + seq.getName());
          System.out.println("seqNumber:\t" + seqNumber);
          System.out.println("na pos:\t" + nameArray.position());
          System.out.println("nt pos:\t" + nameTable.position());
          throw e;
        }
        seqNumber++;
      }
    }
    
    //validateNames(seqCount, nameArray, nameTable);
    
    final MappedByteBuffer rootBuffer = channel.map(
      FileChannel.MapMode.READ_WRITE,
      0,
      structDataSize
    );
    
    rootBuffer.position(0);
    rootBuffer.putInt(hashTablePos);
    rootBuffer.putInt(hitTablePos);
    rootBuffer.putInt(nameArrayPos);
    rootBuffer.putInt(nameTablePos);
    rootBuffer.putInt(wordLength);
    rootBuffer.putInt(packingStream.toByteArray().length);
    rootBuffer.put(packingStream.toByteArray());
    
    rootBuffer.force();
    hashTable_MB.force();
    hitTable.force();
    nameArray_MB.force();
    nameTable.force();
    
    return getDataStore(storeFile);
  }
  
  private void addCount(IntBuffer buffer, int word) {
    int count = buffer.get(word+1);
    count++;
    buffer.put(word+1, count);
  }
  
  private void writeRecord(
    IntBuffer hashTable,
    MappedByteBuffer hitTable,
    int offset,
    int seqNumber,
    int word
  ) {
    int kmerPointer = hashTable.get(word+1);
    if(kmerPointer != -1) {
      kmerPointer += Constants.BYTES_IN_INT;

      int hitCount = hitTable.getInt(kmerPointer);
      int pos = kmerPointer + hitCount * (Constants.BYTES_IN_INT + Constants.BYTES_IN_INT) + Constants.BYTES_IN_INT;
      
      hitTable.position(pos);
      hitTable.putInt(seqNumber);
      hitTable.putInt(offset);
      hitTable.putInt(kmerPointer, hitCount + 1);
    }
  }
  
}
